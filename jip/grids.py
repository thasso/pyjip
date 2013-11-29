#!/usr/bin/env python
import multiprocessing
import os
import subprocess
import signal
import sys

import jip.cluster
from jip.logger import getLogger


class LocalJob(object):
    def __init__(self, cmd=None, cwd=None, stdout=None, stderr=None,
                 dependencies=None, threads=1, job_id=None):
        self.job_id = job_id
        self.dependencies = set([]) if dependencies is None else dependencies
        self.threads = threads
        self.cmd = cmd
        self.working_directory = cwd
        self.stdout = stdout
        self.stderr = stderr
        self.children = set([])
        self.process = None

    def __lt__(self, other):
        num_deps = len(self.dependencies)
        o_deps = len(other.dependencies)
        if num_deps == o_deps:
            return self.job_id < other.job_id
        return num_deps < o_deps

    def __repr__(self):
        return "{Job-%s[deps: %s][children: %s]}" % (
            str(self.job_id), self.dependencies, self.children
        )


class LocalCluster(jip.cluster.Cluster):

    def __init__(self, _start=True):
        self._current_id = 0
        self.log = getLogger("jip.grids.LocalCluster")
        self.master_requests = multiprocessing.Queue()
        self.master_response = multiprocessing.Queue()
        # start the mater process
        self.master_process = multiprocessing.Process(
            target=LocalCluster.master,
            args=[self.master_requests, self.master_response],
            name="worker-master"
        )
        if _start:
            self.start()

    def start(self):
        self.log.info("Starting local cluster master")
        self.master_process.start()

    def shutdown(self):
        self.log.info("Requesting master shutdown")
        self.master_requests.put(['EXIT'])
        self.master_process.join()
        self.master_requests.close()
        self.master_response.close()

    def wait(self):
        self.log.info("Waiting for termination of the master")
        self.master_requests.put(["WAIT"])
        self.master_process.join()
        self.master_requests.close()
        self.master_response.close()

    def _next_id(self):
        """Increment the current id and return the next available job
        id.

            >>> c = LocalCluster(_start=False)
            >>> assert c._current_id == 0
            >>> assert c._next_id() == 1
            >>> assert c._next_id() == 2

        :returns: next id
        """
        self._current_id += 1
        return self._current_id

    def list(self):
        self.master_requests.put(["JOBS"])
        return self.master_response.get()

    def submit(self, job):
        #assign job id
        job_id = self._next_id()
        job.job_id = job_id

        # set log files
        cwd = job.working_directory if job.working_directory is not None \
            else os.getcwd()
        if job.stderr is None:
            job.stderr = os.path.join(cwd, "jip-%d.err" % (job_id))
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "jip-%d.out" % (job_id))

        #collect dependencies
        deps = set([])
        for dep in [d for d in job.dependencies if d.job_id]:
            deps.add(int(dep.job_id))

        # build the local job
        local_job = LocalJob(
            job.get_cluster_command(),
            cwd,
            job.stdout,
            job.stderr,
            dependencies=deps,
            threads=job.threads
        )

        self.log.info("submitting new job with id: %s", job_id)
        self.master_requests.put([
            "SUBMIT", job_id, local_job,
        ])
        return job

    def cancel(self, job):
        id = job.job_id
        self.log.info("Send cancel request for %s", id)
        self.master_requests.put(["CANCEL", id])

    @staticmethod
    def master(requests, response):
        """The master process function.

        This is started as the main background process and reads messages
        from the master queue.

        :param queue: the master queue
        """
        jobs = {}
        running = {}
        cores = multiprocessing.cpu_count()
        slots = [cores, cores]
        log = getLogger("jip.grids.Master")
        log.info("Master | Available slots: %s", slots[1])

        def remove_children(job):
            # recursive remove of all child jobs
            # that are still queued
            log.debug("Master | Removing all children of %s", job.job_id)
            for c in job.children:
                if c in jobs:
                    remove_children(jobs[c])
                if c in jobs:
                    update_dependencies(jobs[c])
                    del jobs[c]

        def update_dependencies(job):
            for c in job.children:
                if c in jobs:
                    log.debug("Master | Removing %s from child %s "
                              "dependencies", job.job_id, c)
                    jobs[c].dependencies.remove(job.job_id)

        def run_job(job):
            job_id = job.job_id
            log.info("Master | Starting job %s", job_id)
            process = multiprocessing.Process(
                target=LocalCluster.execute,
                args=[requests, job_id, job],
                name="worker-%s" % str(job_id)
            )
            job.process = process
            running[job_id] = job
            # remove the job from the queues jobs
            del jobs[job_id]
            # start the process
            process.start()

        def schedule():
            log.info("Master | running scheduler %d/%d", slots[0], slots[1])
            # find next jobs that can be executed and run it
            #
            # sort all queued jobs
            sorted_jobs = sorted(jobs.values())
            # now iterate the jobs and run the first one that fits
            # into the available slots
            for job in sorted_jobs:
                if len(job.dependencies) > 0:
                    # not all dependencies are resolved
                    log.info("Master | No jobs without "
                             "dependencies found: %s", jobs)
                    break
                if job.threads <= slots[0]:
                    log.info("Master | Submitting job for execution %s",
                             job.job_id)
                    # start this one
                    run_job(job)
                    # update slots
                    slots[0] = slots[0] - job.threads
                    # if we still have slots open,
                    # start another scheduling
                    if slots[0] >= 1:
                        schedule()
                    return
        wait = False

        while True:
            msg = requests.get()
            log.debug("Master | received command: %s", msg)
            command = msg[0]
            if command == "EXIT":
                log.info("Master | EXIT request, shutting down")
                for job_id, p in running.iteritems():
                    log.warn("Master | Terminating %s", job_id)
                    p.process.terminate()
                    p.process.join()
                break
            elif command == "JOBS":
                log.info("Master | list jobs")
                # queued and running
                response.put(list(jobs.keys()) + list(running.keys()))
            elif command == "SUBMIT":
                job_id = int(msg[1])
                job = msg[2]
                job.job_id = job_id
                jobs[job_id] = job
                # update children
                for d in job.dependencies:
                    if d in jobs:
                        jobs[d].children.add(job_id)
                    if d in running:
                        running[d].children.add(job_id)
                log.info("Master | Queue new job %s", job_id)
                schedule()
            elif command == "DONE":
                job_id = int(msg[1])
                state = msg[2]
                if job_id in running:
                    job = running[job_id]
                    job.process.join()
                    job.process = None
                    slots[0] = slots[0] + job.threads
                    update_dependencies(job)
                    if state != 0:
                        log.error("Master | Job %s failed with %s",
                                  job_id, state)
                        remove_children(job)
                    else:
                        log.info("Master | Job %s finished with %s",
                                 job_id, state)
                    job.process = None
                    del running[job_id]
                    schedule()
            elif command == "WAIT":
                if not wait:
                    log.info("Master | Wait mode enabled")
                    wait = True
            elif command == "FAILED":
                job_id = int(msg[1])
                error = msg[2]
                log.error("Master | Execution of %s failed: %s", job_id, error)
                if job_id in running:
                    job = running[job_id]
                    job.process = None
                    update_dependencies(job)
                    remove_children(job)
                    slots[0] = slots[0] + job.threads
                    del running[job_id]
                    schedule()
            elif command == "CANCEL":
                job_id = int(msg[1])
                log.warn("Master | Execution of %s canceled", job_id)
                job = None
                if job_id in running:
                    log.warn("Master | Terminating %s", job_id)
                    job = running[job_id]
                    job.process.terminate()
                    job.process.join()
                    job.process = None
                    slots[0] = slots[0] + job.threads
                    del running[job_id]
                elif job_id in jobs:
                    job = jobs[job_id]
                    del jobs[job_id]
                if job:
                    update_dependencies(job)
                    remove_children(job)
                    schedule()

            if wait and (len(jobs) + len(running) == 0):
                break
        log.info("Master | Master loop terminated")

    @staticmethod
    def _terminate_process(process, log):
        """Helper function that terminates the given process"""
        if process is not None:
            log.info("Exec | Sending SIGTERM")
            process.send_signal(signal.SIGTERM)
            #process.terminate()
            # check if the job is dead. if not
            # sleep for a moment and check again.
            if process.poll() is None:
                # give it 5 seconds to cleanup and exit
                import time
                for t in [0.01, 0.01, 0.01,
                          0.02, 0.02, 0.02,
                          0.02, 0.02, 0.02,
                          0.10, 1, 2, 5]:
                    time.sleep(t)
                    if process.poll() is not None:
                        log.info("Exec | Processes terminated after SIGTERM")
                        break
                else:
                    # nothing worked, kill the job
                    log.info("Exec | Processes still running, sending SIGKILL")
                    os.kill(process._popen.pid, signal.SIGKILL)

    @staticmethod
    def execute(requests, job_id, job):
        """Local grid executor method that takes a ``LocalJob`` instance and
        runs its command in a subprocess.

        :param request: the request queue to send completion state back to
                        the master
        :param job_id: the local job id
        :param job: instance of local job
        """
        log = getLogger("jip.grids.Executor")
        process = None

        # setup local signal handling
        def handle_term(sig, frame):
            log.error("Exec | Received termination signal")
            LocalCluster._terminate_process(process, log)
            requests.put(['DONE', job_id, 1])
            sys.exit(1)

        signal.signal(signal.SIGTERM, handle_term)

        log.info("Exec | Start job %s", job_id)
        result = -1
        try:
            cwd = job.working_directory
            stdout = open(job.stdout, 'w')
            stderr = open(job.stderr, 'w')
            cmd = job.cmd
            process = subprocess.Popen("exec " + cmd,
                                       stdout=stdout,
                                       stderr=stderr,
                                       shell=True,
                                       cwd=cwd)
            result = process.wait()
        except Exception as err:
            log.error("Exec | Error : %s", err, exc_info=True)
            requests.put(['FAILED', job_id, str(err)])
        log.info("Exec | Job finished %s with %d", job_id, result)
        requests.put(['DONE', job_id, result])
