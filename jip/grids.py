#!/usr/bin/env python
"""JIP ships with a **small and simple** local queueing system.
"""
import multiprocessing
import os
import subprocess
import signal
import sys
import logging

import jip.cluster
from jip.logger import getLogger


class LocalCluster(jip.cluster.Cluster):

    def __init__(self, _start=True,
                 _remote_ids=False,
                 cores=None, log_level=logging.INFO):
        self._current_id = 0
        self.log = getLogger("jip.grids.LocalCluster")
        self.master_requests = None
        self.master_response = None
        self.master_process = None
        self._current_id
        self._remote_ids = _remote_ids
        self.cores = cores
        self.loglevel = log_level
        # set log level
        self.log.setLevel(log_level)
        # start the mater process
        if _start:
            self.start()

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

    def start(self):
        if self.master_process is None:
            self.log.info("Starting local cluster master")
            self.master_requests = multiprocessing.Queue()
            self.master_response = multiprocessing.Queue()
            self.master_process = multiprocessing.Process(
                target=_GridMaster.create_master,
                args=[self.master_requests, self.master_response],
                kwargs={"cores": self.cores,
                        "loglevel": self.loglevel},
                name="grid-master"
            )
            self.master_process.start()

    def shutdown(self):
        if self.master_process is None:
            return
        self.log.info("Requesting master shutdown")
        self.master_requests.put(['EXIT'])
        self.master_process.join()
        self.master_requests.close()
        self.master_response.close()
        self.master_process = None

    def wait(self):
        if self.master_process is None:
            return
        self.log.info("Waiting for termination of the master")
        self.master_requests.put(["WAIT"])
        self.master_process.join()
        self.master_requests.close()
        self.master_response.close()

    def list(self):
        if self.master_process is None:
            return []
        self.master_requests.put(["JOBS"])
        return self.master_response.get()

    def submit(self, job):
        if self.master_process is None:
            raise jip.cluster.SubmissionError("No Grid master found!")
        local_job = _Job.from_job(job) if not isinstance(job, _Job) else job
        if not self._remote_ids:
            local_job.job_id = self._next_id()

        self.master_requests.put([
            "SUBMIT", local_job,
        ])
        if self._remote_ids:
            job.job_id = self.master_response.get()
        else:
            job.job_id = local_job.job_id
        self.log.info("Submitted new job with id %s", job.job_id)
        return job

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("%J", str(job.id))

    def cancel(self, job):
        if self.master_process is None:
            return
        self.log.info("Send cancel request for %s", job.job_id)
        self.master_requests.put(["CANCEL", job.job_id])


class JIP(jip.cluster.Cluster):

    def __init__(self):
        self.log = getLogger("jip.grids.JIP")
        cfg = jip.config.get('jip_grid', {})
        self.host = cfg.get('host', '127.0.0.1')
        self.port = cfg.get('port', '5556')
        self.log = getLogger("jip.cluster.JIP")

    def _connect(self):
        import zmq
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.setsockopt(zmq.LINGER, 0)
        socket.RCVTIMEO = 1000
        self.log.info("Connecting to remote grid: %s:%s", self.host, self.port)
        socket.connect("tcp://%s:%s" % (self.host, self.port))
        return context, socket

    def list(self):
        try:
            context, socket = self._connect()
            socket.send_json({"cmd": "list"})
            jobs = socket.recv_json()
        except:
            self.log.error("Unable to connect to grid server!")
            raise
        finally:
            socket.close()
            context.term()
        return jobs

    def submit(self, job):
        try:
            context, socket = self._connect()
            local = _Job.from_job(job)
            socket.send_json({"cmd": "submit", "job": {
                "job_id": local.job_id,
                "id": local.id,
                "cmd": local.cmd,
                "stdout": local.stdout,
                "stderr": local.stderr,
                "working_directory": local.working_directory,
                "dependencies": [f for f in local.dependencies],
                "children": [f for f in local.children]

            }})
            submitted = socket.recv_json()
            job.job_id = submitted['job_id']
            job.working_directory = submitted['working_directory']
            job.stdout = submitted['stdout']
            job.stderr = submitted['stderr']
            return job
        except:
            self.log.error("Unable to connect to grid server!", exc_info=True)
            raise jip.cluster.SubmissionError("Unable to connect to cluster")
        finally:
            socket.close()
            context.term()
        return job

    def resolve_log(self, job, path):
        if path is None:
            return None
        return path.replace("%J", str(job.id))

    def cancel(self, job):
        try:
            context, socket = self._connect()
            self.log.info("Send cancel request for %s", job.job_id)
            socket.send_json({"cmd": "cancel", 'id': job.job_id})
            socket.recv()
        except:
            self.log.error("Unable to connect to grid server!")
            raise
        finally:
            socket.close()
            context.term()


class _Job(object):
    """Local Job object that is send to the local grid and contains
    all the information that is needed and supported by the local cluster
    implementation.

    This class is for internal purposes and there should be no need to
    interact with it outside of this module.

    This class is sortable and the order is defined by the number of activate
    dependencies and the job_id.
    """
    def __init__(self, cmd=None, cwd=None, stdout=None, stderr=None,
                 dependencies=None, threads=1, job_id=None):
        self.id = job_id
        self.job_id = job_id
        self.dependencies = [] if dependencies is None else dependencies
        self.threads = threads
        self.cmd = cmd
        self.working_directory = cwd
        self.stdout = stdout
        self.stderr = stderr
        self.children = []
        self.process = None

    @classmethod
    def from_job(cls, job):
        cwd = job.working_directory if job.working_directory is not None \
            else os.getcwd()
        if job.stderr is None:
            job.stderr = os.path.join(cwd, "jip-%J.err")
        if job.stdout is None:
            job.stdout = os.path.join(cwd, "jip-%J.out")

        #collect dependencies
        deps = set([])
        for dep in [d for d in job.dependencies if d.job_id]:
            deps.add(int(dep.job_id))

        # build the local job
        local_job = _Job(
            job.get_cluster_command(),
            cwd,
            job.stdout,
            job.stderr,
            dependencies=deps,
            threads=job.threads
        )
        local_job.id = job.id
        return local_job

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


class _GridMaster(object):
    """The grid master instance"""
    def __init__(self, requests, response, cores=None, loglevel=logging.INFO):
        """Initialize a new grid master with the request and response
        queues and optionally the number of cores or slots available.
        If no cores are specified, the number of cores of the current
        machine are used
        """
        #: the queue this master receives requests on
        self.requests = requests
        #: the queue this master sends responses to
        self.response = response
        #: currently queued jobs
        self.queued = {}
        #: currently running jobs
        self.running = {}
        #: number of currently available slots
        self.slots_available = multiprocessing.cpu_count() \
            if not cores else cores
        #: number of totally available slots
        self.slots_total = self.slots_available
        #: the master logger
        self.log = getLogger("jip.grids.Master")
        #: self wait mode
        self.wait_mode = False
        #: current id
        self._current_id = 0

        self.log.setLevel(loglevel)
        self.log.info("Master | Initialized with %d available slots",
                      self.slots_available)

    def _next_id(self):
        """Increment the current id and return the next available job
        id.

            >>> c = _GridMaster(None, None)
            >>> assert c._current_id == 0
            >>> assert c._next_id() == 1
            >>> assert c._next_id() == 2

        :returns: next id
        """
        self._current_id += 1
        return self._current_id

    def _remove_children(self, job):
        """Given a _Job instance, this recursively removes all child
        jobs from the list of queued jobs. This ONLY works on children
        that are NOT currently running.

        :param job: the parent job instance
        """
        # recursive remove of all child jobs
        # that are still queued
        self.log.debug("Master | Removing all children of %s from queued jobs",
                       job.job_id)
        for c in job.children:
            if c in self.queued:
                child = self.queued[c]
                self._remove_children(child)
                self._update_dependencies(child)
                del self.queued[c]

    def _update_dependencies(self, job):
        """Given a _Job, this updates the dependency list of all
        children of the given job and removes the given job from the
        list of dependencies.

        The dependency list of a jobs is used to determine if the job can
        be executed. It can start if there are no dependencies left. This
        removes the given job from the dependencies and should be called
        after the given job finished successfully.

        :param job: the parent job that will be removed from its child
                    dependencies
        """
        for c in job.children:
            if c in self.queued:
                self.log.debug("Master | Removing %s from child %s "
                               "dependencies", job.job_id, c)
                try:
                    self.queued[c].dependencies.remove(job.job_id)
                except KeyError:
                    self.log.error("Job with id %s no found in the child %s "
                                   "dependencies!", job.job_id, c)

    def _run_job(self, job):
        """Takes the given local job and starts it in a dedicated process
        with a new Executor instance to run the job.
        The jobs process is stored in the job instance and the job instances
        is moved from queued to running

        :param job: the job to send
        """
        if not job.job_id in self.queued:
            self.log.error("Job %s not in the set of queued jobs!", job.job_id)
            return

        job_id = job.job_id
        self.log.info("Master | Starting job %s", job_id)
        process = multiprocessing.Process(
            target=_execute_job,
            args=[self.requests, job_id, job],
            name="worker-%s" % str(job_id)
        )
        job.process = process
        # move the job to running
        self.running[job_id] = job
        # remove the job from the queues jobs
        del self.queued[job_id]
        # start the process
        process.start()

    def schedule(self):
        self.log.info("Master | running scheduler %d/%d",
                      self.slots_available, self.slots_total)
        # find next jobs that can be executed and run it
        #
        # sort all queued jobs
        sorted_jobs = sorted(self.queued.values())
        # now iterate the jobs and run the first one that fits
        # into the available slots
        for job in sorted_jobs:
            if len(job.dependencies) > 0:
                # not all dependencies are resolved
                self.log.info("Master | No jobs without dependencies found")
                break
            if job.threads <= self.slots_available:
                self.log.info("Master | Submitting job for execution %s",
                              job.job_id)
                # start this one
                self._run_job(job)
                # update slots
                self.slots_available -= job.threads
                # if we still have slots open,
                # start another scheduling
                if self.slots_available >= 1:
                    self.schedule()
                return

    def _resolve_log(self, job_id, path):
        if path is None:
            return None
        return path.replace("%J", str(job_id))

    def _handle_exit(self, *args):
        """The EXIT handler"""
        self.log.info("Master | EXIT request, shutting down")
        for job_id, p in self.running.iteritems():
            self.log.warn("Master | Terminating %s", job_id)
            p.process.terminate()
            p.process.join()
        return False

    def _handle_wait(self, *args):
        """The WAIT handler"""
        if not self.wait_mode:
            self.log.info("Master | WAIT requested, entering wait mode!")
            self.wait_mode = True
        return True

    def _handle_jobs(self, *args):
        """The JOBS handler"""
        self.log.info("Master | list jobs")
        # queued and running
        self.response.put(list(self.queued.keys()) +
                          list(self.running.keys()))
        return True

    def _handle_submit(self, *args):
        """The SUBMIT handler"""
        job = args[1]
        create_id = job.job_id is None
        job_id = self._next_id() if create_id else job.job_id
        job.job_id = job_id
        job.stdout = self._resolve_log(job.id, job.stdout)
        job.stderr = self._resolve_log(job.id, job.stderr)
        self.queued[job_id] = job

        # update children
        for d in job.dependencies:
            if d in self.queued:
                self.queued[d].children.append(job_id)
            if d in self.running:
                self.running[d].children.append(job_id)
        self.log.info("Master | Queue new job %s", job_id)
        if create_id:
            self.response.put(job_id)
        if self.slots_available >= 1:
            self.schedule()
        return True

    def _handle_done(self, *args):
        """The DONE handler"""
        job_id = int(args[1])
        state = args[2]
        if job_id in self.running:
            job = self.running[job_id]
            job.process.join()
            job.process = None
            self.slots_available += job.threads
            self._update_dependencies(job)
            if state != 0:
                self.log.error("Master | Job %s failed with %s",
                               job_id, state)
                self._remove_children(job)
            else:
                self.log.info("Master | Job %s finished with %s",
                              job_id, state)
            job.process = None
            del self.running[job_id]
            self.schedule()
        else:
            self.log.warn("Master | Job %s marked as done, but was "
                          "not found in running jobs!", job_id)
        return True

    def _handle_failed(self, *args):
        """The FAILED handler"""
        job_id = int(args[1])
        error = args[2]
        self.log.error("Master | Execution of %s failed: %s", job_id, error)
        if job_id in self.unning:
            job = self.running[job_id]
            job.process = None
            self._update_dependencies(job)
            self._remove_children(job)
            self.slots_available += job.threads

            del self.running[job_id]
            self.schedule()
        else:
            self.log.error("Master | Execution of %s marked as failed, but "
                           "was not found in running jobs!" % job_id)
        return True

    def _handle_cancel(self, *args):
        """The CANCEL handler"""
        job_id = int(args[1])
        self.log.warn("Master | Execution of %s canceled", job_id)
        job = None
        if job_id in self.running:
            self.log.warn("Master | Terminating %s", job_id)
            job = self.running[job_id]
            job.process.terminate()
            job.process.join()
            job.process = None
            self.slots_available += job.threads
            del self.running[job_id]
        elif job_id in self.queued:
            job = self.queued[job_id]
            del self.queued[job_id]
        else:
            self.log.warn("Master | Execution %s cancellation requested, but "
                          "job could not be found!", job_id)
        if job:
            self._update_dependencies(job)
            self._remove_children(job)
            self.schedule()
        return True

    def start(self):
        """Start the master loop and get messages from the requests queue.
        The command messages are lists where the first element is the
        actual command, followed by the command parameters.
        """

        handlers = {
            "EXIT": self._handle_exit,
            "JOBS": self._handle_jobs,
            "SUBMIT": self._handle_submit,
            "WAIT": self._handle_wait,
            "DONE": self._handle_done,
            "FAILED": self._handle_failed,
            "CANCEL": self._handle_cancel,
        }

        while True:
            if self.wait_mode and self._num_jobs() == 0:
                break
            try:
                msg = self.requests.get()
            except KeyboardInterrupt as err:
                self.log.warn("Server process canceled")
                continue
            self.log.debug("Master | received command: %s", msg)
            try:
                command = msg[0]
                handler = handlers[command]
                if not handler(*msg):
                    break
            except Exception as err:
                self.log.error("Master | error while handling message %s : %s",
                               msg, str(err), exc_info=True)
                break
        self.log.info("Master | Master loop terminated")

    def _num_jobs(self):
        return len(self.queued) + len(self.running)

    @staticmethod
    def create_master(request, response, cores=None, loglevel=logging.INFO):
        master = _GridMaster(request, response, cores=cores,
                             loglevel=loglevel)
        master.start()


def _execute_job(requests, job_id, job):
    """Local grid executor method that takes a ``_Job`` instance and
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
        _terminate_process(process, log)
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
        process = subprocess.Popen(
            "exec " + cmd,
            stdout=stdout,
            stderr=stderr,
            shell=True,
            cwd=cwd
        )
        result = process.wait()
    except Exception as err:
        log.error("Exec | Error : %s", err, exc_info=True)
        requests.put(['FAILED', job_id, str(err)])
    log.info("Exec | Job finished %s with %d", job_id, result)
    requests.put(['DONE', job_id, result])


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
