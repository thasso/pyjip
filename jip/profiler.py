#!usr/bin/env python
"""The profiler module is used to profile a single job executed in
the jip environment and store its results.

The profiler is executed using the multiporcessing module. The results are
stored to a per-job file in a table format covering the following information
in the specified order:

    1) time stamp as seconds since epoch
    2) pid
    3) name
    4) cpu usage in percent
    5) cpu user
    6) cpu system
    7) memory percent
    8) memory rss
    9) memory vms
    10) box cpu percent user
    11) box cpu percent system
    12) box cpu percent idle
    13) box disk io read bytes
    14) box disk io write bytes
    15) box net io read bytes
    16) box net io write bytes
    17) box memory used
    18) box memory free
    19) box physical memory free
    20) box total memory
"""
import os
import multiprocessing
import logging
import time

log = logging.getLogger("jip.profiler")


class Profiler(object):
    """The Profilter class takes care of the actual process
    checking. You can initialize the profiler with a process to watch,
    the corresponding job and an optional interval. The interaval is
    given in seconds.

    The profiling will not start automatically. You have to start the
    profiler explicitly using the :py:meth:`start` method.

    :param process: the processes to watch
    :param job: the job that was used to create the process
    :param interval: check interval in seconds
    """
    def __init__(self, process, job, interval=5):
        self.process = process
        self.job = job
        self.interval = interval

    def _collect_and_write(self, process, writer, ts=None, iostat=None,
                           netstat=None):
        import psutil
        if ts is None:
            ts = time.time()
        cpu_percent = process.get_cpu_percent()
        cpu_times = process.get_cpu_times()
        mem_info = process.get_memory_info()
        data = [
            int(ts),
            process.pid,
            process.name,
            cpu_percent,
            cpu_times.user,
            cpu_times.system,
            process.get_memory_percent(),
            mem_info.rss,
            mem_info.vms
        ]

        box_data = None
        if iostat and netstat:
            box_cpu = psutil.cpu_times_percent()
            io = psutil.disk_io_counters()
            netio = psutil.net_io_counters()
            box_mem = psutil.phymem_usage()
            box_data = [
                box_cpu.user,
                box_cpu.system,
                box_cpu.idle,
                io.read_bytes - iostat.read_bytes,
                io.write_bytes - iostat.write_bytes,
                netio.bytes_recv - netstat.bytes_recv,
                netio.bytes_sent - netstat.bytes_sent,
                box_mem.used,
                box_mem.free,
                (psutil.used_phymem() - psutil.cached_phymem()),
                box_mem.total
            ]
        else:
            box_data = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        data.extend(box_data)
        # write data
        print >>writer, "\t".join([str(s) for s in data])
        # process children
        for child in process.get_children():
            self._collect_and_write(child, writer, ts=ts)
        pass

    def _run(self):
        """Internal method that is massed to the multiprocessing
        process"""
        import psutil
        log.info("Profiler | start profiling %s with job id %d",
                 self.job, self.job.id)
        # base directory
        folder = self.job.working_directory
        # create base file
        stats_file = os.path.join(folder, "%s-%s_jip-profile.dat" % (
            str(self.job.id),
            str(self.job.job_id) if self.job.job_id else "0"
        ))
        log.info("Profiler | store profile in %s", stats_file)
        process = psutil.Process(self.process.pid)
        iostat = psutil.disk_io_counters()
        netstat = psutil.net_io_counters()
        with open(stats_file, 'w') as writer:
            # write header
            print >>writer, "#time\tpid\tname\tcpu_percent\tcpu_user" \
                            "\tcpu_system\t" \
                            "mem_percent\tmem_rss\tmem_vms\t" \
                            "box_cpu_user\tbox_cpu_system\tbox_cpu_idle\t" \
                            "disk_read\tdisk_write\tnet_read\tnet_write\t" \
                            "box_mem_used\tbox_mem_free\tbox_phy_mem_free\t" \
                            "box_mem_tot"

            while True:
                try:
                    self._collect_and_write(process, writer, iostat=iostat,
                                            netstat=netstat)
                except psutil.NoSuchProcess:
                    break
                try:
                    process.wait(self.interval)
                    break
                except:
                    pass

    def start(self):
        """Start the profiler in a background task."""
        try:
            import psutil
        except:
            log.error("Unable to import the 'psutil' module. "
                      "You have to install the module before the "
                      "profiler can be used. Try a `pip install psutil`."
                      "!Profiling disabled!")
            return
        p = multiprocessing.Process(target=self._run)
        p.start()
        return p
