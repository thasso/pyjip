#!/usr/bin/env python
"""
Start the JIP grid server.

This allows you to schedule jobs on your local machine. In order to use
this sercice, start the server and make sure you configure your jip client
to use the local cluster. Put this in your ~/.jip/jip.json configuration
file:

    {
        "cluster": "jip.grids.JIP",
        "jip_cluster":{
            "port": 5556
        }
    }

Usage:
   jip-server--help|-h] [-p <port>]

Options:
    -p, --port <port >  The port used for the server
                        [default: 5556]

Other Options:
    -h --help             Show this help message
"""


from jip.logger import getLogger
from . import parse_args
import jip.grids
import sys
import jip.db

log = getLogger("jip.cli.jip_server")


def main():
    args = parse_args(__doc__, options_first=True)

    try:
        import zmq
    except:
        print >>sys.stderr, """
Unable to import the python ZeroMQ binding.
Please make sure that zeromq is installed on your system.

You can install zeroMQ using pip:

    pip intall pyzmq
"""
        sys.exit(1)

    try:
        port = args['--port']
        log.info("Starting JIP grid server no port %s", port)
        context = zmq.Context()
        socket = context.socket(zmq.REP)
        socket.bind("tcp://*:%s" % port)
        log.info("Socket server started")

        cluster = jip.grids.LocalCluster()
        while True:
            msg = socket.recv_json()
            log.info("Received message: %s" % msg)
            if msg['cmd'] == 'list':
                jobs = cluster.list()
                log.debug("Recevied cluster jobs: %s", jobs)
                socket.send_json(jobs)
            elif msg['cmd'] == 'cancel':
                job_id = msg['id']
                job = jip.db.Job()
                job.job_id = job_id
                cluster.cancel(job)
                log.debug("Canceled job: %s", job_id)
                socket.send('ACK')
            elif msg['cmd'] == 'submit':
                job_dict = msg['job']
                job = jip.grids._Job()
                job.__dict__.update(job_dict)
                job = cluster.submit(job)
                log.debug("Submitted job: %s", job.job_id)
                socket.send_json({
                    "job_id": job.job_id,
                    "id": job.id,
                    "cmd": job.cmd,
                    "stdout": job.stdout,
                    "stderr": job.stderr,
                    "working_directory": job.working_directory,
                    "dependencies": [f for f in job.dependencies],
                    "children": [f for f in job.children]
                })
            else:
                log.error("Unknown command: %s" % msg)
                socket.send("UNKNOWN")
    except KeyboardInterrupt:
        log.warn("Shutting down")
        cluster.shutdown()
        socket.close()
        context.term()
    except Exception as e:
        log.error("Error running JIP server: %s", str(e), exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

