#!/usr/bin/env python
"""
The JIP job interpreter command can be used to execute jobs on the local host

Usage:
   jip-interpreter [-p] <file> [<args>...] [-- [submit] [<jip_args>...]]
   jip-interpreter [-p] <file> [--help|-h]
   jip-interpreter [--help|-h]

Options:
    -p,--pipeline   the file contains a pipeline
    <file>          the jip script
    <args>          script arguments
    <jip_args>      additional jip argument that can be passed to
                    the jip runtime. To see a list of options, take
                    a look at the 'jip run' command or the 'jip submit'
                    command.
    submit          forwards the execution to 'jip submit' and sends the script
                    to a cluster. By default, the script is executed locally
    -h,--help       show this help message

"""

from . import parse_args


def split_to_jip_args(args):
    """Check the <args> and search for '--'. If found,
    everything after '--' is put into 'Jip_args'"""
    if args and "<args>" in args:
        try:
            i = args["<args>"].index("--")
            args["<args>"], args["<jip_args>"] = args["<args>"][:i], \
                args["<args>"][i + 1:]
        except ValueError:
            pass


def main():
    initial = parse_args(__doc__)
    args = dict(initial)
    # split args and jip_args
    split_to_jip_args(args)

    script_file = args["<file>"]
    script_args = args["<args>"]
    jip_args = args["<jip_args>"]
    if args['--pipeline']:
        jip_args.append('--pipeline')
    if "--" in jip_args:
        jip_args.remove("--")
    # execute locally or submit
    if "submit" in jip_args:
        from jip.cli.jip_submit import main as jip_submit
        jip_args.remove("submit")
        jip_submit(jip_args + [script_file] + script_args)
    else:
        from jip.cli.jip_run import main as jip_run
        jip_run(jip_args + [script_file] + script_args)


if __name__ == "__main__":
    main()
