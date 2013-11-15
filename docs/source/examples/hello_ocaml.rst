Argument filters and Streams
============================
This small example demonstrates three features of JIP. 

    * Any interpreted language can be used to write executable blocks

    * The ``arg`` and ``else`` :ref:`template filters <template_filters>`
      can be handy when it comes to the tedios task of argument parsing and
      interpretation

    * JIP's :ref:`stream dispatcher <stream_dispatching>` allows you to 
      write streams to files and other processes in parallel.

Other interpreters
------------------
The default interpreter used in JIP templates and scripts is ``bash``, but you 
can change the interpreter easily. For this, you simply specify the name or
the path to the interpreter as a command block argument, for example, 
``#%begin command perl`` will open a *perl* interpreted command block. 

Here is a full example where we use ``OCaml`` to drive our tool:

.. code-block:: ocaml
   :emphasize-lines: 12,12

    #!/usr/bin/env jip
    #Send a greeting from ocaml
    #
    # usage:
    #     hello  -n <name> [-o <output>]
    #
    # Options:
    #   -n, --name <name>      Your name
    #   -o, --output <output>  The output
    #                          [default: stdout]

    #%begin command ocaml
    let o = ${output|arg('open_out "', '"')|else("stdout")};;
    Printf.fprintf o "Dear ${name}\n\n";;
    Printf.fprintf o "From the happy chambers of my spirit,\n";;
    Printf.fprintf o "I send forth to you, my smile.\n";;
    Printf.fprintf o "And, I pray that you shall carry it with you,\n";;
    Printf.fprintf o "as you travel each and every mile.\n";;

In order to switch to OCaml, simply start the command block with 
``#%begin command ocaml``.

Argument filtering
------------------
In the example above, the tool exposes the ``output`` option to the user and
defaults to ``stdout``. Here we use the ``arg`` and ``else`` 
:ref:`template filters <template_filters>` to solve the problem. Somehow we 
have to wither select ``stdout`` or open the user specified file. The line
that solves the problem is::

    let o = ${output|arg('open_out "', '"')|else("stdout")};;

What happens here? First of all, everything outside of ``${}`` is pure OCaml.
Inside the braces, the following logic is applied. Take the value of the
``output`` options. If a value was assigned that is not a file stream and does
not evaluate to false, pass it to ``arg`` (see
:py:func:`jip.templates.arg_filter`) and insert the result. Otherwise, pass it
through ``else``. The ``else`` argument of else is inserted if the passed value
is a file stream or evaluates to false.

In case the user specified an output file, we have to surround it with quotes
and pass it to OCamls ``open_out`` function. Here we solve this by specifying 
the ``arg`` filters ``prefix`` and ``suffix`` arguments.

This results in the following evaluated output. If the user did not specify 
an output file, ``arg`` takes the value and returns it unmodified. In this case,
the ``else`` block is the one that generates the final result::

    let o = stdout;;

If an output file was specified, say ``myfile.txt``, the ``argfilter`` will 
surround it with the specified prefix and suffix and create the final result
as::

    let o = open_out "myfile.txt"


Dispatching streams
-------------------
Lets take the ``hello`` tool from the previous example and create a small, 
indeed not very useful, counting pipeline::

    #!/usr/bin/env jip -p
    # Count words and lines
    #
    # usage:
    #     counter <name>

    greetings = run('hello', name=args['name'])
    line_count = bash('wc -l', input=greetings)
    full_count = bash('wc', input=greetings)

.. note:: In this example, we do not open a block with 
          ``#%begin pipeline`` to start implementing our pipeline. Instead,
          we pass the ``-p`` parameter to the *jip* interpreter in the shebang
          line. This switches to pipeline mode automatically. 

The pipeline in this example is rather straight forward. We take a single 
``name`` argument. Then we call ``hello``, our OCaml tool from the previous
example. Next, the output of ``hello``, called ``greetings``, is passed as
input to two ``bash`` tools, one that does a full count, and one that only 
counts lines.

Lets perform a dry run to see what happens. Don't worry, we will not start the 
pipeline, so there is no need to have OCaml installed::

    $> ./counter.jip Joe -- --dry
    ...
    ############################################################...
    |                                                           ...
    +--------------------------------+--------+-----------------...
    |              Name              | State  |                 ...
    +================================+========+=================...
    | greetings|line_count|full_coun | Hold   |                 ...
    | t                              |        |                 ...
    +--------------------------------+--------+-----------------...
    ####################
    |  Job hierarchy   |
    ####################
    greetings
    ├─line_count
    └─full_count
    ####################

The dry run screen for this pipeline show is two things.

First, the dependencies are correct as you can see from the *Job Hierarchy*.
*Greetings* is the primary job and executed first and the other two jobs depend
on it.

Second, the *Job states* table shows only a single job called
``greetings|line_count|full_count``. This indicates that all three jobs have 
*streams* connecting them and they all have to be executed in parallel. In 
case you want to submit the pipeline to a compute cluster, this means all three
jobs have to be executed on the same node in a single *cluster job*. The reason
for this behavior is obvious. Our ``hello`` tool prints its output to 
``stdout`` and the two counter tools read from ``stdin``.

Lets modify the pipeline a little bit and write the output of our call to 
``hello`` into a result file. The rest of the pipeline stays untouched:

.. code-block:: python
   :emphasize-lines: 1,1

    greetings = run('hello', name=args['name'], output='result.txt')
    line_count = bash('wc -l', input=greetings)
    full_count = bash('wc', input=greetings)

Watch what happens when we perform a dry run::

    $> ./counter.jip Joe -- --dry
    ...
    ###########################################################...
    |                                                          ...
    +--------------------------------+--------+----------------...
    |              Name              | State  |                ...
    +================================+========+================...
    | greetings                      | Hold   |                ...
    | line_count                     | Hold   | result.txt     ...
    | full_count                     | Hold   | result.txt     ...
    +--------------------------------+--------+----------------...
    ####################
    |  Job hierarchy   |
    ####################
    greetings
    ├─line_count
    └─full_count
    ####################

The job hierarchy stays untouched as we did not modify any of the dependencies,
but instead of a single job, all three tools are now executed in dedicated
jobs. The table already shows the reason. The two counter jobs are now operating
on the output file of the ``greetings`` job. ``greetings`` has to finish first,
but the two counters can now be executed in two separate jobs. 

If this is beneficial or not depends on the tasks and on your compute 
infrastructure. For example, if one of the secondary jobs is able to work
multi-threaded while the other one uses only a single CPU, it might be nice 
to split them into two dedicated job on your cluster. On the other hand, if
your data stream is quiet large and the task is not very computationally 
intense, it might be better stream the data through all the jobs and submit
only a single job to you your cluster. 

Now, what if you want the results of the initial job to be stored in 
your ``results.txt`` output file but you still would like to stream the 
data through to your other two jobs. In pure bash, you might solve it with
bashs  ``tee`` command. In JIP, you can achieve the same goal by adding a 
single line to your pipeline that ensures the stream pipe behaviour:

.. code-block:: python
   :emphasize-lines: 5,5

    greetings = run('hello', name=args['name'], output='result.txt')
    line_count = bash('wc -l', input=greetings)
    full_count = bash('wc', input=greetings)

    greetings | (line_count + full_count)

The last line of the pipeline recreates the pipes, but **keeps the results 
file** as part of the stream. The output of ``greetings`` will be **streamed**
to ``results.txt`` and to the two counter jobs.








