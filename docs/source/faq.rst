:tocdepth: 2

=======
Jip FAQ
=======

.. only:: html

   .. contents::


API related questions
---------------------

How can I define and use a tool?
********************************

Have a look at the following code::

    @tool()
    def nop():
        return ""

this defines a tool called ``nop``, by default the function name is used, you
can drop one in the ``@tool`` :ref:`decorator <decorators>`.

The pipeline ``run()`` function takes a string or a tool instance.
The string is used as tool name and a search is performed.

In this example, this would be enough to get the tool and run it::

    >>> p = Pipeline()
    >>> p.run('nop')

The ``run()`` method returns a :py:class:`~jip.pipelines.Node` object and also
exposes all the properties of the node including its options.


How can I access the tool options?
**********************************
Following the example before we can get a :class:`~jip.pipelines.Node` object
from a pipeline and assign options to it::

    >>> p = Pipeline()
    >>> n = p.run(tool)
    >>> n.input = "Makefile"
    >>> n.output = "out.txt"
    >>> n.inter = "inter.out"

The node has a ``set()`` function as a fall-back if direct assignment does not
work, for example, if your option name conflicts with a function name of the
Node object. In that case the Node function are preferred.

You can always access the tool options using the ``self.args`` and
``self.options`` properties. The ``args`` dict contains the raw values while
the ``options`` dict returns an instance of :py:class:`~jip.options.Option`.

If you want to establish dependencies between tools in a pipeline, you have to
use ``options`` but if you just want to access values ``args`` is fine.

.. note:: Note that options contain the raw values and are able to render them.
          The rendering is performed at the latest time possible though.


How can I access local variables within a pipeline?
***************************************************

Rendering values with access to the local context is done explicitly for JIP
scripts. If you want to access the local context in python you need to pass the
``locals()`` to the pipeline.

For example::

    >>> p = jip.Pipeline()
    >>> a = "Makefile"
    >>> p.bash("wc -l ${a}")
    bash
    >>> b = jip.create_jobs(p)[0]
    >>> assert b is not None
    >>> assert b.command == '(wc -l Makefile)'
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    AssertionError
    >>>

Note that the nodes ``cmd`` option references "a", a local variable. To resolve
local variables you should drop the ``locals()`` into the pipeline context 
function before you return the pipeline or do anything with it::

    >>> p = jip.Pipeline()
    >>> a = "Makefile"
    >>> p.bash("wc -l ${a}")
    bash
    >>> p.context(locals())
    >>> b = jip.create_jobs(p)[0]
    >>> assert b is not None
    >>> assert b.command == '(wc -l Makefile)'
