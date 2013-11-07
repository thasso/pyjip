JIP-API reference documentation
===============================

The JIP platform is mostly written in Python, except for the stream dispatcher,
which is written in C and integrated as a Python Extension. This documentation
covers the JIP API and describes the basic modules and classes that make up the
system. 

In a lot of cases it will not be necessary to read and understand the full API
reference. It might come in handy though when you are in the situation of
extending the system, for example, adding support for your own cluster, or if
you want to dig deeper and see how things are created. 

.. toctree::
   :maxdepth: 2
          
   jobs
   tools
   options
   utils
   executils
   cli

