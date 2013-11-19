.. _jip_tools_cmd:

jip tools - list available tools
================================

Synopsis
--------

**jip tools** [--help|-h]

Description
-----------
List all JIP tools/scripts that are available in the search paths.

Tools scripts
-------------
All script found in the current search paths are listed. Please note that there
might be more. Here, we search only for files with the .jip extension!

The following locations are searched:

* Current directory
* Jip configuration (jip_path)
* `JIP_PATH` environment variable

Tools implemented in Python modules
-----------------------------------
The modules must be available in `PYTHONPATH` and must be specified in
the jip configuration or in the `JIP_MODULES` environment variable.
Please note that pipeline scripts that contain python blocks are allowed to
load modules that contain tool implementation. These tools might not be found
by this scan!

The following locations are searched:

* Jip configuration (jip_moduled)
* `JIP_MODULES` environment variable

Options
-------

-h, --help Show this help message

