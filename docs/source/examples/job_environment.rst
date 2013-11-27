.. _tut_job_env:

Modify the execution environment
================================
All executable units in the JIP system are attached to a specific 
:py:class:`~jip.profiles.Profile` that can be modified from within the tool 
implementation as well as at submission or execution time.

The Profile
-----------
The :py:mod:`jip.profiles` module documentation covers the properties that
you can set and modify on the job profile. Here, we will go over some of 
them to explain how and in which order a certain profile is applied to a job
at execution time. 

Job profiled of a tool can be accessed from within the tools *validation* 
block. This is important because the properties of the profile are applied
as submission time properties. That means they have to be set when the tool
is *submitted* to your compute cluster. 
