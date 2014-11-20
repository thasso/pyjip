#!/usr/bin/env python
from jip import *


####################################################
# Implementation as single functions either as tools
# that return a script or as (interpreter, script)
# tuple. You can also decorate a function with
# @pytool and implement the tool execution in python
# directly.
####################################################
@tool()
def fun_hello_world():
    return "echo 'Hello World'"


@tool()
def fun_hello_world_perl():
    return 'perl', "print('Hello World\n');"


@pytool()
def fun_hello_world_py():
    """Prints hello world in a python module"""
    print("Hello python")


####################################################
# Implementation as python classes. For a simple
# hello world, that might seem overhead, but class
# base implementation allow you to customize the
# initialization, setup, and validation
####################################################
@tool()
class cls_hello_world(object):
    def get_command(self):
        return "echo 'Hello World'"


@tool()
class cls_hello_world_perl(object):
    def get_command(self):
        return 'perl', "print('Hello World\n');"


@pytool()
class cls_hello_world_py(object):
    def run(self):
        print("Hello World")
