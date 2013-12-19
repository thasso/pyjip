#!/usr/bin/env python
from jip import *


def test_name_function_single_tool():
    @tool()
    class MyTool(object):
        def setup(self):
            self.name("Custom-Name")

        def get_command(self):
            return "true"
    pipeline = Pipeline()
    pipeline.run('MyTool')
    jobs = create_jobs(pipeline)
    assert jobs[0].name == 'Custom-Name'


def test_name_attr_single_tool():
    @tool()
    class MyTool(object):
        def setup(self):
            self.job.name = "Custom-Name"

        def get_command(self):
            return "true"
    p = Pipeline()
    p.run('MyTool')
    jobs = create_jobs(p)
    assert jobs[0].name == 'Custom-Name'


def test_pipeline_name():
    @tool()
    class MyTool(object):
        def setup(self):
            self.name("Custom-Name")

        def get_command(self):
            return "true"

    @pipeline()
    class MyPipeline(object):
        def setup(self):
            self.name("The-Pipeline")

        def pipeline(self):
            p = Pipeline()
            p.run("MyTool")
            return p

    p = Pipeline()
    p.run('MyPipeline')
    jobs = create_jobs(p)
    assert len(jobs) == 1
    assert jobs[0].pipeline == 'The-Pipeline'
    assert jobs[0].name == 'Custom-Name'


def test_pipeline_and_tool_name_multiplex():
    @tool()
    class MyTool(object):
        """My tool
        usage:
            mytool <input>
        """
        def setup(self):
            self.name("MyName")

        def get_command(self):
            return "true"

    @pipeline()
    class MyPipeline(object):
        def setup(self):
            self.name("The-Pipeline")

        def pipeline(self):
            p = Pipeline()
            p.run("MyTool", input=['A', 'B'])
            return p
    p = Pipeline()
    p.run('MyPipeline')
    jobs = create_jobs(p, validate=False)
    assert jobs[0].pipeline == 'The-Pipeline'
    assert jobs[0].name == 'MyName.0'
    assert jobs[1].pipeline == 'The-Pipeline'
    assert jobs[1].name == 'MyName.1'


def test_dynamic_pipeline_names_fun():
    @pipeline()
    class MyPipeline(object):
        """Pipeline

        usage:
            pipeline --input <input>
        """
        def setup(self):
            self.name("${input|name}")

        def pipeline(self):
            p = Pipeline()
            ls = p.bash("ls ${input}", input=self.input)
            p.context(locals())
            return p

    p = Pipeline()
    p.run('MyPipeline', input='A')
    jobs = create_jobs(p, validate=False)
    assert len(jobs) == 1
    assert jobs[0].pipeline == 'A'
    assert jobs[0].name == 'ls'


def test_dynamic_pipeline_names_attr():
    @pipeline()
    class MyPipeline(object):
        """Pipeline

        usage:
            pipeline --input <input>
        """
        def setup(self):
            self.job.name = "${input|name}"

        def pipeline(self):
            p = Pipeline()
            ls = p.bash("ls ${input}", input=self.input)
            p.context(locals())
            return p

    p = Pipeline()
    p.run('MyPipeline', input='A')
    jobs = create_jobs(p, validate=False)
    assert len(jobs) == 1
    assert jobs[0].pipeline == 'A'
    assert jobs[0].name == 'ls'


def test_dynamic_pipeline_names_fun_multiplex():
    @pipeline()
    class MyPipeline(object):
        """Pipeline

        usage:
            pipeline --input <input>
        """
        def setup(self):
            self.name("${input|name}")

        def pipeline(self):
            p = Pipeline()
            ls = p.bash("ls ${input}", input=self.input)
            p.context(locals())
            return p

    p = Pipeline()
    p.run('MyPipeline', input=['A', 'B', 'C'])
    jobs = create_jobs(p, validate=False)
    assert len(jobs) == 3
    assert jobs[0].pipeline == 'A'
    assert jobs[0].name == 'ls.0'
    assert jobs[1].pipeline == 'B'
    assert jobs[2].name == 'ls.2'
    assert jobs[2].pipeline == 'C'
    assert jobs[1].name == 'ls.1'


def test_dynamic_pipeline_names_attr_multiplex():
    @pipeline()
    class MyPipeline(object):
        """Pipeline

        usage:
            pipeline --input <input>
        """
        def setup(self):
            self.job.name = "${input|name}"

        def pipeline(self):
            p = Pipeline()
            ls = p.bash("ls ${input}", input=self.input)
            p.context(locals())
            return p

    p = Pipeline()
    p.run('MyPipeline', input=['A', 'B', 'C'])
    jobs = create_jobs(p, validate=False)
    assert len(jobs) == 3
    assert jobs[0].pipeline == 'A'
    assert jobs[0].name == 'ls.0'
    assert jobs[1].pipeline == 'B'
    assert jobs[2].name == 'ls.2'
    assert jobs[2].pipeline == 'C'
    assert jobs[1].name == 'ls.1'
