#!/usr/bin/env python
import jip


def test_tool_name():
    @jip.tool()
    class MyTool():
        def validate(self):
            pass

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()
    jobs = jip.create_jobs(p)
    assert len(jobs) == 1
    assert jobs[0].name == "MyTool"
    assert jobs[0].pipeline is None


def test_tool_name_external_profile_on_pipeline():
    @jip.tool()
    class MyTool():
        def validate(self):
            pass

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()

    profile = jip.Profile(name="testname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 1
    assert jobs[0].name == "MyTool"
    assert jobs[0].pipeline == "testname"


def test_tool_name_set_in_validate():
    @jip.tool()
    class MyTool():
        def validate(self):
            self.name("testtool")

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()

    profile = jip.Profile(name="testname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 1
    assert jobs[0].name == "testtool"
    assert jobs[0].pipeline == "testname"


def test_tool_name_set_in_validate_with_job():
    @jip.tool()
    class MyTool():
        def validate(self):
            self.job.name = "testtool"

        def get_command(self):
            return "echo"

    p = jip.Pipeline()
    p.run('MyTool')
    p.expand()

    profile = jip.Profile(name="testname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 1
    assert jobs[0].name == "testtool"
    assert jobs[0].pipeline == "testname"


def test_tool_name_in_pipeline_context():
    @jip.tool()
    class MyTool():
        def validate(self):
            self.job.name = "testtool"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    jobs = jip.create_jobs(p)
    assert len(jobs) == 1
    assert jobs[0].name == "testtool"
    assert jobs[0].pipeline == "thepipeline"


def test_tool_name_in_pipeline_context_with_custom_profile():
    @jip.tool()
    class MyTool():
        def validate(self):
            self.job.name = "testtool"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 1
    assert jobs[0].name == "testtool"
    assert jobs[0].pipeline == "customname"


def test_tool_name_in_pipeline_context_with_custom_profile_and_custom_name():
    @jip.tool()
    class MyTool():
        def validate(self):
            self.job.name = "testtool"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.job('Tool1').run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 1
    assert jobs[0].name == "Tool1"
    assert jobs[0].pipeline == "customname"


def test_tool_name_in_pipelines_with_multiplexing():
    @jip.tool()
    class MyTool():
        """mytool
        usage:
            mytool <data>
        """
        def validate(self):
            self.job.name = "Tool"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.run('MyTool', data=["A", "B"])
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 2
    assert jobs[0].name == "Tool.0"
    assert jobs[0].pipeline == "customname"
    assert jobs[1].name == "Tool.1"
    assert jobs[1].pipeline == "customname"


def test_tool_name_in_pipelines_with_multiplexing_and_custom_name():
    @jip.tool()
    class MyTool():
        """mytool
        usage:
            mytool <data>
        """
        def validate(self):
            self.job.name = "somename"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.job("Tool").run('MyTool', data=["A", "B"])
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 2
    assert jobs[0].name == "Tool.0"
    assert jobs[0].pipeline == "customname"
    assert jobs[1].name == "Tool.1"
    assert jobs[1].pipeline == "customname"


def test_tool_name_in_pipelines_with_multiplexing_and_custom_template_name():
    @jip.tool()
    class MyTool():
        """mytool
        usage:
            mytool <data>
        """
        def validate(self):
            self.job.name = "${data}"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.run('MyTool', data=["A", "B"])
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()
    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 2
    assert jobs[0].name == "A"
    assert jobs[0].pipeline == "customname"
    assert jobs[1].name == "B"
    assert jobs[1].pipeline == "customname"


def test_tool_name_in_pipelines_with_multiplexing_and_custom_template_name_as_job():
    @jip.tool()
    class MyTool():
        """mytool
        usage:
            mytool <data>
        """
        def validate(self):
            self.job.name = "something"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def validate(self):
            self.name("thepipeline")

        def pipeline(self):
            p = jip.Pipeline()
            p.job("${data}").run('MyTool', data=["A", "B"])
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()
    profile = jip.Profile(name="customname")
    jobs = jip.create_jobs(p, profile=profile)
    assert len(jobs) == 2
    assert jobs[0].name == "A"
    assert jobs[0].configuration['data'].get() == "A"
    assert jobs[0].pipeline == "customname"
    assert jobs[1].name == "B"
    assert jobs[1].configuration['data'].get() == "B"
    assert jobs[1].pipeline == "customname"
