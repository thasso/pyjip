#!/usr/bin/env python
import jip
import pytest


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


@pytest.mark.parametrize('data', [1, 3])
def test_pipeline_overwrites_tool(data):
    @jip.tool()
    class MyTool():
        def setup(self):
            self.profile.threads = 2

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job(threads=data).run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == data


def test_pipeline_overwrites_pipeline_from_spec():
    @jip.tool()
    class MyTool():
        def setup(self):
            self.profile.threads = 2
            self.profile.queue = "Org"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job(threads=3, queue="Yeah").run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(threads=10, queue="Test")
    profile.specs['MyTool'] = jip.Profile(threads=5)
    profile.apply_to_pipeline(p)

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == 5
    assert jobs[0].queue == "Yeah"


def test_pipeline_tool_defaults():
    @jip.tool()
    class MyTool():
        def setup(self):
            self.profile.threads = 2
            self.profile.queue = "Org"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job().run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile()
    profile.specs['MyTool'] = jip.Profile()
    profile.apply_to_pipeline(p)

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == 2
    assert jobs[0].queue == "Org"


def test_pipeline_tool_defaults_global():
    @jip.tool()
    class MyTool():
        def setup(self):
            self.profile.threads = 2
            self.profile.queue = "Org"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job().run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(threads=5, queue="yeah")
    profile.specs['MyTool'] = jip.Profile()
    profile.apply_to_pipeline(p)

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == 2
    assert jobs[0].queue == "Org"


def test_pipeline_tool_defaults_global_job():
    @jip.tool()
    class MyTool():
        def setup(self):
            self.profile.threads = 2
            self.profile.queue = "Org"

        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job(threads=3, queue="Intern").run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(threads=5, queue="yeah")
    profile.specs['MyTool'] = jip.Profile()
    profile.apply_to_pipeline(p)

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == 3
    assert jobs[0].queue == "Intern"


def test_pipeline_tool_spec_regexp():
    @jip.tool()
    class MyTool():
        def get_command(self):
            return "echo"

    @jip.pipeline()
    class MyPipeline():
        def pipeline(self):
            p = jip.Pipeline()
            p.job(threads=3, queue="Intern").run('MyTool')
            return p

    p = jip.Pipeline()
    p.run('MyPipeline')
    p.expand()

    profile = jip.Profile(threads=5, queue="yeah", priority="high")
    profile.specs['My*'] = jip.Profile(threads=10, queue="rock")
    profile.apply_to_pipeline(p)

    jobs = jip.create_jobs(p)
    assert jobs[0].threads == 10
    assert jobs[0].queue == "rock"
    assert jobs[0].priority == "high"
