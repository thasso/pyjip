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


#def test_tool_name_set_in_validate_with_job():
    #@jip.tool()
    #class MyTool():
        #def validate(self):
            #self.job.name = "testtool"

        #def get_command(self):
            #return "echo"

    #p = jip.Pipeline()
    #p.run('MyTool')
    #p.expand()

    #profile = jip.Profile(name="testname")
    #jobs = jip.create_jobs(p, profile=profile)
    #assert len(jobs) == 1
    #assert jobs[0].name == "testtool"
    #assert jobs[0].pipeline == "testname"
