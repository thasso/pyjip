#!/usr/bin/env python
"""Helper script to show the change log based on the github milestones"""

import requests

repo_url = "https://api.github.com/repos/thasso/pyjip"

milestones = requests.get(repo_url + "/milestones").json()

for milestone in milestones:
    mid = milestone['number']
    issues = requests.get(repo_url + "/issues", params={
        'milestone': mid,
        'state': 'closed'
    }).json()
    print "Version %s:" % milestone['title']
    for issue in issues:
        print "    * %s [`Issue #%s <%s>`_]" % \
            (issue['title'], issue['number'], issue['html_url'])
