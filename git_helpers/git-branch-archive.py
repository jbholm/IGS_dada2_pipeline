#!/usr/bin/env python2

import subprocess
import sys
import argparse

parser=argparse.ArgumentParser(
    description='''A script that automates the archival of a git branch by tagging it, deleting it, and pushing the changes to origin. ''',
    epilog="""""")
parser.add_argument('branch', nargs=1, help='The branch to archive')
args=parser.parse_args()

print "$ git tag " + sys.argv[1] + " " + sys.argv[1]
status = subprocess.call(["git", "tag", sys.argv[1], sys.argv[1]])
if (status != 0):
    print "Exit status:", status
    exit()
else:
    print "$ git branch -D " + sys.argv[1]
    subprocess.call(["git", "branch", "-D", sys.argv[1]])
    print "$ git push origin :" + sys.argv[1]
    subprocess.call(["git", "push", "origin", ":" + sys.argv[1]])
    print "$ git push origin tag " + sys.argv[1]
    subprocess.call(["git", "push", "origin", "tag", sys.argv[1]])
