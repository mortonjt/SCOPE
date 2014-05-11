#!/usr/bin/python2.7
# wrapper.py: A wrapper for executable sequeces that can't be easily timed.
from subprocess import call
from sys import argv

call(argv[1], shell=True)
