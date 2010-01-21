#!/usr/bin/env python
import sys
import os
import os.path

print "Hello from core"

pyprop_path = "../"
sys.path.insert(0, os.path.abspath(pyprop_path))

import pyprop.build.fabricate as fab
import pyprop.build.config as build_config
import pyprop.build.rules as rules

bld = fab.Builder()
bld_pyste = fab.Builder(depsname=".deps-pyste")
conf = build_config.setup()

sources = []
subdirs = [
	"core",
]

def build():
	for dir in subdirs:
		curdir = path.abspath(path.curdir)
		os.chdir(dir)
		fab.run("./build.py")
		os.chdir(curdir)


fab.main()
