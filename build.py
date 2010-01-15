#!/usr/bin/env python

import pyprop.build.fabricate as fab
import pyprop.build.config as build_config
import pyprop.build.rules as rules
import os.path as path

bld = fab.Builder()
bld_pyste = fab.Builder(depsname=".deps-pyste")
conf = build_config.setup()

subdirs = [
	"core",
]

def build():
	print "Hello"
	for src in sources:
		rules.compile_cxx(src)

	for dir in subdirs:
		fab.run("%s/build.py" % dir)
	

fab.main()
