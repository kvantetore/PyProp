#!/usr/bin/env python

import pyprop.build.fabricate as fab
import pyprop.build.config as build_config
import pyprop.build.rules as rules
import os.path as path
import os

bld = fab.Builder()
bld_pyste = fab.Builder(depsname=".deps-pyste")
conf = build_config.setup()

sources = []
subdirs = [
	"core",
]

def build():
	for src in sources:
		rules.compile_cxx(src)

	for dir in subdirs:
		curdir = path.abspath(path.curdir)
		os.chdir(dir)
		fab.run("./build.py")
		os.chdir(curdir)


fab.main()
