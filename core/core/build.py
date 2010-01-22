#!/usr/bin/env python
import sys
import os
import os.path

pyprop_path = "../../"
sys.path.insert(0, os.path.abspath(pyprop_path))

import pyprop.build.fabricate as fab
import pyprop.build.config as build_config
import pyprop.build.rules as rules

bld = fab.Builder()
bld_pyste = fab.Builder(depsname=".deps-pyste")
conf = build_config.setup()

sources = [
	"representation/overlapmatrix.cpp",
]



def build():
	for s in sources:
		rules.compile_cxx(conf, bld, s)


fab.main()
