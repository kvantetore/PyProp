

def rewrite_extension(filename, newext):
	return "%s%s" % (path.splitext(filename)[0], newext)


def compile_cxx(conf, src, dst=None):
	if not dst:
		dst = rewrite_extension(src, ".o")

	defines = ["-D%s" % c for c in conf.defines]
	include_path = ["-I%s" % c for c in conf.include_path]

	bld.run(conf.compiler_cxx, "-c", src, "-o", dst, defines, include_path)


