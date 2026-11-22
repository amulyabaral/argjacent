"""Builds the C programs `aw`, `awx` and `awv` into the package (argjacent/bin/); everything else is in pyproject.toml.
The compiler and its flags come from the environment (CC, CFLAGS, LDFLAGS), as conda-build sets them; aw needs zlib."""
import os
import shlex
import subprocess

import setuptools
from setuptools import Distribution, setup
from setuptools.command.build_py import build_py

if int(setuptools.__version__.split(".")[0]) < 61:
    # older setuptools ignores pyproject.toml and would install an empty package called UNKNOWN (seen with the pip of
    # Ubuntu 22.04)
    raise SystemExit("argjacent needs setuptools 61 or newer to build. Run:  python -m pip install -U pip setuptools")


class BuildWithAw(build_py):
    def run(self):
        super().run()
        out_dir = os.path.join(self.build_lib, "argjacent", "bin")
        os.makedirs(out_dir, exist_ok=True)
        for name, libs in (("aw", ["-lz", "-lpthread", "-lm"]), ("awx", ["-lpthread"]), ("awv", ["-lpthread"])):
            cmd = (shlex.split(os.environ.get("CC", "cc")) + ["-O3"] + shlex.split(os.environ.get("CPPFLAGS", "")) +
                   shlex.split(os.environ.get("CFLAGS", "")) + ["-o", os.path.join(out_dir, name), os.path.join("src", name + ".c")] +
                   shlex.split(os.environ.get("LDFLAGS", "")) + libs)
            self.announce("building %s: %s" % (name, " ".join(cmd)), level=2)
            subprocess.check_call(cmd)


class BinaryDistribution(Distribution):
    """the wheel holds compiled programs: it is specific to the platform"""

    def has_ext_modules(self):
        return True


setup(cmdclass={"build_py": BuildWithAw}, distclass=BinaryDistribution)
