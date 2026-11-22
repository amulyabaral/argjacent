import os
import subprocess
import sys

import pytest


@pytest.fixture(scope="session")
def data(tmp_path_factory):
    """the four simulated self-test cases (gene.fa, truth.json, *_1.fq, *_2.fq)"""
    from argjacent import selftest
    d = tmp_path_factory.mktemp("sim")
    selftest.simulate(str(d))
    return str(d)


def run_cli(args, expect=0):
    """argjacent in a fresh interpreter (as a user would run it) -> (exit code, stdout, stderr)"""
    r = subprocess.run([sys.executable, "-m", "argjacent"] + list(args), capture_output=True, text=True)
    if expect is not None:
        assert r.returncode == expect, r.stdout + r.stderr
    return r.returncode, r.stdout, r.stderr


@pytest.fixture(scope="session")
def cli():
    return run_cli
