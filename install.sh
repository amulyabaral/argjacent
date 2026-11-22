#!/bin/sh
# One-time setup for use from this folder (without conda or pip): builds the C programs aw, awx and awv into argjacent/bin/
# and makes sure Python has numpy.
# Usage:  sh install.sh          (set PYTHON=/path/to/python3 to choose a Python)
set -e
cd "$(dirname "$0")"
ROOT=$(pwd -P)
echo "== ARGjacent: install in $ROOT"

# 1. C compiler (macOS: Xcode Command Line Tools)
if ! cc --version >/dev/null 2>&1; then
    echo "No working C compiler (cc)."
    echo "On macOS run:   xcode-select --install"
    echo "(a window may already have opened for this), then run:  sh install.sh"
    exit 1
fi

# 2. build
mkdir -p argjacent/bin
echo "-- building argjacent/bin/aw, awx, awv"
if ! { cc -O3 -o argjacent/bin/aw src/aw.c -lz -lpthread -lm && cc -O3 -o argjacent/bin/awx src/awx.c -lpthread && cc -O3 -o argjacent/bin/awv src/awv.c -lpthread; }; then
    echo "Build failed. Please send the messages above."
    exit 1
fi
argjacent/bin/aw -v

# 3. Python 3.8+ with numpy
PY=${PYTHON:-python3}
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "python3 not found. On macOS it comes with the Command Line Tools (xcode-select --install)."
    exit 1
fi
if ! "$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)'; then
    echo "Python 3.8 or newer is needed; found: $("$PY" --version 2>&1)"
    exit 1
fi
if "$PY" -c 'import numpy' >/dev/null 2>&1; then
    PYBIN=$("$PY" -c 'import sys; print(sys.executable)')
    echo "-- using $PYBIN (numpy already installed)"
else
    echo "-- setting up numpy in $ROOT/.venv (needs internet once, about 30 s)"
    if [ ! -x .venv/bin/python ]; then
        "$PY" -m venv .venv
    fi
    .venv/bin/python -m pip install --quiet --disable-pip-version-check --upgrade pip
    .venv/bin/python -m pip install --quiet --disable-pip-version-check "numpy>=1.20"
    PYBIN=$ROOT/.venv/bin/python
fi
"$PYBIN" -c 'import numpy; print("-- numpy", numpy.__version__)'
echo "$PYBIN" > .python

# 4. smoke test
chmod +x bin/argjacent bin/argjacent-cohort
./bin/argjacent --version >/dev/null
./bin/argjacent-cohort --version >/dev/null
echo
echo "Installed. Next:"
echo "  sh demo/run_demo.sh                  # 10-second self-test on simulated data"
echo "  ./bin/argjacent -g genes.fa -1 R1.fq.gz -2 R2.fq.gz -o results/sample"
echo "  ./bin/argjacent-cohort -g genes.fa -i results/s1.flanks.json results/s2.flanks.json -o results/cohort"
