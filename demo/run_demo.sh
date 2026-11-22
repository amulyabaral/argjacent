#!/bin/sh
# Self-test on simulated data with a known answer (about 10 s).  Usage: sh demo/run_demo.sh [out_dir]
# (With the conda or pip package installed, the same test is the command  argjacent-selftest.)
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd -P)
if [ ! -f "$ROOT/.python" ]; then
    echo "Run  sh $ROOT/install.sh  first." >&2
    exit 1
fi
PY=$(cat "$ROOT/.python")
OUT=${1:-$ROOT/demo_out}
PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" exec "$PY" -m argjacent.selftest "$OUT"
