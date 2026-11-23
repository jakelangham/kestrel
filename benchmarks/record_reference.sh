#!/bin/bash
# Record the benchmark errors of a commit (default: the head of main) in
# references/<short hash>.txt, for runbench.jl --reference, which calls this
# when that file is missing.
#
# usage: ./record_reference.sh [git ref]      (run from this directory)
#
# The commit is built, unchanged (./configure --enable-openmp), in a temporary
# git worktree, and runbench.jl runs every case there at its default levels.
# The benchmark harness (runbench.jl, benchlib.jl, the Input files and
# ../tests/testlib.jl) is copied from this checkout, so only the solver comes
# from the commit.
#
# The inputs are frictionless (Drag = None). A commit without Drag = None (such
# as main before the fix-timestep-rejection branch was merged) runs them with
# Drag = Chezy and Chezy Co = 1e-8 instead, the smallest drag it accepts. That
# changes the errors by up to about 7% from frictionless flow, and the
# reference file says so.
#
# Set OMP_NUM_THREADS to choose the number of threads (default here: 4). The
# results do not depend on it.
set -e
cd "$(dirname "$0")"
REF=${1:-main}
R=$(git rev-parse --show-toplevel)
HASH=$(git rev-parse --short --verify "$REF^{commit}")
OUT="references/$HASH.txt"
W=$(mktemp -d)
trap 'git -C "$R" worktree remove --force "$W/tree" >/dev/null 2>&1; rm -rf "$W"' EXIT

echo "Building $REF ($HASH) in $W ..."
git -C "$R" worktree add --detach "$W/tree" "$HASH" >/dev/null 2>&1
T="$W/tree"
( cd "$T" && autoreconf -fi && ./configure --enable-openmp && make ) > "$W/build.log" 2>&1 ||
   { echo "Build failed; the end of the log:"; tail -20 "$W/build.log"; exit 1; }

mkdir -p "$T/benchmarks"
cp runbench.jl benchlib.jl Input_*.txt "$T/benchmarks/"
cp "$R/tests/testlib.jl" "$T/tests/"
NOTE=""
if ! git show "$HASH:src/Parameters.f90" | grep -q "DragClosure => NoDrag"; then
   sed -i 's/^Drag = None$/Drag = Chezy\nChezy Co = 1e-8/' "$T"/benchmarks/Input_*.txt
   NOTE="# Note: run with Drag = Chezy and Chezy Co = 1e-8 (this commit has no Drag = None)"
fi

echo "Running the benchmarks at $HASH ..."
( cd "$T/benchmarks" &&
  OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} julia runbench.jl --reference none --save "$W/errors.txt" ) |
   tr '\r' '\n' | grep -v '^[/|\\-]*$' || true
[ -f "$W/errors.txt" ] || { echo "No errors were recorded."; exit 1; }

mkdir -p references
{ head -2 "$W/errors.txt"
  [ -z "$NOTE" ] || echo "$NOTE"
  tail -n +3 "$W/errors.txt"; } > "$OUT"
echo "Wrote $OUT."
