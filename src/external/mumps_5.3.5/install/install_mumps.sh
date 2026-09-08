#!/usr/bin/env bash
# Build MUMPS 5.3.5 and its MATLAB MEX for DeHNSSo (macOS / Apple Silicon).
#
#   ./install_mumps.sh [build_dir] [matlab_root]
#
# Defaults: ~/mumps_build and /Applications/MATLAB_R2025b.app
# Needs: Homebrew gcc (gfortran) and metis, Xcode command line tools, MATLAB.
#   brew install gcc metis
#
# Only the double-real arithmetic is built, which is all DeHNSSo uses.
set -euo pipefail

BUILD_DIR="${1:-$HOME/mumps_build}"
MATLAB_ROOT="${2:-/Applications/MATLAB_R2025b.app}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER=5.3.5
SRC="$BUILD_DIR/MUMPS_$VER"

echo "==> build dir : $BUILD_DIR"
echo "==> MATLAB    : $MATLAB_ROOT"

for t in gfortran clang make curl; do
    command -v $t >/dev/null || { echo "ERROR: $t not found."; exit 1; }
done
[ -x "$MATLAB_ROOT/bin/mex" ] || { echo "ERROR: no mex at $MATLAB_ROOT/bin/mex"; exit 1; }

mkdir -p "$BUILD_DIR"
if [ ! -d "$SRC" ]; then
    echo "==> downloading MUMPS $VER (public, no registration)"
    curl -fL -o "$BUILD_DIR/MUMPS_$VER.tar.gz" "https://mumps-solver.org/MUMPS_$VER.tar.gz"
    tar xf "$BUILD_DIR/MUMPS_$VER.tar.gz" -C "$BUILD_DIR"
fi

echo "==> installing Makefile.inc"
cp "$HERE/Makefile.inc.macos-arm64" "$SRC/Makefile.inc"

echo "==> building the library (make d)"
( cd "$SRC" && make d )

echo "==> installing MATLAB/make.inc"
sed -e "s#\$(HOME)/mumps_build/MUMPS_$VER#$SRC#" \
    -e "s#/Applications/MATLAB_R2025b.app#$MATLAB_ROOT#" \
    "$HERE/make.inc.matlab.macos-arm64" > "$SRC/MATLAB/make.inc"

echo "==> building the MEX (d only)"
( cd "$SRC/MATLAB" && rm -f dmumpsmex.mexmaca64 dmumpsmex.stamp && make d )

MEXFILE="$SRC/MATLAB/dmumpsmex.mexmaca64"
[ -f "$MEXFILE" ] || { echo "ERROR: mex was not produced."; exit 1; }

echo "==> checking what it linked"
if otool -L "$MEXFILE" | grep -qi openblas; then
    echo "ERROR: linked against OpenBLAS; expected Accelerate."; exit 1
fi
otool -L "$MEXFILE" | grep -qi Accelerate || { echo "ERROR: Accelerate not linked."; exit 1; }
echo "    Accelerate OK, no OpenBLAS"

echo "==> solving a test system in MATLAB"
# matlab -batch takes ONE line, so the test goes in a file
TESTM="$SRC/MATLAB/dehnsso_mumps_selftest.m"
cat > "$TESTM" <<'TESTEOF'
n = 200000;  e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);   b = A*ones(n,1);
id = initmumps;  id.SYM = 0;  id = dmumps(id);
id.JOB = 6;  id.RHS = b;
t = tic;  id = dmumps(id, A);  dt = toc(t);
err = max(abs(id.SOL - 1));
id.JOB = -2;  id = dmumps(id);          % capture the output or the instance leaks
fprintf('    solved %dx%d in %.2f s, max err %.2e\n', n, n, dt, err);
if err > 1e-6
    error('MUMPS self-test: solution is wrong (err %.2e)', err);
end
TESTEOF
OMP_NUM_THREADS=1 "$MATLAB_ROOT/bin/matlab" -batch \
    "addpath('$SRC/MATLAB'); run('$TESTM')" || { echo "ERROR: test solve failed."; exit 1; }
rm -f "$TESTM"

cat <<MSG

Done. Add this to your caller, AFTER any genpath of src/ (genpath would
otherwise shadow it with src/external/mumps_5.3.5/dmumps.m):

    addpath('$SRC/MATLAB', '-begin');
    fprintf('dmumps resolves to: %s\n', which('dmumps'));

and always start MATLAB with OMP_NUM_THREADS=1:

    OMP_NUM_THREADS=1 matlab -batch "maxNumCompThreads(1); my_caller"
MSG
