# MUMPS 5.3.5 on Apple Silicon + R2025b

Built 2026-05-28 against:
- MUMPS 5.3.5 source from https://mumps-solver.org/MUMPS_5.3.5.tar.gz
- gfortran 15.2.0 (Homebrew GCC 15)
- Apple clang 21.0
- OpenBLAS (Homebrew)
- METIS 5.1.0 (Homebrew core)
- Sequential build (no MPI, no ScaLAPACK)

## Required: launch MATLAB with OMP_NUM_THREADS=1

MATLAB pre-loads Intel libiomp5 which conflicts with gfortran's libgomp when
MUMPS tries to spawn threads. Force single-threaded by setting:

    OMP_NUM_THREADS=1 /Applications/MATLAB_R2025b.app/bin/matlab

Or set it from a wrapper. Single-threaded MUMPS on a 10-core M-series chip is
still extremely fast for this problem class.

## Verified on smoke tests
- 3x3 dense system: machine-precision
- 200x200 Laplacian (40k DOFs): 0.16s
- 400x400 Laplacian (160k DOFs): 0.74s

## Single precision (smumps, added 2026-06-01)

`smumps` halves the factor/front memory (~2×), extending the in-core ceiling
by ~one mode. Built with `make s`, plus an `s` arithmetic block added to
`mumpsmex.c` (`#define MUMPS_ARITH_s 1`; map `dmumps_c→smumps_c`,
`double2→float`), built via `mex -R2017b ... -DMUMPS_ARITH=MUMPS_ARITH_s -lsmumps`.

Enable in the solver with `Opt.mumps_single = 1` (run_direct dispatches to
`@smumps` and casts the solution back to double). **Pass the matrix/RHS as
DOUBLE** — MATLAB sparse is always double and smumps casts to float internally;
`single(sparseMatrix)` is invalid and segfaults the MEX. Newton step accuracy
is ~1e-6, fine under continuation.

## Memory ceiling (16 GB laptop)

Factor/front memory grows super-linearly with the number of coupled HNS modes
(M). In-core ceiling: **M=2 double, ~M=3 single**. **M=5 (Scholten spectrum)
monolithic MUMPS is infeasible on 16 GB in any precision** — the dense modal-
coupling front exceeds RAM; BLR (ICNTL35=1/2) and OOC compress/stream the
stored factor, not the front. Use Picard for full M=5; MUMPS for high-amplitude
at M≤3. Full M=5 + high amplitude needs a workstation.

## Files
- `dmumpsmex.mexmaca64` — double MEX (1.9 MB); `dmumps.m`, `initmumps.m` wrappers
- `smumpsmex.mexmaca64` — single MEX (1.9 MB); `smumps.m` wrapper
- `mumps_help.m` — help text
