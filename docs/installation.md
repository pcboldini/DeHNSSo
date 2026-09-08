# Installation

Everything you need before running DeHNSSo, in the order you need it.

## Contents

- [1. Prerequisites](#1-prerequisites)
- [2. The `lapack_band` MEX](#2-the-lapack_band-mex)
- [3. MUMPS, for the Newton solver](#3-mumps-for-the-newton-solver)

---

## 1. Prerequisites

- MATLAB R2020b or later. No additional MATLAB toolboxes are used.
- A base-flow input file `DeHNSSo_input_<case>.mat` holding the `StabGrid` struct.
  Generate one with the scripts in `gridgen/benchmark/`.

Two optional pieces are built once per machine, and only if you need them:

| | needed for | section |
|---|---|---|
| `lapack_band` MEX | `Opt.lu_mode = 'lapack_band'`, the fastest back-end for the compressible operator | [2](#2-the-lapack_band-mex) |
| MUMPS | `Opt.solver = 'newton'`. Nothing else uses it | [3](#3-mumps-for-the-newton-solver) |

Picard runs need neither.

---

## 2. The `lapack_band` MEX

`lu_mode = 'lapack_band'` calls a MEX wrapper around LAPACK's banded LU
(`zgbsv`/`dgbsv`). The MEX binary is platform-specific and is **not** shipped
with the repo — each user builds it once.

```matlab
>> run('src/linsolve/lapack_band/mexbuild_dgbsv.m')
```

Requirements:

- MATLAB R2018a or later (uses the `-R2018a` interleaved-complex API).
- A C compiler configured for MEX (`mex -setup` once if you've never used MEX).
  - macOS: Xcode Command Line Tools.
  - Linux: `gcc`.
  - Windows: MSVC / MinGW.

If `lu_mode = 'lapack_band'` is selected but the MEX is missing, `solve_modes`
errors with a message pointing to `mexbuild_dgbsv.m`. Fall back to
`lu_mode = 'auto'` if you cannot build (e.g. no compiler available).

---

## 3. MUMPS, for the Newton solver

Newton factorizes the Jacobian with MUMPS. It is a separate dependency, and it is the only
part of DeHNSSo that needs a Fortran toolchain.

### macOS (Apple Silicon)

Run the installer:

```bash
src/external/mumps_5.3.5/install/install_mumps.sh          # ~/mumps_build by default
```

It downloads MUMPS 5.3.5, builds the library and the MEX, checks what it linked, and solves
a test system. At the end it prints the two lines to add to your caller. If any step fails it
stops, so you never get a build it did not check. You need `brew install gcc metis` and the
Xcode command line tools.

To build it by hand instead, use the two files the installer uses:

| | |
|---|---|
| `install/Makefile.inc.macos-arm64` | the library: copy to `MUMPS_5.3.5/Makefile.inc`, then `make d` |
| `install/make.inc.matlab.macos-arm64` | the MEX: copy to `MUMPS_5.3.5/MATLAB/make.inc`, then `make d` there |

Four settings break the build if you get them wrong, and the comments in those files point
each one out: `AR` needs a trailing space, `OPTF` needs `-fallow-argument-mismatch`, no file
may use `-fopenmp`, and the MEX needs `-R2017b` and not `-largeArrayDims`. When it is built,
run `otool -L dmumpsmex.mexmaca64`. It should show Accelerate and should not show OpenBLAS.

Linux support is under development. The settings above are for macOS only and were not
tested on Linux, and there is no installer for it yet.

Once it is built, see `docs/DeHNSSo_guide.md`, "Running it: MUMPS, memory and where", for how
to size a run, put the build on the path, and check the answer.
