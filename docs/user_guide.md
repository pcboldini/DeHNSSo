# DeHNSSo — User Guide

**DeHNSSo-v2** | Delft Harmonic Navier-Stokes Solver

## Contents

- [1. What DeHNSSo does](#1-what-dehnsso-does)
- [2. Prerequisites](#2-prerequisites)
- [3. Quick start](#3-quick-start)
    - [Step 1 — Generate base flow and stability grid](#step-1--generate-base-flow-and-stability-grid)
    - [Step 1b — Optional streamwise clustering (all cases)](#step-1b--optional-streamwise-clustering-all-cases)
    - [Step 1c — Cartesian / OpenFOAM CSV input](#step-1c--cartesian--openfoam-csv-input)
    - [Step 1d — Base-flow derivatives (computed or supplied)](#step-1d--base-flow-derivatives-computed-or-supplied)
    - [Step 2 — Configure and run](#step-2--configure-and-run)
- [4. Spectral content](#4-spectral-content)
    - [Single oblique mode](#single-oblique-mode)
- [5. Inflow conditions](#5-inflow-conditions)
    - [Setting initial amplitudes](#setting-initial-amplitudes)
    - [Amplitude convention](#amplitude-convention)
    - [Restarting from a previous solution (`Stab.IC = 'LOAD'`)](#restarting-from-a-previous-solution-stabic--load)
    - [User-defined inflow profiles (`Stab.IC = 'LOAD'` with `phiIC`)](#user-defined-inflow-profiles-stabic--load-with-phiic)
    - [No inflow disturbance (`Stab.IC = 'ZERO'`)](#no-inflow-disturbance-stabic--zero)
- [6. Boundary conditions](#6-boundary-conditions)
    - [Arrays](#arrays)
    - [Wall actuation (blowing/suction)](#wall-actuation-blowingsuction)
    - [Helmholtz-resonator wall (impedance BC)](#helmholtz-resonator-wall-impedance-bc)
    - [Freestream forcing](#freestream-forcing)
    - [Freestream BC type](#freestream-bc-type)
- [7. Linear vs nonlinear mode](#7-linear-vs-nonlinear-mode)
- [8. Solver options](#8-solver-options)
    - [Buffer zone](#buffer-zone)
    - [Convergence](#convergence)
    - [Picard acceleration and under-relaxation](#picard-acceleration-and-under-relaxation)
    - [Pseudo-arclength continuation (`solve_newton`)](#pseudo-arclength-continuation-solvenewton)
    - [LU mode](#lu-mode)
    - [Fill-reducing ordering](#fill-reducing-ordering)
    - [Building the `lapack_band` MEX (one-time setup per machine)](#building-the-lapackband-mex-one-time-setup-per-machine)
    - [Parallel acceleration](#parallel-acceleration)
    - [Output control](#output-control)
    - [Iteration history](#iteration-history)
- [9. Running and reruns](#9-running-and-reruns)
- [10. Output](#10-output)
    - [Saving results to disk](#saving-results-to-disk)
    - [Figures produced during a run](#figures-produced-during-a-run)
    - [Post-processing](#post-processing)
- [11. Common issues](#11-common-issues)

## Option index

Every `Stab` and `Opt` field the guide documents, grouped by what it controls. "Default" is the value used when the field is not set; **required** means the run will not start without it.

**Case definition**

| field | what it does | default | see |
|---|---|---|---|
| `Stab.M` | number of temporal modes (truncation in omega) | required | [Spectral content](#4-spectral-content) |
| `Stab.N` | number of spanwise modes; 0 is 2D | required | [Spectral content](#4-spectral-content) |
| `Stab.omega_0` | fundamental angular frequency | required | [Spectral content](#4-spectral-content) |
| `Stab.beta_0` | fundamental spanwise wavenumber | required | [Spectral content](#4-spectral-content) |
| `Opt.compressible` | `on` selects the five-variable (u,v,w,rho,T) operator | `'off'` | [Spectral content](#4-spectral-content) |
| `Opt.linear` | `on` single linear solve, `off` nonlinear iteration | `'off'` | [Linear vs nonlinear](#7-linear-vs-nonlinear-mode) |

**Inflow**

| field | what it does | default | see |
|---|---|---|---|
| `Stab.IC` | inflow condition: `ILST`, `CLST`, `PSE`, `LOAD` or `ZERO` | required | [Inflow conditions](#5-inflow-conditions) |
| `Stab.phaseRef` | phase reference: `pwall` (wall pressure) or `umax` | required | [Inflow conditions](#5-inflow-conditions) |
| `Stab.A0_fund` | fundamental amplitude, the physical peak 2·abs(u_hat) | required | [Amplitude convention](#amplitude-convention) |
| `Stab.A0` | per-mode amplitude vector; use instead of `A0_fund` to seed several modes | from `A0_fund` | [Setting initial amplitudes](#setting-initial-amplitudes) |
| `Stab.ICfile` | file read when `Stab.IC='LOAD'`: a checkpoint, or a `phiIC` profile file | — | [Restarting](#restarting-from-a-previous-solution-stabic--load) |

**Boundary conditions**

| field | what it does | default | see |
|---|---|---|---|
| `Opt.bc_top` | freestream BC type per channel; the channels differ incompressible vs compressible | all `H_DR` | [Freestream BC type](#freestream-bc-type) |
| `Stab.bcx` | streamwise coordinate the `bcw`/`bct` arrays are given on | — | [Wall actuation](#wall-actuation-blowingsuction) |
| `Stab.bcw` | wall BC values, per channel and mode | zeros (no-slip) | [Wall actuation](#wall-actuation-blowingsuction) |
| `Stab.bct` | freestream values, per channel and mode; used by the `IH_DR` type | zeros | [Freestream forcing](#freestream-forcing) |
| `Opt.hr` | Helmholtz-resonator wall: admittance, mouth mask, taper | — | [Resonator wall](#helmholtz-resonator-wall-impedance-bc) |

**Nonlinear iteration**

| field | what it does | default | see |
|---|---|---|---|
| `Opt.solver` | nonlinear solver: `picard` or `newton` | `'picard'` | [Picard acceleration](#picard-acceleration-and-under-relaxation) |
| `Opt.anderson` | Anderson acceleration depth (0 is plain Picard) | `0` | [Picard acceleration](#picard-acceleration-and-under-relaxation) |
| `Opt.relax` | Picard under-relaxation factor; use with `anderson = 0` | off | [Picard acceleration](#picard-acceleration-and-under-relaxation) |
| `Opt.Conv` | nonlinear convergence criterion | `1e-4` | [Convergence](#convergence) |
| `Opt.ConvF` | convergence criterion on the forcing | `100` | [Convergence](#convergence) |
| `Opt.TH` | nonlinear-term activation threshold | `1e-11` | [Convergence](#convergence) |
| `Opt.AFg` | amplitude-factor ramping increment | — | [Convergence](#convergence) |
| `Opt.AMAX` | maximum amplitude factor | `0.1` | [Convergence](#convergence) |
| `Opt.buffer` | outflow treatment: `on` damping, `para` parabolisation, `off` | `'on'` | [Buffer zone](#buffer-zone) |
| `Opt.xb` | buffer start, in % of the domain | `85` | [Buffer zone](#buffer-zone) |

**Newton pseudo-arclength continuation**

| field | what it does | default | see |
|---|---|---|---|
| `Opt.newton_target` | target `A0_fund` for the continuation | required | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_start` | continuation start amplitude | `Stab.A0_fund` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_seed2` | second seed point | `1.04*start` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_reuse` | steps a factorisation is reused before a forced refactor | `8` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_it_thresh` | refactor next step if the corrector exceeded this | `8` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_maxit` | corrector back-solves before the factor is called stale | `25` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_dsmax` | maximum arclength step | `3` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_tol` | corrector residual tolerance | `1e-6` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_corrector` | `modified`, or `full` to refactorise every corrector iteration | `'modified'` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_predictor` | `curvature`, or `tangent` for the first-order Keller step | `'curvature'` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_arclength` | arclength metric | `'absolute'` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_seedfile` | cache file for the two seed points | `[]` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_ckfile` | checkpoint written each continuation step | `[]` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_resume` | `on` resumes from `newton_ckfile` | `'off'` | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |
| `Opt.newton_mumps_relax` | MUMPS working-space relaxation schedule | — | [Arclength continuation](#pseudo-arclength-continuation-solvenewton) |

**Linear algebra**

| field | what it does | default | see |
|---|---|---|---|
| `Opt.lu_mode` | back-end: `auto`, `full`, `seq`, `lapack_band`, `block_thomas`, `none` | `'auto'` | [LU mode](#lu-mode) |
| `Opt.lu_max_cache` | maximum cached factorisations for the bounded cache | `max(8, nphys)` | [LU mode](#lu-mode) |
| `Opt.lu_cache` | `off` forces one-shot factor+solve | on | [LU mode](#lu-mode) |
| `Opt.lu_ordering` | fill-reducing ordering: `auto`, `symrcm` or `colamd` | `'auto'` | [Fill-reducing ordering](#fill-reducing-ordering) |
| `Opt.bt_np` | block-Thomas partition count | `2` | [Fill-reducing ordering](#fill-reducing-ordering) |
| `Opt.bt_ooc` | `on` streams the block-Thomas factor to disk | `'off'` | [Fill-reducing ordering](#fill-reducing-ordering) |
| `Opt.bt_apply` | block-Thomas apply mode | `'lu'` | [Fill-reducing ordering](#fill-reducing-ordering) |
| `Opt.parfor` | `on` distributes the per-mode back-solve | `'off'` | [Parallel acceleration](#parallel-acceleration) |

**Run and output**

| field | what it does | default | see |
|---|---|---|---|
| `Opt.rerun` | `on` reuses the assembled operator from a previous run | `'off'` | [Running and reruns](#9-running-and-reruns) |
| `Opt.plot` | `on` / `off` run-time figures | `'on'` | [Output control](#output-control) |
| `Opt.iter_history` | `on` records amplitudes, residual and energy each iteration | `'off'` | [Iteration history](#iteration-history) |
| `Opt.xloc` | station for the local energy measure; unset skips it | `[]` | [Iteration history](#iteration-history) |
| `Opt.save_output` | results file written after post-processing | `[]` | [Saving results to disk](#saving-results-to-disk) |
| `Opt.save_checkpoint` | pre-post-process state, the file `Stab.IC='LOAD'` consumes | `[]` | [Restarting](#restarting-from-a-previous-solution-stabic--load) |
---

## 1. What DeHNSSo does

DeHNSSo solves the nonlinear stability problem for laminar boundary layers. It computes the amplitude evolution of wave-like disturbances (Tollmien-Schlichting waves, cross-flow vortices) as they propagate downstream and interact nonlinearly.

The equations are the compressible Navier-Stokes equations, written for perturbations on top of a steady base flow. The perturbations are expanded in a Fourier series in time and spanwise direction (harmonic balancing), and the resulting coupled system is solved iteratively until convergence. The incompressible four-variable form `(u,v,w,p)` is solved by default; `Opt.compressible = 'on'` selects the five-variable form `(u,v,w,rho,T)`, which reduces to the incompressible one in the low-Mach limit and has been validated against published stability data up to Mach 4.5.

---

## 2. Prerequisites

- MATLAB R2020b or later (no additional toolboxes required)
- A base-flow input file `DeHNSSo_input_<case>.mat` containing the `StabGrid` struct

---

## 3. Quick start

### Step 1 — Generate base flow and stability grid

For a Blasius flat-plate boundary layer:

```matlab
cd baseflow/IBL/
caller        % produces bf_blasius.mat in baseflow/output/
```

Then build the stability grid via the per-case script in `gridgen/benchmark/`:

```matlab
cd gridgen/benchmark/
blasius       % reads bf_blasius.mat from baseflow/output/benchmark/, writes DeHNSSo_input_blasius.mat
```

For other supported cases (`blasius_csv`, `sweptwing_flat`, `sweptwing_hump`, `m3j` / OpenFOAM): see [README.md](../README.md).

### Step 1b — Optional streamwise clustering (all cases)

By default the gridgen scripts distribute `nx` stations uniformly in ξ (wall arc length for curvilinear cases, `x` for the flat plate). For cases with a localized feature (hump, leading-edge curvature, roughness element) — or simply to concentrate stations where the disturbance grows fastest — you can cluster stations around one or more centers via a Gaussian step-size modulation. Configured in the per-case gridgen script via three `params` fields:

| Parameter | Meaning |
|---|---|
| `params.mug` | Clustering center(s) in physical x (scalar or vector for multiple centers). Taken from `params.forcing.xc` when a forcing is defined, so set it explicitly only for a feature that is not a wall source, such as a hump |
| `params.spread` | How far the refinement reaches: half-width in physical x. Defaults to the `params.forcing` half-width. Larger means gentler stretching |
| `params.refine` | How many times finer Δξ is at the peak than on an equidistant mesh |

All three are in physical units, so they read the same way as the feature they follow. Internally `spread` becomes the Gaussian width `sig = 3·spread / half-domain` and `refine` becomes `ag = 1/refine`; `params.sig` and `params.ag` are **not** accepted as inputs and raise an error naming the replacement.

These fields work for both the curvilinear (`resample_curvilinear.m`) and flat-plate (`resample_cartesian.m`) paths, which share one implementation in `gridgen/src/resample/cluster_xi.m`.

Algorithm (`gridgen/src/resample/cluster_xi.m`):

1. Map each `mug` into normalized ξ ∈ [-1, 1].
2. At every candidate ξ, compute `G_k(ξ) = exp(−½·((ξ − μ_n)/σ)²)` for each hump `k`.
3. Local modulation: `fg(ξ) = 1 − (1 − ag) · max_k G_k(ξ)` → `fg = ag` at each peak.
4. Build the new ξ grid by integrating `∫ dξ / fg` and renormalizing — denser where `fg` is small.

Example (heritage `sweptwing_hump.m`):

```matlab
x_m = 0.1837410879 / 2.1394e-4;   % physical hump center / lref ≈ 859
params.mug    = [x_m];            % single hump (no forcing here, so set the centre)
params.spread = 55.137096;        % half-width the refinement reaches
params.refine = 2;                % dxi twice as fine at the peak
```

Multi-hump example:

```matlab
params.mug    = [x_m1, x_m2];     % two humps
params.spread = [55.1, 41.4];     % per-hump half-widths
params.refine = [2, 1.43];        % per-hump strengths
```

Gridgen reports what the mesh actually does, not just what was asked for: the finest Δξ as a multiple of the equidistant spacing, and the stretching, the worst ratio between neighbouring Δξ (1.00 means none). The measured refinement is always milder than the nominal `refine`, because compressing the mesh at the peak stretches it elsewhere. Aim for stretching near 1.0; widen `spread` or lower `refine` if it grows.

Leave `params.mug = []` and define no `params.forcing`, or set `params.refine = 1`, to fall back to uniform spacing (wall arc length for curvilinear, `x` for the flat plate).

### Step 1c — Cartesian / OpenFOAM CSV input

The cartesian path (`gridgen/src/resample/resample_cartesian.m`) reads either a heritage `.mat` (e.g. `bf_blasius.mat`) or a flat-lattice CSV from a CFD code (e.g. OpenFOAM `midPlane_*.csv`). For CSV input a few extra `params` fields are usually needed in the gridgen script:

| Parameter | Meaning |
|---|---|
| `params.rescale` | `true` if the CSV is dimensional (m, m/s). When set, X/Y are divided by `lref` and U/V/W by `Uref` before resampling |
| `params.Uref`, `params.lref` | Required when `rescale = true`. Together with kinematic viscosity ν they fix `params.Re = Uref·lref/ν` |
| `params.interp_method` | `interp2` method used on the rescaled grid. Default `'makima'` — Modified Akima, handles non-uniform wall-clustered CFD spacing without overshoot. Alternatives: `'spline'` (smoother interior but boundary ringing), `'linear'` (safest, lowest order), `'cubic'` (silently falls back to spline on non-uniform grids) |
| `params.xi_trim_inflow` | Fraction (0–1) of the x-span to drop at the inflow before resampling. Useful for CFD sources with noisy inflow BCs |
| `params.xi_trim_outflow` | Same at the outflow |

The path detects whether the CSV is structured (`X`, `Y` are 2-D arrays) or a flat lattice scatter (`X`, `Y` are paired column vectors of length `nx·ny`). The reshape from flat scatter to canonical `(ny × nx)` form is automatic — no manual preprocessing needed.

Example (heritage `blasius_csv.m`):

```matlab
params.rescale      = true;
params.Uref         = 20.9;                    % m/s
params.lref         = 2.688803650654011e-04;   % m
params.Re           = 3.719126161394363e+02;   % Uref·lref/ν_air
% params.interp_method   = 'makima';   % default; uncomment to override
% params.xi_trim_inflow  = 0.02;       % drop first 2% of x
% params.xi_trim_outflow = 0.0;
```

Leaving `params.n_eta_new` / `params.n_xi_new` empty keeps the source resolution. Setting them resamples to the requested size.

### Step 1d — Base-flow derivatives (computed or supplied)

`gridgen` computes the physical base-flow derivatives that the stability operator needs and stores them on the grid, so the solver never differentiates the base flow itself:

- First derivatives `∂/∂x`, `∂/∂y` of `U, V, W` (and of `rho, T` for a compressible base flow).
- Compressible only: second streamwise `∂²/∂x²` (`dxxU/V/W/T`), second wall-normal `∂²/∂y²` (`ddyU/V/W/T`), and mixed `∂²/∂x∂y` (`dxyU`, `dxyV`).

By default these come from finite differences of the (interpolated) base flow. If you have an **analytic** base flow you can instead supply the derivatives directly on the BF struct — gridgen detects them automatically and skips the corresponding FD differentiation. No flag is needed:

```matlab
BF.dxU = ...;  BF.dyU = ...;   % supply both halves of a pair
```

Rules:

- Detection is **per field** and **per pair**: supply *both* `dx<f>` and `dy<f>` to replace that field's first-derivative differencing. Supply only one and gridgen warns (`<f>: dx<f> provided without its pair ...`) and recomputes both by FD.
- Recognized fields: `dxU,dyU, dxV,dyV, dxW,dyW` (incompressible); plus `dxrho,dyrho, dxT,dyT` and the second derivatives `dxxU/V/W/T, ddyU/V/W/T, dxyU,dxyV` (compressible). The single second-derivative fields have no pair, so they are taken as-is when present.
- gridgen reports what it used: `compute_bf_derivatives: using N provided base-flow derivative field(s)`.

Supplying analytic derivatives avoids the finite-difference noise that interpolated CFD base flows can introduce near the wall (which otherwise inflates the wall shear by `O(N²)`). For a smooth analytic base flow the FD path is already accurate to ~0.1%, so this is mainly useful for noisy or analytically-known base flows.

### Step 2 — Configure and run

Per-case stability scripts live in `callers/benchmark/` (`blasius.m`, `blasius_compressible.m`, `sweptwing_flat.m`, `sweptwing_hump.m`, `m3j.m`, `resonator.m`, `synthetic_jet.m`, `blasius_newton.m`). They share an identical sectioned layout — only case-specific values differ:

| Section | What to set |
|---|---|
| **Spectral content** | `Stab.M`, `Stab.N`, `Stab.omega_0`, `Stab.beta_0` |
| **Inflow conditions** | `Stab.IC`, `Stab.phaseRef`, `Opt.linear`, `Stab.A0_fund` |
| **Boundary conditions** | `Opt.bc_top`; optional `Stab.bcw`/`Stab.bct` for actuation/forcing |
| **Buffer** | `Opt.buffer` (`'on'`/`'para'`/`'off'`), `Opt.xb` |
| **Solver options** | `Opt.plot`, `Opt.Conv`, `Opt.TH` |
| **Performance** | `Opt.lu_mode`, `Opt.parfor` |
| **Nonlinear solver choice** | `Opt.solver` (`'picard'`/`'newton'`) |
| **Run** | `Opt.rerun`, `Opt.save_output`, call to `main()` |
| **Post-processing** | `plot_amplitudes`, `plot_stability` |

Run a caller from the MATLAB command window:

```matlab
cd callers/benchmark/
blasius           % or sweptwing_flat / sweptwing_hump / m3j
```

---

## 4. Spectral content

```matlab
Stab.M = 5;            % temporal modes: indices -M ... +M (total 2M+1)
Stab.N = 0;            % spanwise modes: indices -N ... +N (total 2N+1)
Stab.omega_0 = 0.034;  % fundamental angular frequency (non-dimensional)
Stab.beta_0  = 0;      % fundamental spanwise wavenumber (0 = 2D)
```

Total number of modes: `nf = (2M+1)(2N+1)`. Each mode `(m,n)` has frequency `m*omega_0` and wavenumber `n*beta_0`.

### Single oblique mode

To track a single oblique wave without spanwise harmonics, set `N = 0` with `beta_0 > 0`:

```matlab
Stab.M = 1;
Stab.N = 0;
Stab.beta_0 = 0.1;  % oblique wavenumber
```

The solver assigns `betavec = beta_0 * sign(Mvec)` for correct conjugate symmetry.

To get the index of mode `(m,n)` in the solution arrays:

```matlab
j = ModeToModeNumber(m, n, Stab.M, Stab.N);
```

---

## 5. Inflow conditions

```matlab
Stab.IC       = 'ILST';   % local stability eigenvalue problem at inflow
Stab.phaseRef = 'pwall';  % 'pwall' = wall pressure reference | 'umax' = max u
```

`Stab.IC` selects how the inlet disturbance is built:

| value | inflow |
|---|---|
| `'ILST'` | local incompressible LST eigenmode |
| `'CLST'` | local compressible LST eigenmode (requires `Opt.compressible='on'`) |
| `'PSE'` | non-parallel linear-PSE seed instead of the local eigenmode |
| `'LOAD'` | warm-start from a checkpoint, or impose user-defined profiles from a `phiIC` file |
| `'ZERO'` | no inlet disturbance; the field is generated by wall/freestream forcing alone |

### Setting initial amplitudes

The simplest way is to set the fundamental amplitude as a scalar:

```matlab
Stab.A0_fund = 0.00125*sqrt(2);  % fundamental mode amplitude
```

`init.m` automatically builds the full amplitude vector with conjugate pairs. For per-mode control, set `Stab.A0` directly:

```matlab
Stab.A0 = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Stab.A0(ModeToModeNumber( 1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
Stab.A0(ModeToModeNumber(-1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
```

Always set both the mode and its conjugate (index `-m`, `-n`) to ensure real physical fields.

### Amplitude convention

`Stab.A0_fund`, `Stab.A0`, and the wall/freestream forcing `Stab.bcw`/`Stab.bct` are the **physical peak amplitude** `2|û|` of the real disturbance `q'(x,y,t) = q̂(x,y) e^{-iωt} + c.c.`, normalized by the freestream velocity. Internally the solver seeds the complex amplitude `û = A0/2`, and post-process reports `StabRes.A = 2|û|` — so **what you put in is what you get out**:

- `StabRes.A` at the inlet equals `A0_fund`;
- the reconstructed wall velocity `StabRes.v .* StabRes.A` over an actuator equals the `bcw` you imposed.

So set `A0_fund`/`bcw` to the physical amplitude you want and compare directly.

**Porting a case from DeHNSSo-v1 (Westerbeek et al. 2024).** The *stored* meaning of `Stab.A0` changed, so do not paste a v1 literal verbatim. In v1 the code seeded `û = A0` directly (no halving) and doubled only at output, so `Stab.A0` held the complex amplitude `|û|`; the user was expected to type **half** the physical amplitude in the caller, which is why every v1 caller reads e.g. `Stab.A0(...) = 0.00125*sqrt(2); % Half the desired inflow amplitude` or `= 3.5e-3/2`. v2 keeps the seeded `û` identical but moves that `/2` into the solver, so `Stab.A0` now holds the physical `2|û|`. Carry over the **desired physical** amplitude (the value *before* the v1 `/2`), not the raw v1 number:

| case | v1 caller (`A0 = \|û\|`) | v2 caller (`A0 = 2\|û\|`) |
|------|------------------------|--------------------------|
| Blasius TS | `0.00125*sqrt(2)` | `0.0025*sqrt(2)` |
| swept-wing CFI | `3.5e-3/2` | `3.5e-3` |

Pasting the v1 literal into v2 unchanged halves it a second time (a 4× weaker fundamental for nonlinear runs). The physics and results are otherwise identical between v1 and v2 for the same physical amplitude — only the input label moved.

### Restarting from a previous solution (`Stab.IC = 'LOAD'`)

`Stab.IC='LOAD'` warm-starts the whole domain from a converged solution instead of the ILST inflow. Its main use is **amplitude continuation**: converge at a low `A0_fund`, then restart at a higher one. The loaded field is rescaled by the amplitude ratio, giving Picard a far better initial guess than a cold start.

**The file must be a checkpoint** (`Opt.save_checkpoint`), *not* a results file (`Opt.save_output`) — see §10. Post-processing strips `q`/`phi`/`R`, drops the conjugate modes and normalizes the shapes, so a results file cannot restore the solver state. Loading one raises `inflow:ICfileNotCheckpoint`.

Run 1 — converge at the seed amplitude and write the checkpoint:

```matlab
Stab.IC             = 'ILST';
Stab.A0_fund        = 0.00125*sqrt(2);
Opt.anderson        = 5;                       % see §8
Opt.save_checkpoint = '<...>/ckpt_seed.mat';   % NOT save_output
```

Run 2 — continue at the higher amplitude:

```matlab
Stab.IC      = 'LOAD';
Stab.ICfile  = '<...>/ckpt_seed.mat';
Stab.A0_fund = 0.002;                          % target; field is rescaled by 0.002/A0_seed
Opt.anderson = 5;
```

`inflow` reports what it restored:

```
LOAD: full-domain warm-start (11 modes -> 11 incl. conjugates, 40 x 800), A0 scale = 1.131
```

If it instead prints `LOAD: inflow-only (station 1) — no downstream warm-start`, the checkpoint's `nx` does not match the current grid and only the inflow station was used — there is no continuation benefit.

**Anderson is not optional here.** Near saturation plain Picard can diverge at the target amplitude even though the continuation gives it a good initial guess; set `Opt.anderson = 5` on both runs (§8).

A complete working example is `callers/benchmark/blasius.m` with a `restart` flag. On the Blasius benchmark (800×40, M=5), continuing `A0 = 0.00125√2 → 0.002` reaches the same solution as a cold run at 0.002 (`max(A)` agreeing to 7e-05) in 43 iterations instead of 57.

The same mechanism serves parameter sweeps: restart at a new `omega_0`, `beta_0` or grid-independent parameter from a nearby converged state. The checkpoint must come from a run with the same `M`/`N` and grid.

### User-defined inflow profiles (`Stab.IC = 'LOAD'` with `phiIC`)

The same `'LOAD'` switch also accepts a file holding **only** a variable named `phiIC`, in which case it does not warm-start the domain: it imposes your own profiles at the inflow station, exactly where `ILST`/`CLST` would put the computed eigenmode. Use it when the inflow shape comes from somewhere else, for example a DNS slice, an experimental profile, an external stability solver, or a hand-built superposition.

Give `phiIC` one column per variable, with the wall-normal index down the column, following the grid convention index 1 = freestream, index `ny` = wall. For a single mode that is an `ny x nvar` array (`nvar = 4` incompressible, `5` compressible), and `phiIC_modes` names the mode it belongs to:

```matlab
phiIC       = [u, v, w, p];      % ny x 4   (compressible: [u, v, w, rho, T], ny x 5)
phiIC_modes = [1 0];             % this profile is mode (m,n) = (1,0)
save('my_inflow.mat','phiIC','phiIC_modes')

Stab.IC     = 'LOAD';
Stab.ICfile = 'my_inflow.mat';
```

**Several modes** stack along the third dimension, `ny x nvar x k`, one `(m,n)` row of `phiIC_modes` per slice. The layout of each slice never changes: it is always the same `nvar` columns. Only `phiIC_modes` says which mode a slice is, so any modes may be given, in any order:

```matlab
phiIC        = zeros(ny, 4, 2);   % ny x nvar x 2
phiIC(:,:,1) = [u, v, w, p];      % this slice is (2,0)
phiIC(:,:,2) = [u, v, w, p];      % this slice is (0,0)
phiIC_modes  = [2 0;              % row 1 -> (m,n) = (2,0)
                0 0];             % row 2 -> (m,n) = (0,0)
save('my_inflow.mat','phiIC','phiIC_modes')
```

The solver echoes what it mapped, which is worth checking in the log:

```
  LOAD: 2 user profile(s) mapped to mode(s) (2,0) (0,0)
```

**Set `Stab.A0` for every mode you supply.** Only modes with `Stab.A0 > 0` are read, so a profile for a mode left at zero amplitude loads and is then silently ignored. `Stab.A0_fund` seeds the fundamental *alone*, so supplying more than one mode means setting the `Stab.A0` vector yourself instead (`init` builds `A0` from `A0_fund` only when you have not supplied `A0`). Use `ModeToModeNumber(m, n, Stab.M, Stab.N)` for the indices rather than counting by hand:

```matlab
nf  = (2*Stab.M+1)*(2*Stab.N+1);
j20 = ModeToModeNumber( 2, 0, Stab.M, Stab.N);
jm20= ModeToModeNumber(-2, 0, Stab.M, Stab.N);
j00 = ModeToModeNumber( 0, 0, Stab.M, Stab.N);
Stab.A0 = zeros(1, nf);
Stab.A0([j20 jm20]) = 5e-5;   % (2,0) AND its conjugate
Stab.A0(j00)        = 2e-5;   % (0,0) is its own conjugate: one entry, one slice
```

An oscillating mode needs its conjugate `(-m,-n)` given a nonzero `Stab.A0` as well; the profile itself is rebuilt by conjugation, so no extra slice is required. The mean-flow distortion `(0,0)` is its own conjugate, so it takes a single entry. Since the stored amplitude is discarded, `Stab.A0` is what sets the relative weighting between the modes you supply.

**Only the shape is taken from the file.** Each mode is renormalized by `max|u|` and then rescaled by `Stab.A0`, so the amplitude you saved is discarded and the run's amplitude is still set by `Stab.A0_fund` as usual. The phase is likewise re-referenced according to `Stab.phaseRef`. This means a profile does not need to arrive normalized or phase-aligned.

Two further conveniences: the conjugate half of the mode list is rebuilt automatically by conjugation, so those rows are overwritten whatever you put in them; and in the compressible case `p'` is derived from the equation of state, `p' = (rho'*Tbar + rhobar*T')/(gamma*Ma^2)`, so supply `rho` and `T` and not `p`.

`phiIC_modes` is a `k x 2` list of `(m,n)` pairs, one per profile, so you supply **only the modes you have**. Modes not listed are left undisturbed, and the solver reports what it mapped:

```
  LOAD: 1 user profile(s) mapped to mode(s) (1,0)
```

Two older layouts are still accepted. A **flat row stacking**, `k x nvar*ny`, concatenates the variables along each row instead of using columns, `[u(:).', v(:).', w(:).', p(:).']`; it behaves identically. And **omitting `phiIC_modes` entirely** requires `phiIC` to span the whole mode list, exactly `Mode.nf = (2M+1)(2N+1)` rows with conjugates included, row `j` being mode `j`. For `M=1, N=0` those three rows are `(m,n) = (-1,0), (0,0), (1,0)`. That is the form a solver template comes in, but the template must come from a **checkpoint** (`Opt.save_checkpoint`), not a results file: `postprocess` deletes `phi` and drops the conjugate modes, so a `save_output` file has neither the field nor the right number of rows.

```matlab
% run once with Opt.save_checkpoint set, then:
ref   = load('<...>/ckpt_seed.mat');
phiIC = ref.StabRes.phi(:,:,1);        % nf x nvar*ny, correct layout by construction
save('my_inflow.mat','phiIC')
```

Whichever layout, the wall-normal length must match the run's `ny`.

Only the modes with `Stab.A0 > 0` are read. A **linear** run solves one mode: for `M=1, N=0` the log reads `M = 1, N = 0 (1 mode)` and a single `Solving 1 modes`, because `(0,0)` is never seeded (there is no nonlinear forcing to generate a mean-flow distortion) and the conjugate `(-1,0)` is filled by conjugation rather than solved. So one profile, for `(1,0)`, is all a linear run needs.

**Every block of that row is required, including the pressure.** The inflow is imposed as a Dirichlet condition on the whole state vector at station 1, so a missing block is imposed as zero rather than ignored. Filling only `u, v, w` and leaving `p` zero shifts the downstream amplitude by 15%, and supplying `u` alone by 39%, whereas the complete row reproduces the reference run to 1e-12. If the external data has no pressure, reconstruct it rather than zero-fill; note also that `phaseRef='pwall'` references the wall pressure, which degenerates if `p` is zero.

Building the array from scratch is otherwise forgiving, and three things help: the conjugate rows may be left as zeros, the overall amplitude is irrelevant, and so is an overall phase. Reloading a template with the conjugate row zeroed, scaled by 1000, or rotated by an arbitrary phase reproduces the original run to round-off in each case.

A file containing a `StabRes` variable takes the warm-start path above instead; `phiIC` is checked first. If the file has neither, `inflow` raises `ICfile must contain a variable named 'StabRes' or 'phiIC'`.

### No inflow disturbance (`Stab.IC = 'ZERO'`)

`Stab.IC='ZERO'` seeds **no** disturbance at the inflow — the entire field is generated by wall/freestream forcing (`Stab.bcw` / `Stab.bct`). Use it when the quantity of interest is the receptivity of a wall source *in isolation*: a synthetic-jet or plasma actuator, a roughness element, or the resonator impedance BC, with **no incoming instability wave**. It works in both incompressible and compressible mode, and needs no `Stab.A0_fund` (there is no inflow mode to scale).

```matlab
Stab.IC  = 'ZERO';
Stab.bcx = x;
Stab.bcw = zeros(3, numel(x), (2*Stab.N+1)*(2*Stab.M+1));
Stab.bcw(2, :, j_fund) = vjet;   % e.g. a wall-normal velocity actuator on the fundamental
```

The inflow plane is clamped to zero, so the active modes come from the wall forcing alone. `inflow` warns (`inflow:zeroNoForcing`) if there is no forcing at all, since the field would then be identically zero. A minimal example is `callers/benchmark/synthetic_jet.m` — a wall synthetic jet driving a Blasius BL with no incoming TS wave.

---

## 6. Boundary conditions

By default all BCs are homogeneous (zero perturbation at wall and freestream). You can prescribe inhomogeneous BCs to model actuation or external disturbances.

### Arrays

```matlab
Stab.bcx          % (1 x nx_bc)          x-locations of BC profile
Stab.bcw(k,ix,j)  % (n_bc x nx_bc x nf)  wall BC
Stab.bct(k,ix,j)  % (n_bc x nx_bc x nf)  freestream BC
```

`n_bc` is 3 incompressible and 4 compressible, and the channel `k` maps differently at the two boundaries (see `bc_rows`):

| `k` | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| wall, incompressible | u | v | w | — |
| wall, compressible | u | v | w | T |
| freestream, incompressible | u | v | w | — |
| freestream, compressible | u | w | rho | T |

`add_wall_forcing` sizes `Stab.bcw` from `Opt.compressible` for you, so this matters only when writing the arrays by hand.

### Wall actuation (blowing/suction)

Wall forcing follows a **grid-first** workflow: you describe what is at the wall **in the gridgen script**, gridgen refines the mesh to resolve it and stores the profile on the grid, and the caller applies it in one line. Location, width and shape are defined **once**, in gridgen; only the strength (amplitude, phase) lives in the caller. The reference pair is `gridgen/benchmark/synthetic_jet.m` + `callers/benchmark/synthetic_jet.m`.

This works on **any grid** — cartesian, curvilinear, and unstructured. On a flat plate `xc`/`x_data` are the physical `x`; on a curved or unstructured wall (swept-wing, hump, M3J/CSV) they are still physical `x` — gridgen samples the profile at the physical wall-x per station, so it lands geometrically correctly even though the grid runs in arc-length.

**Step 1 — gridgen: define the forcing.** `params.forcing` describes what sits at the wall. Two choices decide the form: where the **shape** comes from (the built-in cosine, or an array you supply), and **how many** sources there are (one, or an array of them). That gives seven forms:

| # | shape | how many | spec |
|---|---|---|---|
| 1 | cosine | one | `struct('xc', xc, 'w', w)` |
| 2 | cosine | N, same width | `struct('xc', [xc1 xc2 xc3], 'w', w)` |
| 3 | cosine | N, own widths | `struct('xc', [xc1 xc2], 'w', [w1 w2])` |
| 4 | tabulated | one | `struct('xc', xc, 'w', w, 'v', v_data)` |
| 5 | tabulated | one, uneven samples | `struct('xc', xc, 'w', w, 'x', x_data, 'v', v_data)` |
| 6 | tabulated | N | `struct('xc', [xc1 xc2 xc3], 'w', w, 'v', v_data)` |
| 7 | tabulated | one, **placed by the data** | `struct('x', x_data, 'v', v_data)` |

`'shape','cosine'` may be given explicitly in 1-3, but it is the default and `'cosine'` is the only value. In 4-6 `v_data` may be complex, carrying phase.

**Forms 1-6 place the source themselves.** `xc` is the centre and `w` the HALF-width, so the footprint spans `2w`; a supplied array is stretched onto `[xc-w, xc+w]` and normalised to unit peak. Its own units, origin and sample order therefore never enter, and in form 5 `x_data` is read for the *relative spacing* of the samples alone. Pass a profile tabulated in millimetres about its own centre exactly as it comes out of the file.

**Form 7 is the exception**: with no `xc`/`w`, the centre and half-width are taken from where `|v|` first and last exceeds `1e-6` of its peak, so there `x_data` **is** the position and must already be in the grid's coordinate. Use it when the distribution is the specification and has no natural centre and width to quote. Two consequences: a profile that is exactly zero at its end points reports a slightly narrow `w`, and a units mistake moves the source rather than being ignored.

Whichever form, a source lying entirely off the domain raises `main_gridgen:forcingOffDomain`, and one hanging over an edge warns with `forcingClipped`.

```matlab
% one actuator
params.forcing = struct('xc', 1379.1, 'w', 18);

% a measured neck profile, as it comes out of the file
D = load('hAbs_xProfile_TM295hz.mat');           % xh in metres, h in [0,1]
params.forcing = struct('xc', Xc, 'w', s_nd/2, 'v', D.h);

% an array of three identical necks
params.forcing = struct('xc', Xc + (0:2)*(s_nd + d_nd), 'w', s_nd/2, 'v', D.h);

params.spread = 30;     % how far the refinement reaches
params.refine = 5;      % dxi at the centre vs equidistant
```

`refine` is the friendlier form of `ag` (see the clustering table above). `spread`, how far the refinement reaches, defaults to the forcing half-width, which is right for an actuator; set it explicitly only when the refinement must reach further than the feature itself, as the resonator neck does.

**Sources that differ from each other** go in a cell array, one struct each. This composes with the table: a vector `xc` inside a cell entry still expands, so the example below is four sources.

```matlab
params.forcing = { struct('xc', [xc1 xc2 xc3], 'w', 18, 'v', v_data), ...   % three alike
                   struct('xc', xc4, 'w', 24) };                            % one different
```

gridgen stores `StabGrid.forcing` (the unit-peak profile(s) sampled on the refined grid) plus `StabGrid.cluster.xc/.w`. The distribution's `v` values need only encode the **shape**; the amplitude is set in the caller. Note the fixed station budget is shared across all centers, so raise `params.n_xi_new` when refining several spots.

**Step 2 — caller: apply the strength.** `add_wall_forcing` reads `StabGrid.forcing`, scales it, and does all the `Stab.bcw` setup — component, mode index, `Stab.bcx`, sizing, and the conjugate mode (added automatically for nonlinear runs):

```matlab
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 1e-4));   % xc/w/shape from gridgen
```

`spec` fields: `A` (amplitude, physical `2|û|` like `Stab.A0`, required), optional `phi` (phase), `comp` (1 u / 2 v / 3 w, and 4 T compressible; default 2), `mode` (`[m n]`, default `[1 0]`), and `which` to pick one forcing when gridgen stored several. The forcing is **added** to any existing `Stab.bcw`, so call once per source to stack several:

```matlab
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 2e-4, 'which', 1));
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 5e-4, 'which', 2, 'phi', pi/2));
```

`which` is required once the grid holds more than one forcing (the call errors otherwise) and unnecessary for a single one. You can also override the profile here for a quick test — pass `x`/`v` (distribution) or `xc`/`w` (actuator) — but the mesh is only refined where gridgen put it, so an override off the refined patch is coarser (`add_wall_forcing` warns if it lands on too few stations).

Because the mesh is static, the width is capped by what you refined for: a narrower forcing than `spread` stays resolved; a **wider** one, or a moved center, needs a new grid. Size the gridgen `w`/support for the widest case you intend to run, then sweep amplitude/phase freely with no re-mesh.

**Low-level fallback.** `Stab.bcw`/`Stab.bcx` is ultimately just an array over the wall, cubic-interpolated onto `Grid.xun` (zero outside `bcx`). You can set it by hand if you need full control — but you are then responsible for the mode/conjugate/convention handling and for refining the grid where the profile varies:

```matlab
j_fund = ModeToModeNumber(1, 0, Stab.M, Stab.N);
j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
Stab.bcx             = x_data;                        % your x locations
Stab.bcw             = zeros(3, numel(x_data), (2*Stab.N+1)*(2*Stab.M+1));
Stab.bcw(2,:,j_fund) = v_data;                        % k=2 = wall-normal; physical 2|û|
Stab.bcw(2,:,j_conj) = conj(v_data);                 % nonlinear runs only
```

`k=1/2/3` = `u/v/w`.

### Helmholtz-resonator wall (impedance BC)

A Helmholtz-resonator neck in the wall is imposed as a reactive impedance BC through `Opt.hr`. Over the mouth the wall `v'`-row of the operator is replaced by

```
v'(x) − Y·p'(x) = 0,
```

coupling the mouth wall-normal velocity to the wall pressure through the resonator admittance `Y` (the closed-loop resonator response). It is applied to the base operator only, so the `(0,0)` mean-flow distortion keeps no-slip while the fluctuations see the resonator.

Give the admittance either as a magnitude and phase, which keeps the phase visible, or directly as `Y`:

```matlab
Opt.hr = struct('mag', 12.7688, 'phi', pi);   % |Y| and arg Y [rad]
```

`setup_boundary` resolves this once per run: it forms `Y = mag·exp(i·phi)`, takes the mouth from `StabGrid.forcing` when no `mask` is given, and reports the result. The full set of fields:

```matlab
Opt.hr = struct( ...
    'Y',     Ynd, ...          % admittance: a scalar, or one value per mouth node for a profile Y(x)
    'mask',  in_neck, ...      % optional, logical nx: the mouth stations (else from StabGrid.forcing)
    'x',     StabGrid.xun, ... % optional: grid-exact flux-preserving edge taper
    'taper', 0.1, ...          % optional: raised-cosine mouth edge taper (0 = hard)
    'uslip', false);           % optional: free-slip u' over the mouth
```

Cluster the streamwise mesh at the neck in gridgen as for a wall actuator, but note that `params.forcing.w` is the **physical mouth half-width**, not a clustering knob: use `params.spread` to set how far the refinement reaches. Keeping the two apart means `StabGrid.forcing` carries the true mouth, so the same grid also serves the actuator path. The reference pair is `gridgen/benchmark/resonator.m` + `callers/benchmark/resonator.m` (validated against Michelis 2023).

### Freestream forcing

Prescribe a perturbation at `y=H` (e.g., an acoustic wave or free-stream turbulence):

```matlab
Stab.bct(1,:,j_fund) = A_fst * ones(1, numel(x_dom));  % u at y=H
Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));
```

**Note:** for the `(0,0)` mean-flow distortion mode, the `v` freestream row carries the transpiration boundary condition and is not overridden by `bct`.

### Freestream BC type

`Opt.bc_top` is a cell array with one entry per freestream channel:

```matlab
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % incompressible default: all zero at y=H

% BC types:
%   'H_DR'  = homogeneous Dirichlet  (q = 0)
%   'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%   'H_NM'  = homogeneous Neumann    (dq/dy = 0)
%   'free'  = no imposed BC; the governing equation is kept on that row
```

**The channels are not the same variables in the two formulations**, because the two boundaries constrain different states. `Opt.bc_top` must have exactly this many entries or `init` raises `init:bcTopSize`:

| | channel 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| incompressible freestream (`Opt.bc_top`) | `u` | `v` | `w` | — |
| compressible freestream (`Opt.bc_top`) | `u` | `w` | `rho` | `T` |
| incompressible wall (`Stab.bcw`) | `u` | `v` | `w` | — |
| compressible wall (`Stab.bcw`) | `u` | `v` | `w` | `T` |

`Stab.bct` (inhomogeneous freestream values) follows the freestream row, so compressibly `Stab.bct(2,...)` forces `w`, not `v`. Note the compressible freestream leaves `v` free: it carries the displacement, matching the heritage far-field treatment. Incompressibly the pressure has no freestream BC, and the MFD `v` is always free there (displacement velocity). The single source of truth for these maps is `src/boundary/bc_rows.m`, which `LHS`, `LHS_compressible` and `RHS` all read.

**Receptivity example** — acoustic wave forcing v at freestream:
```matlab
Opt.bc_top = {'H_DR', 'IH_DR', 'H_DR'};  % u=0, v=forced, w=0
Stab.bct(2,:,j_fund) = 1e-4 * ones(1, numel(x_dom));
Stab.bct(2,:,j_conj) = conj(Stab.bct(2,:,j_fund));
```

For backward compatibility, scalar strings `'dirichlet'` and `'neumann'` are auto-converted to `{'H_DR','H_DR','H_DR'}` and `{'H_NM','H_NM','H_NM'}` respectively.

The ILST inflow computation uses the same BC type.

---

## 7. Linear vs nonlinear mode

```matlab
Opt.linear = 'on';   % linear: NLT off, single iteration, AF = 1
Opt.linear = 'off';  % nonlinear (default): iterative solve with NLT
```

In linear mode the solver converges in one iteration (equivalent to linear PSE). Use this for parametric sweeps, N-factor computations, or any analysis where nonlinear interactions are not needed.

---

## 8. Solver options

### Buffer zone

```matlab
Opt.buffer = 'on';    % amplitude damping (prevents end reflections)
Opt.buffer = 'para';  % parabolization (zeroes streamwise diffusion, LPSE-like)
Opt.buffer = 'off';   % no buffer
Opt.xb     = 85;      % buffer start as % of domain length
```

### Convergence

```matlab
Opt.Conv  = 1e-4;   % convergence criterion (relative amplitude change)
Opt.ConvF = 100;    % relaxation: use Opt.Conv*Opt.ConvF until AF=1
Opt.AFg   = 1.1;    % amplitude factor growth rate
Opt.AMAX  = 0.1;    % maximum amplitude factor
Opt.TH    = 1e-11;  % nonlinear term activation threshold
```

### Picard acceleration and under-relaxation

The nonlinear convergence loop is Picard iteration with AF (amplitude factor) continuation. Two optional knobs help it converge: Anderson acceleration (speed) and under-relaxation (stability). There is no Newton fallback: `Opt.solver` selects Picard or Newton up front and they never mix (see below).

**Anderson acceleration** (Walker & Ni, 2011). A least-squares fixed-point accelerator that combines the last `m` Picard iterates to produce a faster update.

```matlab
Opt.anderson = 0;   % off (default — plain Picard)
Opt.anderson = m;   % history depth m (typical: 3-6)
```

How it behaves:

- Maintains a rolling buffer of the last `m` `(q_old, residual)` pairs.
- At each iteration, solves a small `m × m` least-squares problem and forms an accelerated step from the difference history. The solve is damped, so when the residuals become nearly parallel near saturation it falls back to a plain Picard step.
- Safety fallback: if the accelerated step diverges (magnitude > 1.5× plain Picard), the iteration keeps the plain Picard step instead.
- The history is reset automatically when `AF` changes or the set of active modes changes (the fixed point itself has moved).

Often turns the linear convergence of plain Picard into a superlinear one near saturation, halving the iteration count on M=5 nonlinear cases.

**Under-relaxation** (`Opt.relax`). A fixed damping of the Picard step, `q ← q_old + ω·(q_new − q_old)` with `0 < ω < 1`.

```matlab
Opt.relax    = 0.4;   % omega; omit (or >= 1) to disable — default: off
Opt.anderson = 0;     % turn Anderson OFF: the two interfere
```

Use it when plain Picard **oscillates** near saturation — the maximum amplitude and the residual bounce and never settle (the unstable Picard eigenvalue is oscillatory, which Anderson cannot damp). Under-relaxation shrinks the effective step so a marginally-repelling fixed point becomes attracting, giving a monotone descent. It is a *stabilizer*, not an accelerator: convergence is slower (more iterations) but reaches tolerance where plain Picard would not. Use with `Opt.anderson = 0` (they fight). A typical value is `ω ≈ 0.4`, needed for high-amplitude nonlinear cases such as the compressible oblique breakdown that saturates near `max|u'| ≈ 0.16`.

**No Newton takeover.** `Opt.solver` is either `'picard'` or `'newton'`, never a mix — Picard no longer switches to Newton when it stalls. If Picard stalls at AF < 1 (typically near saturation) it reports once and keeps going; for high amplitudes use `Opt.solver='newton'` (pseudo-arclength continuation, below) instead. `Opt.newton_takeover` is gone.

**There is no fixed-amplitude Newton.** The old backends (JFNK, `direct`, `matlab_direct`, `lu_krylov`) and the fixed-amplitude Newton driver have been removed: a cold-start solve at a fixed high amplitude does not converge, which is precisely what the continuation exists to solve. Setting those `Opt.solver` values now raises an error. `Opt.solver='newton'` means `src/newton/solve_newton.m`.

### Pseudo-arclength continuation (`solve_newton`)

For high-amplitude K-type cases (Scholten), a direct solve at the target amplitude does not converge from a cold start. `src/newton/solve_newton.m` continues the solution along the branch in the inflow amplitude `A0_fund` instead of trying to reach the target in one go. Call it directly (it does its own `init`/`inflow`/`setup_boundary`/`newton_context`) rather than through `main`:

```matlab
[StabRes, info] = solve_newton(Stab, StabGrid, Opt);
```

The default scheme is a **modified pseudo-arclength continuation**. The predictor is a Keller step (tangent plus secant curvature); the corrector is a modified Newton on the full coupled real-form Jacobian, factorized with MUMPS. The factor is **held across steps** (the corrector then needs only back-solves) and refactored when it ages or the corrector stalls, since one factorization costs about 230× a back-solve and the state moves little between steps.

Both departures from the textbook scheme of Dijkstra et al. (2014, Sec. 2.2) can be switched off independently:

```matlab
Opt.newton_corrector = 'modified';   % default; 'full' refactorizes every corrector iteration
Opt.newton_predictor = 'curvature';  % default; 'tangent' is the first-order Keller step
```

Setting both alternatives together (`'full'` + `'tangent'`) reproduces classical Dijkstra continuation. It converges quadratically and needs only one or two corrector iterations per step, but pays a factorization for each: measured on a 1050×40 M=2 Blasius branch it took 9 factorizations and 52.9 min against 2 and 23.9 min for the default, for solutions agreeing to five significant figures. Use it to check the default, not for production. Under `'full'`, `newton_reuse` and `newton_it_thresh` are meaningless and warn once if set.

```matlab
Opt.newton_target    = 2.828e-3;  % target A0_fund (required) — here u_rms = 0.004
Opt.newton_start      = 1.20e-3;   % continuation start (default Stab.A0_fund)
Opt.newton_seed2      = 1.25e-3;   % second seed point (default 1.04*lam0)
Opt.newton_reuse      = 8;         % reuse steps before a forced refactor; 1 = refactor every step
Opt.newton_it_thresh = 8;         % refactor next step if the corrector exceeded this
Opt.newton_maxit     = 25;        % max back-solves before the factor is declared stale
Opt.newton_dsmax     = 3;         % max arclength step
Opt.newton_tol       = 1e-6;      % corrector residual tolerance
Opt.newton_seedfile  = '...';     % cache for the two seed points ([] = recompute)
Opt.newton_ckfile    = '...';     % checkpoint written each step ([] = none)
Opt.newton_resume    = 'on';      % resume from newton_ckfile
Opt.newton_mumps_relax = [100 200 40];  % MUMPS ICNTL(14) retry ladder
Opt.newton_arclength = 'absolute'; % arclength metric; see below
```

**Arclength metric.** `'absolute'` (default) measures `ds^2 = |dq|^2 + dlam^2`, which is badly conditioned when the state norm dwarfs the parameter: on the M3J x505 branch `|q|~8.2` against `lam~6.5e-4`, so `ds` tracks the state and the amplitude crawls at `dA0/ds = 2.2e-5`, needing about 275 steps to traverse. `'relative'` measures `ds^2 = |dq|^2/|qc|^2 + (dlam/lc)^2`, so both contribute comparably. Same branch, different parameterisation — it changes the step count severalfold, not the answer.

`newton_reuse` is the reuse policy and the main performance knob: `newton_reuse=1` refactors every step (most robust), larger values trade factorizations for more corrector back-solves. Which wins is case-dependent — measure it, and only trust a timing from a log whose first step is `STEP 1` (a resumed run's reported total covers only its own segment).

`lambda` enters the residual only through the inflow Dirichlet BC, hence affinely: `F(q,λ) = newton_residual(q,ctx0) + (λ−λ0)·F_λ`. `F_λ` is computed once by finite difference and its linearity is checked and reported at startup.

**Memory.** The coupled Jacobian is factorized whole (2N ≈ 1.7M for M=4 on a 1050×40 grid). Peak RSS during factorization reaches ~150 GB at M=4 and falls back to ~35 GB once the factor is formed — each refactor re-spikes. M=4 needs a workstation; a 16 GB laptop is limited to roughly M ≤ 3.

### LU mode

```matlab
Opt.lu_mode = 'full';         % cache every LU factor — fastest reuse, unbounded memory
Opt.lu_mode = 'auto';         % LRU-cached LUs, bounded memory (Opt.lu_max_cache, default max(8, n_modes)) (default)
Opt.lu_mode = 'seq';          % one LU at a time: factor, solve its modes, discard — minimal memory
Opt.lu_mode = 'none';         % plain backslash each time, no caching (slowest, baseline)
Opt.lu_mode = 'lapack_band';  % banded LU MEX (zgbtrf/zgbtrs) with per-mode factor caching —
                              %   fastest single solve (~3x); one-time MEX build (see below)
Opt.lu_mode = 'block_thomas'; % dual-Schur block-Thomas (Paredes 2024) — COMPRESSIBLE operator only
% Opt.lu_max_cache = 3;       % cache size for 'auto' (and the banded / block-Thomas factor caches)
```

### Fill-reducing ordering

Applies to the sparse-LU back-ends (`full`, `auto`, `seq`); `lapack_band` and
`block_thomas` have their own storage and ignore it.

```matlab
Opt.lu_ordering = 'auto';     % (default) symrcm for the compressible operator,
                              %   COLAMD for the incompressible one
Opt.lu_ordering = 'symrcm';   % force reverse Cuthill-McKee
Opt.lu_ordering = 'colamd';   % force MATLAB's 4-output lu (COLAMD + column pivoting)
```

The operator is a chain of station blocks, so it is banded and bandwidth reduction is
what limits the fill. Measured fill, `nnz(L+U)/nnz(M)`, at 40/60/80 stations:

| operator | COLAMD | symrcm |
|---|---|---|
| incompressible | **6.7 / 6.8 / 6.8** | 7.1 / 7.3 / 7.3 |
| compressible | 25 / 43 / 50 | **36.4 / 37.1 / 37.7** |

`symrcm` stays flat in `Nx` while COLAMD grows, crossing near `Nx = 55`; at `Nx = 200`
COLAMD gave nearly double the fill. Hence the default. Forcing `colamd` is worth it on
small compressible cases, where it is both leaner and several times faster. AMD and
nested dissection were measured at up to 79x and 174x on these operators and are not
offered.



`full`, `auto`, `seq`, and `none` all use MATLAB's sparse LU and differ only in how many factors they keep in memory: `full` keeps all of them, `auto` keeps the last `lu_max_cache`, `seq` keeps one at a time, `none` keeps none. `lapack_band` and `block_thomas` replace that sparse LU with a banded or block factorization of the same operator, then cache the factors the same way.

**`block_thomas`** is for the **compressible** operator only. The incompressible continuity rows have a zero diagonal, which makes the unpivoted block factorization silently wrong (verified ~2× error), so the solver raises an error if you select it on an incompressible run. Sub-options:

```matlab
Opt.bt_np    = 2;         % number of partitions (default 2)
Opt.bt_ooc   = 'on';      % out-of-core: stream the factor to disk, minimal RAM (default off)
Opt.bt_apply = 'inverse'; % how the diagonal station block is applied (default 'lu')
Opt.lu_cache = 'off';     % one-shot factor+solve, no caching (default: cache like the other modes)
```

`bt_apply` chooses what is stored for each diagonal station block `G` of the recursion:

| Value       | Stored           | Applied by          |
| ----------- | ---------------- | ------------------- |
| `'lu'`      | pivoted LU of `G` (default) | triangular solves |
| `'inverse'` | `G^-1`                      | matrix product    |

`'lu'` is the default. It keeps the pivoting that forming the inverse discards, which matters
because a few station blocks are close to singular (RCOND down to 5e-18 on the Sec. 5.2
verification case). Measured there, over 80 iterations on 800x40 with M=5: the solution differs
from the cached sparse-LU path by 8.1e-9 with `'lu'` and 4.8e-6 with `'inverse'`, and the residual
of the per-mode solves falls from 8.9e-4 to 2.8e-6, for about 4% more run time. `'inverse'`
follows Paredes et al. (2024). Memory is the same either way, since the LU factors pack into one
block, and both forms stream out of core, the LU record carrying one extra vector of pivots
per station.

**Picking quickly:**

| Run                                  | Recommended                                |
| ------------------------------------ | ------------------------------------------ |
| Linear single-mode (TS)              | `lapack_band` (or `full`)                  |
| Nonlinear multi-mode, plenty of RAM  | `full` (caches every reusable factor)      |
| Memory-constrained                   | `seq`, or `auto` with small `lu_max_cache` |
| Compressible, large grid, low RAM    | `block_thomas` with `bt_ooc='on'`          |
| Validation / debugging               | `none` (reference, no caching)             |

### Building the `lapack_band` MEX (one-time setup per machine)

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

### Parallel acceleration

Requires the MATLAB Parallel Computing Toolbox. If omitted, auto-detected at startup.

```matlab
Opt.parfor = 'off';  % 'on'  parfor over LHS stations, mode solves, and the NLT forcing
```

**`parfor`:** parallelizes three loops — the LHS station assembly (Phase 1), the per-mode solve/factorization (Phase 3), and the **nonlinear forcing** (`NLT_HNS`, the mode-to-mode coupling). Each output mode's forcing is independent, so the forcing loop scales with the number of physical modes. That forcing is the dominant cost of many-mode nonlinear runs, because the cubic three-wave term grows as the mode count squared. So `Opt.parfor='on'` is the main speedup for large `M`/`N` cases, for example an M=N=4 oblique breakdown. Use a Threads pool: its workers share the broadcast fields and factor objects instead of copying them. This matters most for the direct solve, where each process worker would otherwise hold its own copy of the LU/factor objects and the memory cost would grow with the worker count.

### Output control

```matlab
Opt.plot = 'on';   % 'on'  show amplitude and convergence plots during iteration (default)
Opt.plot = 'off';  % 'off' suppress all intermediate plots (recommended for batch runs)
```

### Iteration history

`Opt.iter_history = 'on'` records, for every iteration, the modal amplitudes, the
convergence measure `max_dal`, the amplitude factor, and the mode-summed
perturbation kinetic energy of Eq. (24) of DeHNSSo-v1. The record is returned in
`StabRes.iter_history`.

```matlab
Opt.iter_history = 'on';   % record the per-iteration history ('off' by default)
Opt.xloc         = 1700;   % station for the local energy measure
```

Two energies are stored. `E_tot` integrates over the buffer-free domain and is
always available. `E_loc` is the same quantity at a single station, and it is
computed only if you name one through `Opt.xloc`; left unset, `E_loc` is `NaN`.
There is no default station, because a fixed one cannot suit every domain: the
swept-wing benchmark spans x = 219 to 1360, so a station at 1700 would lie
outside it. A station you do name must lie inside the buffer-free range, and the
run stops if it does not, rather than reporting the energy somewhere else. The
value is interpolated to `Opt.xloc` rather than snapped to the nearest node.

---

## 9. Running and reruns

```matlab
Opt.rerun = 'off';  % full run: builds LHS matrices and returns them in LHS struct
Opt.rerun = 'on';   % reuse LHS matrices from a previous run (faster)
```

For a single run:

```matlab
[StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
```

For a parameter sweep (same base flow, different omega):

```matlab
omega_vec = Stab.omega_0 * linspace(0.8, 1.2, 5);
for k = 1:numel(omega_vec)
    Stab.omega_0 = omega_vec(k);
    if k == 1
        [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt);
    else
        Opt.rerun = 'on';
        [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt, LHS);
    end
end
```

When `Opt.rerun = 'on'`:
- LHS matrices (`MLA1`, `MLA2`, `MFD_diff`) are reused — no reassembly.
- LU factorizations are reused when `omega` and `beta` are unchanged (actuation/amplitude sweeps).
- LU factorizations are recomputed when `omega` or `beta` changes (frequency/wavenumber sweeps).

**With a resonator (`Opt.hr`).** The admittance sits in `MLA1` (the wall `v'`-row),
so a factorization belongs to one admittance. The assembly is still reused — the
previous mouth rows are reset and the new ones written — but the factorization is
redone whenever `Y`, the mouth or `uslip` changes. This covers a clean run followed
by a forced one, a sweep over `Y`, a moved mouth, and a return to a rigid wall:

```matlab
Opt.rerun = 'off';
[StabRes0, StabGrid, LHS] = main(Stab, StabGrid, Opt);   % clean: capture LHS
Opt.hr = hr;  Opt.rerun = 'on';
StabResA = main(Stab, StabGrid, Opt, LHS);               % forced: reuses assembly
```

Expect ~1.3x (assembly only) when the wall rows change, and ~3.8x when they do not
(4000x40 compressible: assembly ~20 s, factorization ~50 s). Every case reproduces
the from-scratch result bit-for-bit.

Indicative cost on the nonlinear Blasius TS case (M=5, Nx=800, Ny=40, `lu_mode='full'`):

| call | what runs | wall-clock |
|---|---|---|
| first (`rerun='off'`) | assembly + 6 LU factorizations + Picard | ~40 s |
| reuse (`rerun='on'`, amplitude sweep) | restored blocks + 6 cached factorizations, Picard only | ~6–8 s (~6× faster) |

The reused solve reproduces the from-scratch result to machine precision.

---

## 10. Output

The main output is `StabRes`:

| Field | Size | Description |
|---|---|---|
| `StabRes.A(j,i)` | `nf × nx` | Mode amplitude at each station |
| `StabRes.alpha(j,i)` | `nf × nx` | Streamwise wavenumber |
| `StabRes.u(j,y,i)` | `nf × ny × nx` | Streamwise velocity shape function |
| `StabRes.v(j,y,i)` | `nf × ny × nx` | Wall-normal velocity shape function |
| `StabRes.w(j,y,i)` | `nf × ny × nx` | Spanwise velocity shape function |
| `StabRes.p(j,y,i)` | `nf × ny × nx` | Pressure shape function |

Index convention: `y=1` is the freestream (`y=H`), `y=ny` is the wall (`y=0`).

### Saving results to disk

By default a run returns `StabRes` to the workspace but writes nothing to disk. Set `Opt.save_output` to a path to save the final (post-processed) results:

```matlab
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_blasius.mat');   % [] to skip
```

This writes `StabRes`, `StabGrid` and `Stab` to that `.mat` after post-processing. The benchmark callers set it by default to `output/StabRes_<case>.mat`; the `output/` folder is created automatically and its `.mat` files are git-ignored (like `input/` and `gridgen/output/`). Every run ends with a completion line:

```text
=== DeHNSSo finished — output saved to output/StabRes_blasius.mat ===
```

(or `=== DeHNSSo finished (output not saved; set Opt.save_output to save) ===` when the option is empty).

`Opt.save_output` is distinct from `Opt.save_checkpoint`, which dumps the *pre*-post-process state (full mode set with `.q` intact) for restarting via `Stab.IC='LOAD'`:

| Option | Saved | Contents | Use |
|---|---|---|---|
| `Opt.save_output` | after post-process | `StabRes, StabGrid, Stab` | final results for analysis/plotting |
| `Opt.save_checkpoint` | before post-process | `StabRes, StabGrid, Stab, Mode` (`.q` intact) | restart / amplitude continuation (`IC='LOAD'`) |

### Figures produced during a run

| Figure | Content |
|---|---|
| 1 | Inflow: base flow profiles, eigenvalue spectrum, selected shape function |
| 2 | Amplitude evolution during convergence iterations |
| 3 | Convergence history (relative error per mode, log scale) |
| 4 | Final amplitude development (optionally overlaid with PSE/DNS reference) |
| 5 | N-factor, growth rates, wavenumber, phase speed |
| 6 | Eigenfunctions at selected streamwise stations |

### Post-processing

The unified callers in `callers/benchmark/` already wire post-processing for you:

```matlab
% Reference for overlay (set ref_file = [] to skip)
ref_file  = fullfile(rootdir, 'literature', 'StabRes_<case>_<METHOD>.mat');
ref_label = '<METHOD>';

if strcmpi(Opt.plot, 'on')
    if ~isempty(ref_file) && isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label));
    else
        plot_amplitudes(StabGrid, StabRes);
    end
    plot_stability(StabGrid, StabRes, Opt);
end
```

Reference data shipped in `literature/`:

| File | Used by | Source |
|---|---|---|
| `StabRes_Blasius_NPSE.mat` | `blasius.m` | NPSE reference (Westerbeek 2024 Fig. 8) |
| `StabRes_SweptWing_flat_NPSE.mat` | `sweptwing_flat.m` | NPSE reference (Westerbeek 2024 Fig. 9) |
| `StabRes_SweptWing_Hump_AHLNS.mat` | `sweptwing_hump.m` | AHLNS reference (Westerbeek 2024 Fig. 10) |

To overlay a different reference, point `ref_file` at any `.mat` file with a `StabRes` struct in the format that `plot_amplitudes` expects (see `src/plot/plot_amplitudes.m`).

---

## 11. Common issues

| Problem | Likely cause | Fix |
|---|---|---|
| ILST does not find a physical mode | Domain too short/tall, wrong `omega`/`beta` | Check eigenvalue spectrum in Figure 1 |
| Solver does not converge | Amplitude too high, too few modes | Reduce `Stab.A0`, increase `Stab.M` |
| Memory error | Large `nx × ny`, `lu_mode = 'full'` | Use `lu_mode = 'auto'` with `lu_max_cache`; disable `parfor` to avoid worker copies |
| Results differ from PSE | Buffer too short, nonlinear effects | Increase `Opt.xb`, check `Stab.M` |
| Freestream BC oscillations | Domain truncated, Dirichlet BC | Use `'H_NM'` on affected components: e.g. `Opt.bc_top = {'H_NM','H_NM','H_NM'}` |
