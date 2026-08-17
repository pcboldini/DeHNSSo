# DeHNSSo — User Guide

**Version 1.2** | Delft Harmonic Navier-Stokes Solver

---

## 1. What DeHNSSo does

DeHNSSo solves the nonlinear stability problem for laminar boundary layers. It computes the amplitude evolution of wave-like disturbances (Tollmien-Schlichting waves, cross-flow vortices) as they propagate downstream and interact nonlinearly.

The equations are the incompressible Navier-Stokes equations, written for perturbations on top of a steady base flow. The perturbations are expanded in a Fourier series in time and spanwise direction (harmonic balancing), and the resulting coupled system is solved iteratively until convergence.

---

## 2. Prerequisites

- MATLAB R2020b or later (no additional toolboxes required)
- A base-flow input file `DeHNSSo_input.mat` containing the `StabGrid` struct

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
blasius       % reads bf_blasius.mat from baseflow/output/benchmark/, writes DeHNSSo_input.mat
```

For other supported cases (`blasius_csv`, `sweptwing_flat`, `sweptwing_hump`, `m3j` / OpenFOAM): see [README.md](../README.md).

### Step 1b — Optional streamwise clustering (all cases)

By default the gridgen scripts distribute `nx` stations uniformly in ξ (wall arc length for curvilinear cases, `x` for the flat plate). For cases with a localized feature (hump, leading-edge curvature, roughness element) — or simply to concentrate stations where the disturbance grows fastest — you can cluster stations around one or more centers via a Gaussian step-size modulation. Configured in the per-case gridgen script via three `params` fields:

| Parameter | Meaning |
|---|---|
| `params.mug` | Clustering center(s) in physical x (scalar or vector for multiple centers) |
| `params.sig` | Gaussian width as a fraction of the half-domain (default 0.2). Scalar broadcasts to all centers; vector must match `mug` length |
| `params.ag` | Step-size modulation at the peak: `<1` denser, `=1` no refinement, `>1` sparser (default 1). Scalar or vector |

The same three fields work for both the curvilinear (`resample_curvilinear.m`) and flat-plate (`resample_cartesian.m`) paths, which share one implementation in `gridgen/src/resample/cluster_xi.m`.

Algorithm (`gridgen/src/resample/cluster_xi.m`):

1. Map each `mug` into normalized ξ ∈ [-1, 1].
2. At every candidate ξ, compute `G_k(ξ) = exp(−½·((ξ − μ_n)/σ)²)` for each hump `k`.
3. Local modulation: `fg(ξ) = 1 − (1 − ag) · max_k G_k(ξ)` → `fg = ag` at each peak.
4. Build the new ξ grid by integrating `∫ dξ / fg` and renormalizing — denser where `fg` is small.

Example (heritage `sweptwing_hump.m`):

```matlab
x_m = 0.1837410879 / 2.1394e-4;   % physical hump center / lref ≈ 859
params.mug = [x_m];               % single hump
params.sig = [0.2];               % Gaussian half-width fraction
params.ag  = [0.5];               % step size halved at the peak
```

Multi-hump example:

```matlab
params.mug = [x_m1, x_m2];        % two humps
params.sig = [0.2, 0.15];         % per-hump widths
params.ag  = [0.5, 0.7];          % per-hump strengths
```

Leave `params.mug = []` or `params.ag = 1` everywhere to fall back to uniform spacing (wall arc length for curvilinear, `x` for the flat plate).

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

Per-case stability scripts live in `callers/benchmark/` (`blasius.m`, `sweptwing_flat.m`, `sweptwing_hump.m`, `m3j.m`). All four share an identical sectioned layout — only case-specific values differ:

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

| case | v1 caller (`A0 = |û|`) | v2 caller (`A0 = 2|û|`) |
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

A complete working example is `tests/warmstart/test_load_warmstart.m` — `callers/benchmark/blasius.m` plus a `restart` flag. On the Blasius benchmark (800×40, M=5), continuing `A0 = 0.00125√2 → 0.002` reaches the same solution as a cold run at 0.002 (`max(A)` agreeing to 7e-05) in 43 iterations instead of 57.

The same mechanism serves parameter sweeps: restart at a new `omega_0`, `beta_0` or grid-independent parameter from a nearby converged state. The checkpoint must come from a run with the same `M`/`N` and grid.

### No inflow disturbance (`Stab.IC = 'ZERO'`)

`Stab.IC='ZERO'` seeds **no** disturbance at the inflow — the entire field is generated by wall/freestream forcing (`Stab.bcw` / `Stab.bct`). Use it when the quantity of interest is the receptivity of a wall source *in isolation*: a synthetic-jet or plasma actuator, a roughness element, or the resonator impedance BC, with **no incoming instability wave**. It works in both incompressible and compressible mode, and needs no `Stab.A0_fund` (there is no inflow mode to scale).

```matlab
Stab.IC  = 'ZERO';
Stab.bcx = x;
Stab.bcw = zeros(3, numel(x), (2*Stab.N+1)*(2*Stab.M+1));
Stab.bcw(2, :, j_fund) = vjet;   % e.g. a wall-normal velocity actuator on the fundamental
```

The inflow plane is clamped to zero, so the active modes come from the wall forcing alone. `inflow` warns (`inflow:zeroNoForcing`) if there is no forcing at all, since the field would then be identically zero. A minimal example is `tests/test_zero_ic.m` — a wall synthetic jet driving a Blasius BL with no incoming TS wave.

---

## 6. Boundary conditions

By default all BCs are homogeneous (zero perturbation at wall and freestream). You can prescribe inhomogeneous BCs to model actuation or external disturbances.

### Arrays

```matlab
Stab.bcx          % (1 x nx_bc)       x-locations of BC profile
Stab.bcw(k,ix,j)  % (3 x nx_bc x nf)  wall BC:       k=1 u, k=2 v, k=3 w
Stab.bct(k,ix,j)  % (3 x nx_bc x nf)  freestream BC: k=1 u, k=2 v, k=3 w
```

### Wall actuation (blowing/suction)

Wall forcing follows a **grid-first** workflow: you describe what is at the wall **in the gridgen script**, gridgen refines the mesh to resolve it and stores the profile on the grid, and the caller applies it in one line. Location, width and shape are defined **once**, in gridgen; only the strength (amplitude, phase) lives in the caller. The reference pair is `gridgen/benchmark/synthetic_jet.m` + `callers/benchmark/synthetic_jet.m`.

This works on **any grid** — cartesian, curvilinear, and unstructured. On a flat plate `xc`/`x_data` are the physical `x`; on a curved or unstructured wall (swept-wing, hump, M3J/CSV) they are still physical `x` — gridgen samples the profile at the physical wall-x per station, so it lands geometrically correctly even though the grid runs in arc-length.

**Step 1 — gridgen: define the forcing.** `params.forcing` is either a parameterized actuator or a measured/CFD distribution; gridgen sizes the clustering to it automatically:

```matlab
% (a) raised-cosine actuator, center xc, half-width w (physical x)
params.forcing = struct('xc', xc, 'w', w, 'shape', 'cosine');

% (b) an arbitrary v-distribution — any shape; xc/w inferred from its support
params.forcing = struct('x', x_data, 'v', v_data);   % v_data may be complex (phase)

params.ag = 0.3;   % ~3x denser over the forcing (density knob)
```

For **several forcings**, give `xc` a vector (many same-type actuators) or pass a cell array (mixed actuators and distributions); gridgen refines for all of them and stores one profile each:

```matlab
params.forcing = struct('xc', [xc1, xc2, xc3], 'w', 18);          % w scalar broadcasts, or a vector
params.forcing = { struct('xc', xc1, 'w', 18), struct('x', x_data, 'v', v_data) };
```

gridgen stores `StabGrid.forcing` (the unit-peak profile(s) sampled on the refined grid) plus `StabGrid.cluster.xc/.w`. The distribution's `v` values need only encode the **shape**; the amplitude is set in the caller. Note the fixed station budget is shared across all centers, so raise `params.n_xi_new` when refining several spots.

**Step 2 — caller: apply the strength.** `add_wall_forcing` reads `StabGrid.forcing`, scales it, and does all the `Stab.bcw` setup — component, mode index, `Stab.bcx`, sizing, and the conjugate mode (added automatically for nonlinear runs):

```matlab
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 1e-4));   % xc/w/shape from gridgen
```

`spec` fields: `A` (amplitude, physical `2|û|` like `Stab.A0`, required), optional `phi` (phase), `comp` (1 u / 2 v / 3 w, default 2), `mode` (`[m n]`, default `[1 0]`), and `which` to pick one forcing when gridgen stored several. The forcing is **added** to any existing `Stab.bcw`, so call once per source to stack several:

```matlab
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 2e-4, 'which', 1));
Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('A', 5e-4, 'which', 2, 'phi', pi/2));
```

`which` is required once the grid holds more than one forcing (the call errors otherwise) and unnecessary for a single one. You can also override the profile here for a quick test — pass `x`/`v` (distribution) or `xc`/`w` (actuator) — but the mesh is only refined where gridgen put it, so an override off the refined patch is coarser (`add_wall_forcing` warns if it lands on too few stations).

Because the mesh is static, the width is capped by what you refined for: a narrower forcing than `w_act` stays resolved; a **wider** one, or a moved center, needs a new grid. Size the gridgen `w`/support for the widest case you intend to run, then sweep amplitude/phase freely with no re-mesh.

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

```matlab
Opt.hr = struct( ...
    'Y',     Ynd, ...          % non-dim resonator admittance (a fixed value in the caller)
    'mask',  in_neck, ...      % logical, nx: the mouth streamwise stations
    'x',     StabGrid.xun, ... % optional: grid-exact flux-preserving edge taper
    'taper', 0.1, ...          % optional: raised-cosine mouth edge taper (0 = hard)
    'uslip', false);           % optional: free-slip u' over the mouth
```

Cluster the streamwise mesh at the neck in gridgen exactly as for a wall actuator (`params.forcing = struct('xc', Xc, 'w', w)`), so the mouth carries enough stations for the impedance row. The reference pair is `gridgen/benchmark/resonator.m` + `callers/benchmark/resonator.m` (validated against Michelis 2023).

### Freestream forcing

Prescribe a perturbation at `y=H` (e.g., an acoustic wave or free-stream turbulence):

```matlab
Stab.bct(1,:,j_fund) = A_fst * ones(1, numel(x_dom));  % u at y=H
Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));
```

**Note:** for the `(0,0)` mean-flow distortion mode, the `v` freestream row carries the transpiration boundary condition and is not overridden by `bct`.

### Freestream BC type

`Opt.bc_top` is a cell array specifying the BC type per velocity component {u, v, w}:

```matlab
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};   % default: all zero at y=H

% BC types:
%   'H_DR'  = homogeneous Dirichlet  (q = 0)
%   'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%   'H_NM'  = homogeneous Neumann    (dq/dy = 0)
```

Pressure has no BC at the freestream (determined by the equations). The MFD `v` component is always free at the freestream (displacement velocity).

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

### Picard acceleration & Newton fallback

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

**There is no fixed-amplitude Newton.** The old backends (JFNK, `direct`, `matlab_direct`, `lu_krylov`) and `solve_newton.m` have been removed: a cold-start solve at a fixed high amplitude does not converge, which is precisely what the continuation exists to solve. Setting those `Opt.solver` values now raises an error. `Opt.solver='newton'` means `src/newton/solve_newton.m`.

### Pseudo-arclength continuation (`solve_newton`)

For high-amplitude K-type cases (Scholten), a direct solve at the target amplitude does not converge from a cold start. `src/newton/solve_newton.m` continues the solution along the branch in the inflow amplitude `A0_fund` instead of trying to reach the target in one go. Call it directly (it does its own `init`/`inflow`/`setup_boundary`/`build_ctx`) rather than through `main`:

```matlab
[StabRes, info] = solve_newton(Stab, StabGrid, Opt);
```

The predictor is a Keller step (tangent plus secant curvature); the corrector is a modified Newton on the full coupled real-form Jacobian, factorized with MUMPS. The factor is **held across steps** (the corrector then needs only back-solves) and refactored when it ages or the corrector stalls, since one factorization costs about 230× a back-solve and the state moves little between steps.

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
```

`newton_reuse` is the reuse policy and the main performance knob: `newton_reuse=1` refactors every step (most robust), larger values trade factorizations for more corrector back-solves. Which wins is case-dependent — measure it, and only trust a timing from a log whose first step is `STEP 1` (a resumed run's reported total covers only its own segment).

`lambda` enters the residual only through the inflow Dirichlet BC, hence affinely: `F(q,λ) = residual_hns(q,ctx0) + (λ−λ0)·F_λ`. `F_λ` is computed once by finite difference and its linearity is checked and reported at startup.

**Memory.** The coupled Jacobian is factorized whole (2N ≈ 1.7M for M=4 on a 1050×40 grid). Peak RSS during factorization reaches ~150 GB at M=4 and falls back to ~35 GB once the factor is formed — each refactor re-spikes. M=4 needs a workstation; a 16 GB laptop is limited to roughly M ≤ 3.

### LU mode

```matlab
Opt.lu_mode = 'full';         % cache every LU factor — fastest reuse, unbounded memory (default)
Opt.lu_mode = 'auto';         % LRU-cached LUs, bounded memory (Opt.lu_max_cache, default 3)
Opt.lu_mode = 'seq';          % one LU at a time: factor, solve its modes, discard — minimal memory
Opt.lu_mode = 'none';         % plain backslash each time, no caching (slowest, baseline)
Opt.lu_mode = 'lapack_band';  % banded LU MEX (zgbtrf/zgbtrs) with per-mode factor caching —
                              %   fastest single solve (~3x); one-time MEX build (see below)
Opt.lu_mode = 'block_thomas'; % dual-Schur block-Thomas (Paredes 2024) — COMPRESSIBLE operator only
% Opt.lu_max_cache = 3;       % cache size for 'auto' (and the banded / block-Thomas factor caches)
```

`full`, `auto`, `seq`, and `none` all use MATLAB's sparse LU and differ only in how many factors they keep in memory: `full` keeps all of them, `auto` keeps the last `lu_max_cache`, `seq` keeps one at a time, `none` keeps none. `lapack_band` and `block_thomas` replace that sparse LU with a banded or block factorization of the same operator, then cache the factors the same way.

**`block_thomas`** is for the **compressible** operator only. The incompressible continuity rows have a zero diagonal, which makes the unpivoted block factorization silently wrong (verified ~2× error), so the solver raises an error if you select it on an incompressible run. Sub-options:

```matlab
Opt.bt_np    = 2;      % number of partitions (default 2)
Opt.bt_ooc   = 'on';   % out-of-core: stream the factor to disk, minimal RAM (default off)
Opt.lu_cache = 'off';  % one-shot factor+solve, no caching (default: cache like the other modes)
```

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
