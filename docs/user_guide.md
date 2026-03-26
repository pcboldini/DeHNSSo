# DeHNSSo — User Guide

**Version 2** | Delft Harmonic Navier-Stokes Solver

**Contact:** Pietro Carlo Boldini — P.C.Boldini@tudelft.nl — for questions, bug reports, or collaboration requests.

> **Note:** The accompanying paper is currently in preparation. Please contact the author for questions about the method or implementation.

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
cd baseflow/
caller        % produces BF_Blasius.mat in baseflow/output/
```

Then convert to DeHNSSo format:

```matlab
cd interpolation/
edit config.m   % set Re, grid size, domain height
main_interpol   % reads ../baseflow/output/, saves DeHNSSo_input.mat to ../stability/
```

For a curvilinear geometry (e.g., from OpenFOAM): see [README.md](../README.md).

### Step 2 — Configure and run

Open `callers/example_BlasiusBL.m`. The file is organized into sections:

| Section | What to set |
|---|---|
| **Spectral content** | `Stab.M`, `Stab.N`, `Stab.omega_0`, `Stab.beta_0` |
| **Inflow conditions** | `Stab.IC`, `Stab.A0`, `Stab.phaseRef` |
| **Boundary conditions** | Wall actuation (`Stab.bcw`) or freestream forcing (`Stab.bct`) |
| **Solver options** | Buffer, convergence criteria, LU cache, freestream BC type |
| **Run** | `Opt.rerun`, call to `main()` |
| **Post-processing** | Amplitude plots, reference data overlay |

Run the caller from the MATLAB command window:

```matlab
cd callers/
example_BlasiusBL
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

```matlab
Stab.A0 = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Stab.A0(ModeToModeNumber( 1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
Stab.A0(ModeToModeNumber(-1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);
```

Always set both the mode and its conjugate (index `-m`, `-n`) to ensure real physical fields. The amplitude is the max streamwise velocity perturbation, normalized by the freestream velocity.

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

Model a wall-normal velocity slot (e.g., to control TS waves):

```matlab
j_fund = ModeToModeNumber(1, 0, Stab.M, Stab.N);
j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
x_dom  = StabGrid.xun;
x0_act = x_dom(1) + 0.3*(x_dom(end) - x_dom(1));  % actuator centre
sig    = 0.01*(x_dom(end) - x_dom(1));              % actuator half-width
A_act  = 1e-4;
Stab.bcx             = x_dom;
Stab.bcw             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
Stab.bcw(2,:,j_fund) = A_act * exp(-(x_dom - x0_act).^2 / sig^2);
Stab.bcw(2,:,j_conj) = conj(Stab.bcw(2,:,j_fund));
```

`k=2` is the wall-normal velocity. Set `k=1` for streamwise, `k=3` for spanwise.

### Freestream forcing

Prescribe a perturbation at `y=H` (e.g., an acoustic wave or free-stream turbulence):

```matlab
Stab.bct(1,:,j_fund) = A_fst * ones(1, numel(x_dom));  % u at y=H
Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));
```

**Note:** for the `(0,0)` mean-flow distortion mode, the `v` freestream row carries the transpiration boundary condition and is not overridden by `bct`.

### Freestream BC type

```matlab
Opt.bc_top = 'dirichlet';  % u=v=w=0 at y=H (default; requires tall domain)
Opt.bc_top = 'neumann';    % d/dy=0  at y=H (use for truncated domains)
```

The same BC is applied in the ILST inflow computation.

---

## 7. Solver options

### Buffer zone

```matlab
Opt.buffer = 'on';    % amplitude damping (prevents end reflections)
Opt.buffer = 'para';  % parabolisation (zeroes streamwise diffusion, LPSE-like)
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

### LU factorisation cache

The LHS matrix is factorised once per mode and reused across iterations. Three caching strategies:

```matlab
Opt.lu_mode = 'full';  % cache all LU factorisations (fastest, most memory)
Opt.lu_mode = 'auto';  % cache up to lu_max_cache, evict oldest (bounded memory)
Opt.lu_mode = 'none';  % solve with backslash each time (slowest, least memory)
% Opt.lu_max_cache = 3;  % only relevant for 'auto' mode
```

---

## 8. Running and reruns

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
- LU factorisations are reused when `omega` and `beta` are unchanged (actuation/amplitude sweeps).
- LU factorisations are recomputed when `omega` or `beta` changes (frequency/wavenumber sweeps).

---

## 9. Output

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

### Figures produced during a run

| Figure | Content |
|---|---|
| 1 | Inflow: base flow profiles, eigenvalue spectrum, selected shape function |
| 2 | Amplitude evolution during convergence iterations |
| 3 | Convergence history (relative error per mode, log scale) |
| 4 | Final amplitude development (optionally overlaid with PSE/DNS reference) |

### Post-processing

```matlab
ref = struct('file', fullfile(root, 'stability', 'literature', 'StabRes_Blasius_NPSE.mat'), 'label', 'NPSE');
plot_amplitudes(StabGrid, StabRes, ref);
```

---

## 10. Common issues

| Problem | Likely cause | Fix |
|---|---|---|
| ILST does not find a physical mode | Domain too short/tall, wrong `omega`/`beta` | Check eigenvalue spectrum in Figure 1 |
| Solver does not converge | Amplitude too high, too few modes | Reduce `Stab.A0`, increase `Stab.M` |
| Memory error | Large `nx × ny`, `lu_mode = 'full'` | Use `lu_mode = 'auto'` with `lu_max_cache` |
| Results differ from PSE | Buffer too short, nonlinear effects | Increase `Opt.xb`, check `Stab.M` |
| Freestream BC oscillations | Domain truncated, Dirichlet BC | Switch to `Opt.bc_top = 'neumann'` |
