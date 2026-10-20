# Sample run

A complete nonlinear test run of DeHNSSo-v2: the development of TS waves in a Blasius
boundary layer, the DeHNSSo-v1 test case used for the timing comparison of Table 4 in the
DeHNSSo-v2 paper. The run starts from the local-stability solution at the inflow and
converges the nonlinear interactions of the fundamental and its harmonics with the Picard
iteration. The result is compared with the DeHNSSo-v1 solution on the same grid.

| | |
|---|---|
| Grid | 800 streamwise × 40 Chebyshev points, R = 400 to 1067, identical to the grid of DeHNSSo-v1 |
| Disturbance | ω₀ = 0.034576, β = 0, M = 5, A₀ = 0.0025√2 |
| Solver | Picard iteration, `lu_mode = 'full'`, buffer from 85 % of the domain |

## Running it

From MATLAB:

```matlab
run sample/run_sample.m
```

The run needs MATLAB R2020b or later and no toolbox.

## Expected output

The run converges in 32 Picard iterations, saves the solution to `output/StabRes_sample.mat`,
and plots the amplitude of each mode against the DeHNSSo-v1 solution, as in
`reference/sample.pdf`. The two solutions agree to within 0.2 % of the peak of each mode.

## Files

| file | content |
|---|---|
| `run_sample.m` | sample run |
| `input/DeHNSSo_input_sample.mat` | stability grid and base flow (`StabGrid`), built by `literature/DeHNSSo-v1/TS_Blasius/gridgen_ts_blasius.m` from the base flow of DeHNSSo-v1 |
| `reference/StabRes_sample_v1.mat` | DeHNSSo-v1 solution: mode amplitudes and streamwise wavenumbers (`StabResV1`) and streamwise coordinate (`StabGridV1`) |
| `reference/sample.log` | console output of the run |
| `reference/sample.pdf` | mode amplitudes of the run and of DeHNSSo-v1 |
