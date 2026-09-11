# CFI over a step in a swept-wing boundary layer — NOT AVAILABLE IN v2

This v1 case cannot currently be reproduced. It is kept here with its reference data so
the gap is visible and the case can be restored when the capability returns.

**v1 caller:** `Callers/CFI over a step in a swept-wing BL/Example_Caller_CFI_step_in_SweptWing_BL.m`
**Reference data:** `StabRes_SweptWing_step_DNS.mat` (DNS)

## Why it is missing

v1 did not mesh the step. It kept the wall flat and imposed the step through an embedded
boundary method:

```matlab
xw = linspace(0.1464,0.2125,5000)/BF.lref;
yw = 0*xw;            % wall is flat
Grid.wall = [xw;yw];
Grid.StepX = 0.1837/BF.lref;
Grid.StepH = 7.4710e-04/BF.lref;
Grid.StepType = "FFS";
```

with the note *"Current implementation features EBM. Supply flat wall data and separate
step specifications"*. That is how a sharp corner was handled without meshing it.

DeHNSSo-v2 has no `StepX`, `StepH`, `StepType` or embedded-boundary code in `src/` or
`gridgen/`. Its body-fitted path cannot take their place: the compressible operator needs
the second derivatives of the mapping (`build_stab_grid.m` builds `xixx`/`etaxx` and
divides by `J^3`), a corner makes those unbounded, and wall normals from either side of a
step cross so the mapping folds. `solve_elliptic_mesh.m` only warns on a negative
Jacobian, so a folded mesh runs and gives a wrong answer rather than failing.

## To restore it

Re-add the embedded boundary method, then port the v1 caller. The same machinery also
covers the second v1 step case,
`Callers/TS-Waves in a Blasius BL/Example_Caller_Blasius_BL_TS_step.m`.
