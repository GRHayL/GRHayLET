# One-Dimensional Hydrodynamic Tests

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page records visible hydrodynamic dispatch, state constants, coordinate
selection, velocity rotation, pointwise assignments, and EOS-library call
sites for one-dimensional test data. It does not establish error-macro
termination, runtime reachability, EOS semantics, physical correctness, or
numerical validity.

## Summary

The visible pre-loop dispatch has seven uncommented arms: Balsara1 through
Balsara5, equilibrium, and shock tube. A sound-wave arm is present only inside
a block comment. The final invalid-name `else` precedes a point loop that
contains a sound-wave sine branch. Step data select left or right states at
`discontinuity_position`, with coordinate and velocity rotation controlled by
`shock_direction`.

## Mode Applicability

| Applicability | Visible hydrodynamic behavior |
| --- | --- |
| HydroTest1D | Seven uncommented step-state arms plus the separately located, visibly dead sound-wave dispatch arm and later loop-side sine branch. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-HYDRO-01` | Visible function order places seven uncommented state arms and a final error arm before rotation and the point loop. | visible-implementation | Function body | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` |
| `ID-HYDRO-02` | Declared `sound wave` choice has a commented-out pre-loop dispatch arm, while a loop-side sine branch appears after the final error arm. | unresolved | Function body and commented dispatch text | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` |
| `ID-HYDRO-03` | `initial_data_1D` declares eight choices, including `sound wave`, and defaults to Balsara1. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `ID-HYDRO-04` | `shock_direction` declares x, y, and z choices and defaults to x. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=shock_direction` |
| `ID-HYDRO-05` | `discontinuity_position` declares an unrestricted real and defaults to zero. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=discontinuity_position` |
| `ID-HYDRO-06` | `wave_amplitude` declares a nonnegative real and defaults to `1.0e-3`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` |
| `ID-HYDRO-07` | Function calls `ghl_hybrid_find_polytropic_index`; its lookup semantics and returned index are external. | out-of-scope | Call site within function | `c:GRHayLID/src/1D_tests_hydro_data.c#call=ghl_hybrid_find_polytropic_index?function=GRHayLID_1D_tests_hydro_data` |

## Details

### EOS guard and visible dispatch

The function first compares `EOS_type` with `Tabulated` and calls
`CCTK_ERROR` when equal. Its message says the data is defined for "hybrid or
ideal fluid EOS" and calls the ideal-fluid EOS the standard comparison. This
is narrower local evidence than any general EOS behavior: only one keyword
comparison and one message are visible, and `CCTK_ERROR` termination semantics
remain external.

The seven uncommented pre-loop arms assign these left/right constants before
direction rotation:

| Arm | Left `(rho, press; vx, vy, vz)` | Right `(rho, press; vx, vy, vz)` |
| --- | --- | --- |
| `Balsara1` | `(1.0, 1.0; 0, 0, 0)` | `(0.125, 0.1; 0, 0, 0)` |
| `Balsara2` | `(1.0, 30.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |
| `Balsara3` | `(1.0, 1000.0; 0, 0, 0)` | `(1.0, 0.1; 0, 0, 0)` |
| `Balsara4` | `(1.0, 0.1; 0.999, 0, 0)` | `(1.0, 0.1; -0.999, 0, 0)` |
| `Balsara5` | `(1.08, 0.95; 0.4, 0.3, 0.2)` | `(1.0, 1.0; -0.45, -0.2, 0.2)` |
| `equilibrium` | `(1.0, 1.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |
| `shock tube` | `(2.0, 2.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |

These are seven live textual arms, not eight dataset branches. The declaration
has an eighth choice, `sound wave`, whose corresponding pre-loop arm is
inside `/* ... */`. That comment says the case is handled in the loop. The
uncommented chain therefore has no sound-wave arm before its final
`CCTK_VERROR` call.

### Direction and pointwise selection

Pre-loop velocity rotation maps an x-oriented state as follows:

- y direction maps `(vx, vy, vz)` to `(vz, vx, vy)`.
- z direction maps `(vx, vy, vz)` to `(vy, vz, vx)`.
- x direction leaves assigned components unchanged.

Inside the loop, `step` selects `x`, `y`, or `z` using the same direction
parameter. Non-sound-wave points use left state when
`step <= discontinuity_position` and right state otherwise.

The later sound-wave branch instead assigns `rho=1`, `press=1`, x-velocity
`wave_amplitude * sin(M_PI * step)`, and zero y/z velocity. Its pressure line
contains `should add kinetic energy here`. This page records only visible
ordering and text. It does not infer whether control reaches that loop after
the earlier `CCTK_VERROR` call.

### Internal-energy expression and delegated call

After either pointwise branch, code reads
`ghl_eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(ghl_eos,
rho[index])]` into `Gamma`, then assigns
`eps = press / (rho * (Gamma - 1))`. Array contents, lookup behavior, valid
indices, EOS consistency, and arithmetic validity are GRHayLib behavior or
runtime properties and remain out of scope.

## Caveats

- "Live arm" means uncommented source arm only, not observed runtime
  execution.
- Sound-wave control-order mismatch remains open under
  [GID-0005](../contradictions.md#gid-0005); no termination or reachability
  conclusion is made.
- Sound-wave pressure comment records acknowledged incomplete intent under
  [GID-0006](../contradictions.md#gid-0006), not proof of a numerical defect.
- Cylindrical explosion documentation mismatch is tracked by
  [GID-0003](../contradictions.md#gid-0003).
- README and guard/error wording differ on supported EOS terminology; see
  [GID-0013](../contradictions.md#gid-0013).

## Sources

- [Hydrodynamic implementation](../../../GRHayLID/src/1D_tests_hydro_data.c)
- [Parameter declarations](../../../GRHayLID/param.ccl)

## Related Pages

- [One-Dimensional Magnetic Tests](one-d-tests-magnetic.md)
- [GRHayLib Contract](../integration/grhaylib-contract.md)
- [GID-0003](../contradictions.md#gid-0003)
- [GID-0005](../contradictions.md#gid-0005)
- [GID-0006](../contradictions.md#gid-0006)
- [GID-0013](../contradictions.md#gid-0013)
