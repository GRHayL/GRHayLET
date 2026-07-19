# One-Dimensional Magnetic Tests

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page records the visible magnetic keyword guard, left/right magnetic
states, rotation, `Bvec` and `Avec` assignments, coordinate expressions, and
checked-in parameter consumption. It does not establish magnetic divergence,
staggering at runtime, vector-potential reconstruction semantics, error-macro
termination, or physical validity.

## Summary

The function's guard accepts either `initial_Avec="GRHayLID"` or
`initial_Bvec="GRHayLID"`, while the subsequent point loop visibly assigns
both `Bvec` and `Avec` without a second keyword guard. Five Balsara arms provide
nonzero left/right states; equilibrium, sound wave, and shock tube share zero
states. Direction blocks rotate components and select the step coordinate.
Variables named `x_stag`, `y_stag`, and `z_stag` are assigned directly from
`x[index]`, `y[index]`, and `z[index]`, with no explicit half-cell offset, and
no checked-in GRHayLID source consumes `stagger_A_fields`.

## Mode Applicability

| Applicability | Visible magnetic behavior |
| --- | --- |
| HydroTest1D+Magnetic | Magnetic schedule guard and function cover Balsara, equilibrium, sound-wave, and shock-tube selections. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-MAG-01` | Function visibly accepts either selection keyword, constructs test states, and unconditionally assigns both `Bvec` and `Avec` in its loop. | visible-implementation | Function body | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` |
| `ID-MAG-02` | `initial_Avec` is extended with the `GRHayLID` keyword value. | declared | Keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` |
| `ID-MAG-03` | `initial_Bvec` is extended with the `GRHayLID` keyword value. | declared | Keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_Bvec` |
| `ID-MAG-04` | `stagger_A_fields` declares a default-enabled staggering control, but no checked-in source file reads it. | unresolved | Parameter declaration plus complete local search | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` |

## Details

### Selection guard and magnetic states

The opening condition calls `CCTK_VERROR` only when neither
`initial_Avec` nor `initial_Bvec` equals `GRHayLID`. Thus visible boolean
text accepts either choice. After that guard, no separate keyword condition
surrounds the point loop, whose body assigns all three components of both
`Bvec` and `Avec`. Whether the error macro terminates is external; the local
mismatch is between accepted textual selection and visible writes.

Before direction rotation, the dispatch assigns:

| Arm | Left `(Bx, By, Bz)` | Right `(Bx, By, Bz)` |
| --- | --- | --- |
| `Balsara1` | `(0.5, 1.0, 0.0)` | `(0.5, -1.0, 0.0)` |
| `Balsara2` | `(5.0, 6.0, 6.0)` | `(5.0, 0.7, 0.7)` |
| `Balsara3` | `(10.0, 7.0, 7.0)` | `(10.0, 0.7, 0.7)` |
| `Balsara4` | `(10.0, 7.0, 7.0)` | `(10.0, -7.0, -7.0)` |
| `Balsara5` | `(2.0, 0.3, 0.3)` | `(2.0, -0.7, 0.5)` |
| `equilibrium`, `sound wave`, `shock tube` | `(0.0, 0.0, 0.0)` | `(0.0, 0.0, 0.0)` |

### Rotation and field assignments

For y direction, code maps `(Bx, By, Bz)` to `(Bz, Bx, By)`. For z
direction it maps the tuple to `(By, Bz, Bx)`. Pointwise `step` selects the
matching x, y, or z coordinate; `step <= discontinuity_position` selects the
left tuple, otherwise the right tuple, for direct `Bvec` writes.

For `Avec`, code uses these side-specific expressions, where `(Bx, By, Bz)`
means either selected left or right tuple:

| Direction | `(Ax, Ay, Az)` visible expression |
| --- | --- |
| x | `(By*z_stag - Bz*y_stag, 0, Bx*y_stag)` |
| y | `(By*z_stag, Bz*x_stag - Bx*z_stag, 0)` |
| z | `(0, Bz*x_stag, Bx*y_stag - By*x_stag)` |

These formulas are visible assignments only. No local evidence establishes
their curl, gauge, staggering, or consistency with separately assigned
`Bvec`.

### Staggering declaration versus visible coordinates

A file comment says the routine assumes staggered vector potential. However,
the body assigns `x_stag=x[index]`, `y_stag=y[index]`, and
`z_stag=z[index]`; no half-cell offset appears in these expressions. Complete
search of checked-in GRHayLID files finds `stagger_A_fields` only in its
parameter declaration. This is checked-in non-consumption, not proof about
generated arguments, external frameworks, or runtime layout.

## Caveats

- Guard/write mismatch is tracked by
  [GID-0012](../contradictions.md#gid-0012); no execution claim follows.
- ThornGuide's statement that only `Avec` is set conflicts with visible
  `Bvec` writes; see [GID-0002](../contradictions.md#gid-0002).
- Declared staggering intent, absent consumer, direct coordinate-array
  expressions with no visible offset arithmetic, and source comment remain
  unresolved under
  [GID-0004](../contradictions.md#gid-0004).
- Cactus coordinate, vector-gridfunction, and indexing semantics are external.

## Sources

- [Magnetic implementation](../../../GRHayLID/src/1D_tests_magnetic_data.c)
- [Parameter declarations](../../../GRHayLID/param.ccl)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)

## Related Pages

- [One-Dimensional Hydrodynamic Tests](one-d-tests-hydro.md)
- [Parameters and Configurations](../integration/parameters-and-configurations.md)
- [GID-0002](../contradictions.md#gid-0002)
- [GID-0004](../contradictions.md#gid-0004)
- [GID-0012](../contradictions.md#gid-0012)
