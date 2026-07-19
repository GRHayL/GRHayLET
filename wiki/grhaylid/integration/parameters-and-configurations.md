# Parameters and Configurations

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page owns all 21 restricted parameters declared locally, their declared
domains/defaults, and their visible consumers or checked-in non-consumption.
It also records the absence of checked-in GRHayLID parfiles. Shared HydroBase
and GRHayLib parameter semantics, runtime values, parser behavior, successful
configuration, and external files remain out of scope.

## Summary

`param.ccl` declares 21 restricted local parameters: two general magnetic
controls, four one-dimensional-test controls, three IsotropicGas values, ten
ConstantDensitySphere values, and two beta-equilibrium controls. Eleven real
setup values use `-1` as a forbidden default checked by `CHECK_PARAMETER`;
three interior velocities instead default to zero and accept any real value.
`stagger_A_fields` defaults to `yes` but has no visible consumer in the 14-file
source tree. No parfile, test declaration, or oracle is checked in below
`GRHayLID/**`.

## Mode Applicability

| Applicability | Local controls |
| --- | --- |
| Common | `initialize_magnetic_quantities` and `stagger_A_fields`; only the former is visibly consumed. |
| HydroTest1D | `initial_data_1D`, `shock_direction`, `discontinuity_position`, and `wave_amplitude`. |
| HydroTest1D+Magnetic | General and one-dimensional controls apply; magnetic body also reads the test, direction, and discontinuity controls. |
| IsotropicGas | Three required nonnegative thermodynamic parameters. |
| ConstantDensitySphere | Radius, interior/exterior thermodynamics, and three interior velocity components. |
| BetaEquilibrium | Boolean schedule switch and required nonnegative temperature. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-PAR-01` | `initialize_magnetic_quantities` is a boolean, defaults to `yes`, and guards magnetic scheduling. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=initialize_magnetic_quantities` |
| `INT-PAR-02` | `stagger_A_fields` is a boolean, defaults to `yes`, and declares +1/2 staggering intent, but has no visible consumer in the checked-in source tree. | unresolved | Local declaration plus full-tree search | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` |
| `INT-PAR-03` | `initial_data_1D` is a keyword with eight admitted values and default `Balsara1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `INT-PAR-04` | `shock_direction` admits `x`, `y`, or `z` and defaults to `x`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=shock_direction` |
| `INT-PAR-05` | `discontinuity_position` accepts any real value and defaults to `0.0`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=discontinuity_position` |
| `INT-PAR-06` | `wave_amplitude` is nonnegative and defaults to `1.0e-3`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` |
| `INT-PAR-07` | `IsotropicGas_rho` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_rho` |
| `INT-PAR-08` | `IsotropicGas_Y_e` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_Y_e` |
| `INT-PAR-09` | `IsotropicGas_temperature` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_temperature` |
| `INT-PAR-10` | `ConstantDensitySphere_sphere_radius` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_sphere_radius` |
| `INT-PAR-11` | `ConstantDensitySphere_rho_interior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_interior` |
| `INT-PAR-12` | `ConstantDensitySphere_Y_e_interior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_interior` |
| `INT-PAR-13` | `ConstantDensitySphere_T_interior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_interior` |
| `INT-PAR-14` | `ConstantDensitySphere_vx_interior` accepts any real value and defaults to `0`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vx_interior` |
| `INT-PAR-15` | `ConstantDensitySphere_vy_interior` accepts any real value and defaults to `0`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vy_interior` |
| `INT-PAR-16` | `ConstantDensitySphere_vz_interior` accepts any real value and defaults to `0`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vz_interior` |
| `INT-PAR-17` | `ConstantDensitySphere_rho_exterior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_exterior` |
| `INT-PAR-18` | `ConstantDensitySphere_Y_e_exterior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_exterior` |
| `INT-PAR-19` | `ConstantDensitySphere_T_exterior` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_exterior` |
| `INT-PAR-20` | `impose_beta_equilibrium` is a boolean and defaults to `no`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=impose_beta_equilibrium` |
| `INT-PAR-21` | `beq_temperature` is nonnegative with forbidden/default value `-1`. | declared | Local declaration | `ccl:GRHayLID/param.ccl#parameter=beq_temperature` |
| `INT-PAR-22` | `CHECK_PARAMETER(name)` visibly compares its argument with `-1` and expands to a `CCTK_VERROR` call on equality. | visible-implementation | Local macro definition | `macro:GRHayLID/src/GRHayLID.h#name=CHECK_PARAMETER` |

## Details

### Exact 21-parameter inventory

| Parameter | Type; declared range/values; default | Declared meaning | Visible consumer or checked-in absence |
| --- | --- | --- | --- |
| `initialize_magnetic_quantities` | `CCTK_BOOLEAN`; `yes` | Whether Avec should be initialized | Nested schedule guard for magnetic 1D data |
| `stagger_A_fields` | `CCTK_BOOLEAN`; `yes` | Whether Avec is staggered +1/2 following named convention | No match outside `param.ccl`; [GID-0004](../contradictions.md#gid-0004) |
| `initial_data_1D` | `KEYWORD`; equilibrium, sound wave, shock tube, Balsara1-5; `Balsara1` | One-dimensional test selection | Both 1D C bodies; hydro dispatch has seven live arms and one commented arm; `initial_hydro` description instead names `test_1D_initial_data` ([GID-0015](../contradictions.md#gid-0015)) |
| `shock_direction` | `KEYWORD`; x/y/z; `x` | Discontinuity direction | Both 1D C bodies rotate state and select coordinate |
| `discontinuity_position` | `CCTK_REAL`; `*:*`; `0.0` | Shock location | Both 1D C bodies compare selected coordinate |
| `wave_amplitude` | `CCTK_REAL`; `0.0:*`; `1.0e-3` | Sound-wave amplitude | Hydro loop sine expression ordered after the invalid-name `CCTK_VERROR` arm; reachability unverified |
| `IsotropicGas_rho` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Gas density | Sentinel macro, EOS call, uniform density write |
| `IsotropicGas_Y_e` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Electron fraction | Sentinel macro, EOS call, uniform electron-fraction write |
| `IsotropicGas_temperature` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Temperature | Sentinel macro, EOS call, uniform temperature write |
| `ConstantDensitySphere_sphere_radius` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Sphere radius | Sentinel macro and `r[index]` branch |
| `ConstantDensitySphere_rho_interior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Interior density | Sentinel macro, interior EOS call/write |
| `ConstantDensitySphere_Y_e_interior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Interior electron fraction | Sentinel macro, interior EOS call/write |
| `ConstantDensitySphere_T_interior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Interior temperature | Sentinel macro, interior EOS call/write |
| `ConstantDensitySphere_vx_interior` | `CCTK_REAL`; `*:*`; `0` | Interior x velocity | Interior x-velocity assignment; no sentinel macro |
| `ConstantDensitySphere_vy_interior` | `CCTK_REAL`; `*:*`; `0` | Interior y velocity | Interior y-velocity assignment; no sentinel macro |
| `ConstantDensitySphere_vz_interior` | `CCTK_REAL`; `*:*`; `0` | Interior z velocity | Interior z-velocity assignment; no sentinel macro |
| `ConstantDensitySphere_rho_exterior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Exterior density | Sentinel macro, exterior EOS call/write |
| `ConstantDensitySphere_Y_e_exterior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Exterior electron fraction | Sentinel macro, exterior EOS call/write |
| `ConstantDensitySphere_T_exterior` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Exterior temperature | Sentinel macro, exterior EOS call/write |
| `impose_beta_equilibrium` | `CCTK_BOOLEAN`; `no` | Impose beta equilibrium | Schedule guard; no direct C read |
| `beq_temperature` | `CCTK_REAL`; `0:*`, `-1` forbidden; `-1` | Temperature used for beta equilibrium | Bounds comparison, sentinel macro, helper calls, point writes |

### Sentinel pattern

The macro `CHECK_PARAMETER(par)` visibly compares `par == -1` and invokes
`CCTK_VERROR` with a parfile-setting message. IsotropicGas applies it to all
three setup values. ConstantDensitySphere applies it to radius and six
thermodynamic values, but not the three velocities. BetaEquilibrium applies it
to `beq_temperature` after a table-bounds comparison; that ordering is
[GID-0007](../contradictions.md#gid-0007). Macro expansion does not establish
termination semantics.

### Checked-in configuration absence

The registered 14-file tree contains README, documentation, CCL/build inputs,
a header, and six C units. It contains no `.par`, `test.ccl`, test directory,
or numeric oracle. This is a checked-in coverage gap, not evidence that no
external configuration or testing exists. ThornGuide's nonexistent
`compute_entropy` parameter is tracked in [GID-0001](../contradictions.md#gid-0001).
Separately, the `initial_hydro="HydroTest1D"` description points to
`test_1D_initial_data`, but this inventory's declared selector is
`initial_data_1D`; see [GID-0015](../contradictions.md#gid-0015).

## Caveats

- Declared defaults are not observed runtime selections.
- CCL ranges and keyword values do not prove external parser behavior.
- A visible consumer does not prove that its containing schedule runs.
- Non-consumption is scoped to the complete checked-in 14-file tree and must
  be reconciled if files are added.
- The HydroTest1D extension description and restricted selector use different
  parameter names; external alias behavior is unverified under
  [GID-0015](../contradictions.md#gid-0015).

## Sources

- [Local parameter declarations](../../../GRHayLID/param.ccl)
- [Sentinel macro](../../../GRHayLID/src/GRHayLID.h)

## Related Pages

- [HydroBase Keyword Extensions](hydrobase-keyword-extensions.md)
- [GRHayLib Contract](grhaylib-contract.md)
- [Coverage Gaps](../validation/coverage-gaps.md)
- [GID-0001: entropy-control mismatch](../contradictions.md#gid-0001)
- [GID-0004: unconsumed staggering control](../contradictions.md#gid-0004)
- [GID-0007: beta-temperature check ordering](../contradictions.md#gid-0007)
- [GID-0015: HydroTest1D parameter-name mismatch](../contradictions.md#gid-0015)
