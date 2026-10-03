# Parameters and Configurations

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

param.ccl declares 23 restricted local parameters, including the added
native entropy-proxy opt-in and beta residual tolerance. The original -1
required-input convention remains, with finite/domain checks in the bodies.
Staggering now has a visible consumer. No local parfile/test/oracle is checked in.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | Shared selector checks and explicit precision/representation contracts. |
| HydroTest1D | Test/direction/discontinuity/amplitude controls. |
| HydroTest1D+Magnetic | Enabled both-selector magnetics and staggering. |
| IsotropicGas | Finite effective-bounds primitive triple. |
| ConstantDensitySphere | Two triples, radius, and subluminal Cartesian velocity. |
| BetaEquilibrium | Temperature and residual tolerance. |
| Entropy/Hybrid | Explicit native-proxy opt-in. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `PARAMETERS-AND-CONFIGURATIONS-01` | initialize_magnetic_quantities has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initialize_magnetic_quantities` |
| `PARAMETERS-AND-CONFIGURATIONS-02` | stagger_A_fields has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` |
| `PARAMETERS-AND-CONFIGURATIONS-03` | initial_data_1D has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `PARAMETERS-AND-CONFIGURATIONS-04` | shock_direction has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=shock_direction` |
| `PARAMETERS-AND-CONFIGURATIONS-05` | discontinuity_position has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=discontinuity_position` |
| `PARAMETERS-AND-CONFIGURATIONS-06` | wave_amplitude has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` |
| `PARAMETERS-AND-CONFIGURATIONS-07` | IsotropicGas_rho has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_rho` |
| `PARAMETERS-AND-CONFIGURATIONS-08` | IsotropicGas_Y_e has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_Y_e` |
| `PARAMETERS-AND-CONFIGURATIONS-09` | IsotropicGas_temperature has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_temperature` |
| `PARAMETERS-AND-CONFIGURATIONS-10` | ConstantDensitySphere_sphere_radius has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_sphere_radius` |
| `PARAMETERS-AND-CONFIGURATIONS-11` | ConstantDensitySphere_rho_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-12` | ConstantDensitySphere_Y_e_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-13` | ConstantDensitySphere_T_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-14` | ConstantDensitySphere_vx_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vx_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-15` | ConstantDensitySphere_vy_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vy_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-16` | ConstantDensitySphere_vz_interior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vz_interior` |
| `PARAMETERS-AND-CONFIGURATIONS-17` | ConstantDensitySphere_rho_exterior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_exterior` |
| `PARAMETERS-AND-CONFIGURATIONS-18` | ConstantDensitySphere_Y_e_exterior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_exterior` |
| `PARAMETERS-AND-CONFIGURATIONS-19` | ConstantDensitySphere_T_exterior has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_exterior` |
| `PARAMETERS-AND-CONFIGURATIONS-20` | impose_beta_equilibrium has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=impose_beta_equilibrium` |
| `PARAMETERS-AND-CONFIGURATIONS-21` | beq_temperature has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=beq_temperature` |
| `PARAMETERS-AND-CONFIGURATIONS-22` | allow_native_entropy_proxy has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=allow_native_entropy_proxy` |
| `PARAMETERS-AND-CONFIGURATIONS-23` | beq_residual_tolerance has a local parameter declaration. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=beq_residual_tolerance` |

## Details

| Parameter | Declared type | Default |
| --- | --- | --- |
| initialize_magnetic_quantities | CCTK_BOOLEAN | "yes" |
| stagger_A_fields | CCTK_BOOLEAN | "yes" |
| initial_data_1D | KEYWORD | "Balsara1" |
| shock_direction | KEYWORD | "x" |
| discontinuity_position | CCTK_REAL | 0.0 |
| wave_amplitude | CCTK_REAL | 1.0e-3 |
| IsotropicGas_rho | CCTK_REAL | -1 |
| IsotropicGas_Y_e | CCTK_REAL | -1 |
| IsotropicGas_temperature | CCTK_REAL | -1 |
| ConstantDensitySphere_sphere_radius | CCTK_REAL | -1 |
| ConstantDensitySphere_rho_interior | CCTK_REAL | -1 |
| ConstantDensitySphere_Y_e_interior | CCTK_REAL | -1 |
| ConstantDensitySphere_T_interior | CCTK_REAL | -1 |
| ConstantDensitySphere_vx_interior | CCTK_REAL | 0 |
| ConstantDensitySphere_vy_interior | CCTK_REAL | 0 |
| ConstantDensitySphere_vz_interior | CCTK_REAL | 0 |
| ConstantDensitySphere_rho_exterior | CCTK_REAL | -1 |
| ConstantDensitySphere_Y_e_exterior | CCTK_REAL | -1 |
| ConstantDensitySphere_T_exterior | CCTK_REAL | -1 |
| impose_beta_equilibrium | CCTK_BOOLEAN | "no" |
| beq_temperature | CCTK_REAL | -1 |
| allow_native_entropy_proxy | CCTK_BOOLEAN | "no" |
| beq_residual_tolerance | CCTK_REAL | 1.0e-8 |

ParamCheck validates producer selections. Native wave amplitude is checked
finite and in [0,1) despite CCL's inclusive upper endpoint. Gas/sphere input
triples require positive finite density/temperature and effective bounds;
CCL nonnegative ranges do not replace these checks. Sphere radius must be
finite and nonnegative; its unrestricted velocity component declarations are
bounded by the body's finite, combined-speed check.

Beta checks the sentinel before temperature use; finite positive temperature
and residual tolerance are required. Native entropy-proxy opt-in defaults to
no. Magnetic staggering defaults to yes and now adds component-specific
half-cell offsets. GID-0001/0004/0007/0015 record resolved historical issues.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [param.ccl](../../../GRHayLID/param.ccl)

## Related Pages

- [Keyword Extensions](hydrobase-keyword-extensions.md)
- [GRHayLib Contract](grhaylib-contract.md)
- [Coverage Gaps](../validation/coverage-gaps.md)
- [GID-0001](../contradictions.md#gid-0001)
- [GID-0004](../contradictions.md#gid-0004)
- [GID-0007](../contradictions.md#gid-0007)
- [GID-0015](../contradictions.md#gid-0015)
