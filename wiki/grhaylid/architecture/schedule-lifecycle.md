# Declared Schedule Lifecycle

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

Eight schedule declarations include an unconditional parameter check at
CCTK_PARAMCHECK. Native hydro routines run in HydroBase_Initial and read
ADMBase::metric as well as Grid coordinates. Magnetic scheduling requires
both local selectors. Beta follows HydroBase_Initial; entropy follows both
that group and the beta alias, before HydroBase_Prim2ConInitial.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | ParamCheck validates the producer matrix. |
| HydroTest1D | HydroBase_Initial hydro writes. |
| HydroTest1D+Magnetic | Enabled both-selector magnetic writes after hydro. |
| IsotropicGas | Three-dimensional gas writes. |
| ConstantDensitySphere | Three-dimensional sphere writes. |
| BetaEquilibrium | Conditional Ye/T storage; five primitive writes including rho. |
| Entropy/Hybrid | Simple or Hybrid selection invokes the native-proxy routine. |
| Entropy/Tabulated | Tabulated entropy writes six fields. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `SCHEDULE-LIFECYCLE-01` | GRHayLID_ParamCheck has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_ParamCheck` |
| `SCHEDULE-LIFECYCLE-02` | GRHayLID_1D_tests_hydro_data has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_hydro_data` |
| `SCHEDULE-LIFECYCLE-03` | GRHayLID_1D_tests_magnetic_data has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_magnetic_data` |
| `SCHEDULE-LIFECYCLE-04` | GRHayLID_IsotropicGas has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_IsotropicGas` |
| `SCHEDULE-LIFECYCLE-05` | GRHayLID_ConstantDensitySphere has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_ConstantDensitySphere` |
| `SCHEDULE-LIFECYCLE-06` | GRHayLID_BetaEquilibrium has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |
| `SCHEDULE-LIFECYCLE-07` | GRHayLID_compute_entropy_hybrid has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` |
| `SCHEDULE-LIFECYCLE-08` | GRHayLID_compute_entropy_tabulated has a local schedule declaration. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` |

## Details

| Routine | Reads | Writes | Guard / order |
| --- | --- | --- | --- |
| ParamCheck | Parameter values in body | None | CCTK_PARAMCHECK |
| 1D hydro | Coordinates, metric | rho, press, eps, vel | HydroTest1D in HydroBase_Initial |
| 1D magnetic | Coordinates, metric | Avec, Bvec | HydroTest1D plus enabled magnetics and both GRHayLID selectors; after hydro |
| Gas / sphere | Coordinates, metric | rho, press, eps, vel, Ye, T | Corresponding hydro family in HydroBase_Initial |
| Beta | rho | rho, press, eps, Ye, T | impose_beta_equilibrium; after HydroBase_Initial; alias impose_beta_equilibrium |
| Native entropy | rho, press | entropy | initial_entropy=GRHayLID and EOS_type=Simple or Hybrid |
| Tabulated entropy | rho, Ye, T | rho, press, eps, entropy, Ye, T | initial_entropy=GRHayLID and EOS_type=Tabulated |

Beta and both entropy routines are declared at CCTK_INITIAL before
HydroBase_Prim2ConInitial. Entropy declares after
(HydroBase_Initial impose_beta_equilibrium). Beta requests one timelevel of
Ye/T storage independently of hydro family or optional selectors. Actual
storage and alias resolution remain framework boundaries. Unsupported entropy
EOS values are diagnosed by ParamCheck rather than an unselected schedule arm.

The corrected beta description omits entropy and includes atmosphere density
reset. Gas/sphere descriptions now name three-dimensional data. Historical
issues GID-0008, GID-0010, and GID-0011 are resolved by these local changes.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [schedule.ccl](../../../GRHayLID/schedule.ccl)

## Related Pages

- [Purpose and Build Surface](purpose-and-build-surface.md)
- [Keyword Extensions](../integration/hydrobase-keyword-extensions.md)
- [GID-0008](../contradictions.md#gid-0008)
- [GID-0010](../contradictions.md#gid-0010)
- [GID-0011](../contradictions.md#gid-0011)
