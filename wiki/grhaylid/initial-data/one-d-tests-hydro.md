# One-Dimensional Hydrodynamic Tests

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

All eight declared selectors have live pre-loop dispatch arms. The sound
wave sets a direction-selected sinusoidal velocity on rho=P=1. Simple energy
uses Gamma_ppoly[0]; Hybrid calls the cold-plus-thermal energy helper and
rejects pressures below its computed cold pressure. Native tests check an
identity spatial metric supplied by ADMBase.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| HydroTest1D | Simple/Hybrid hydro setups with live sound-wave dispatch and direction rotation. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ONE-D-TESTS-HYDRO-01` | The hydro body initializes all dispatch states, rotates step velocities, and selects the sound-wave velocity component. | visible-implementation | Named local source | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` |
| `ONE-D-TESTS-HYDRO-02` | The local selector admits eight choices, default Balsara1. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `ONE-D-TESTS-HYDRO-03` | Local code requires finite wave amplitude in [0,1). | visible-implementation | Named local source | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` |
| `ONE-D-TESTS-HYDRO-04` | Hybrid energy helper semantics remain external. | out-of-scope | Named local source | `c:GRHayLID/src/1D_tests_hydro_data.c#call=ghl_hybrid_compute_epsilon?function=GRHayLID_1D_tests_hydro_data` |
| `ONE-D-TESTS-HYDRO-05` | wave_amplitude defaults to 1e-3 in CCL. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` |

## Details

### Step-state inventory

| Arm | Left `(rho, press; vx, vy, vz)` | Right `(rho, press; vx, vy, vz)` |
| --- | --- | --- |
| `Balsara1` | `(1.0, 1.0; 0, 0, 0)` | `(0.125, 0.1; 0, 0, 0)` |
| `Balsara2` | `(1.0, 30.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |
| `Balsara3` | `(1.0, 1000.0; 0, 0, 0)` | `(1.0, 0.1; 0, 0, 0)` |
| `Balsara4` | `(1.0, 0.1; 0.999, 0, 0)` | `(1.0, 0.1; -0.999, 0, 0)` |
| `Balsara5` | `(1.08, 0.95; 0.4, 0.3, 0.2)` | `(1.0, 1.0; -0.45, -0.2, 0.2)` |
| `equilibrium` | `(1.0, 1.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |
| `shock tube` | `(2.0, 2.0; 0, 0, 0)` | `(1.0, 1.0; 0, 0, 0)` |

For y, velocity rotation maps (vx,vy,vz) to (vz,vx,vy); for z it maps
to (vy,vz,vx). Left state includes equality at discontinuity_position.
Sound-wave pre-loop values are fully initialized before rotation. The loop
sets rho=P=1 and all velocities to zero, then fills the selected component
with wave_amplitude*sin(pi*step). README describes this as a longitudinal
velocity perturbation rather than a pure traveling eigenmode. Thermal
pressure excludes kinetic energy; the former incompleteness comment is removed.

For Hybrid, a cold-pressure/energy lookup precedes the epsilon helper.
Pressure below the returned cold pressure is diagnosed. Simple retains
P/[rho*(Gamma_ppoly[0]-1)]. Nonfinite epsilon is diagnosed after either branch.
External helper semantics and initialized EOS metadata remain unverified here.
Historical dispatch, pressure-comment, unsupported-capability, and EOS-name
issues GID-0003/0005/0006/0013 have local resolution evidence.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [1D_tests_hydro_data.c](../../../GRHayLID/src/1D_tests_hydro_data.c)
- [param.ccl](../../../GRHayLID/param.ccl)

## Related Pages

- [Magnetic Tests](one-d-tests-magnetic.md)
- [GRHayLib Contract](../integration/grhaylib-contract.md)
- [GID-0003](../contradictions.md#gid-0003)
- [GID-0005](../contradictions.md#gid-0005)
- [GID-0006](../contradictions.md#gid-0006)
- [GID-0013](../contradictions.md#gid-0013)
