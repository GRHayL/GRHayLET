# One-Dimensional Magnetic Tests

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

Both magnetic selectors must name GRHayLID. The body checks storage for
Avec/Bvec/metric, checks the identity metric, rotates benchmark fields, and
writes Bvec at base coordinates. Avec uses component-specific transverse
half-cell shifts when stagger_A_fields=yes.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| HydroTest1D+Magnetic | Enabled both-selector magnetic production on Cartesian identity-metric data. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ONE-D-TESTS-MAGNETIC-01` | The body requires both selectors, checks storage/metric, and writes both destinations. | visible-implementation | Named local source | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` |
| `ONE-D-TESTS-MAGNETIC-02` | stagger_A_fields is consumed in component-specific coordinate expressions. | visible-implementation | Named local source | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` |
| `ONE-D-TESTS-MAGNETIC-03` | Both magnetic selectors extend HydroBase with GRHayLID. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` |
| `ONE-D-TESTS-MAGNETIC-04` | README declares normalized Balsara coefficients and the rescale_magnetics=no import requirement. | declared | Named local source | `doc:GRHayLID/README#section=1. Purpose` |
| `ONE-D-TESTS-MAGNETIC-05` | stagger_A_fields defaults to yes in CCL. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` |

## Details

### Benchmark states

| Arm | Left `(Bx, By, Bz)` | Right `(Bx, By, Bz)` |
| --- | --- | --- |
| `Balsara1` | `(0.5, 1.0, 0.0)` | `(0.5, -1.0, 0.0)` |
| `Balsara2` | `(5.0, 6.0, 6.0)` | `(5.0, 0.7, 0.7)` |
| `Balsara3` | `(10.0, 7.0, 7.0)` | `(10.0, 0.7, 0.7)` |
| `Balsara4` | `(10.0, 7.0, 7.0)` | `(10.0, -7.0, -7.0)` |
| `Balsara5` | `(2.0, 0.3, 0.3)` | `(2.0, -0.7, 0.5)` |
| `equilibrium`, `sound wave`, `shock tube` | `(0.0, 0.0, 0.0)` | `(0.0, 0.0, 0.0)` |

Direction rotation maps (Bx,By,Bz) to (Bz,Bx,By) for y and (By,Bz,Bx)
for z. Direct Bvec uses the base coordinate to select the side.

### Component placement

| Component | Coordinates when staggered |
| --- | --- |
| Ax | (x,y+dy/2,z+dz/2) |
| Ay | (x+dx/2,y,z+dz/2) |
| Az | (x+dx/2,y+dy/2,z) |

With staggering disabled all three components use base coordinates. Each
component chooses its side using its shifted longitudinal coordinate.
For x-oriented data the potential is (By*z-Bz*y,0,Bx*y); for y it is
(By*z,Bz*x-Bx*z,0); for z it is (0,Bz*x,Bx*y-By*x), evaluated separately
at each component's coordinates. Runtime curl, boundary stencil, and gauge
compatibility are not established by these assignments.

README requires Cartesian coordinates and an identity spatial metric from
ADMBase, and records the importer normalization requirement. Its external
consumer statements are documented prerequisites, not sibling-thorn evidence
in this branch. Historical GID-0002/0004/0012 are resolved by output documentation,
visible shifts, and matching both-selector guards.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [1D_tests_magnetic_data.c](../../../GRHayLID/src/1D_tests_magnetic_data.c)
- [param.ccl](../../../GRHayLID/param.ccl)
- [README](../../../GRHayLID/README)

## Related Pages

- [Hydro Tests](one-d-tests-hydro.md)
- [Parameters](../integration/parameters-and-configurations.md)
- [GID-0002](../contradictions.md#gid-0002)
- [GID-0004](../contradictions.md#gid-0004)
- [GID-0012](../contradictions.md#gid-0012)
