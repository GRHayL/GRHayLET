# ConstantDensitySphere Initial Data

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page records ConstantDensitySphere precondition text, ten local parameter
declarations, seven visible sentinel checks, two tabulated-EOS calls, radius
selection, and pointwise HydroBase assignments. It does not establish
error-macro termination, EOS-table behavior, returned values, successful
scheduling, geometry semantics, or physical and numerical validity.

## Summary

`GRHayLID_ConstantDensitySphere` visibly checks for Tabulated `EOS_type` and
`GRHayLID` selections for electron fraction and temperature. Seven
non-velocity controls pass through `CHECK_PARAMETER`; three interior velocity
controls do not. Separate EOS calls populate interior and exterior pressure
and internal-energy variables. Points with `r` greater than the declared
radius receive exterior state and zero velocity; all other points receive
interior state and its three declared velocity components.

## Mode Applicability

| Applicability | Visible initial-data behavior |
| --- | --- |
| ConstantDensitySphere | Two delegated thermodynamic states selected by radius, with configured interior velocity and zero exterior velocity. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-SPHERE-01` | Function visibly performs three precondition comparisons, seven sentinel checks, two tabulated-EOS calls, and radius-selected pointwise writes. | visible-implementation | Function body | `c:GRHayLID/src/ConstantDensitySphere.c#symbol=GRHayLID_ConstantDensitySphere` |
| `ID-SPHERE-02` | Sphere radius is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_sphere_radius` |
| `ID-SPHERE-03` | Interior density is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_interior` |
| `ID-SPHERE-04` | Interior electron fraction is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_interior` |
| `ID-SPHERE-05` | Interior temperature is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_interior` |
| `ID-SPHERE-06` | Interior x-velocity is unrestricted and defaults to zero. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vx_interior` |
| `ID-SPHERE-07` | Interior y-velocity is unrestricted and defaults to zero. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vy_interior` |
| `ID-SPHERE-08` | Interior z-velocity is unrestricted and defaults to zero. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vz_interior` |
| `ID-SPHERE-09` | Exterior density is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_exterior` |
| `ID-SPHERE-10` | Exterior electron fraction is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_exterior` |
| `ID-SPHERE-11` | Exterior temperature is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_exterior` |
| `ID-SPHERE-12` | ThornGuide prints the EOS selection as lowercase `"tabulated"` and the electron-fraction and temperature selections as `"GRHayLID"`. | declared | ConstantDensitySphere subsection | `doc:GRHayLID/doc/documentation.tex#section=ConstantDensitySphere` |

## Details

### Preconditions and parameter checks

The opening comparisons call `CCTK_ERROR` unless `EOS_type` equals
`Tabulated`, `initial_Y_e` equals `GRHayLID`, and `initial_temperature`
equals `GRHayLID`. ThornGuide prints the EOS selection as lowercase
`"tabulated"`. Equivalence of those EOS spellings is external; see
[GID-0014](../contradictions.md#gid-0014).

`CHECK_PARAMETER` is visibly applied to sphere radius and interior/exterior
density, electron fraction, and temperature: seven controls total. All seven
declare nonnegative values, a forbidden `-1`, and default `-1`. The three
interior velocities permit any real, default to zero, and have no visible
`CHECK_PARAMETER` call.

### Two delegated thermodynamic states

Before its point loop, the function calls
`ghl_tabulated_compute_P_eps_from_T` twice:

- interior density, electron fraction, and temperature populate
  `P_interior` and `eps_interior`;
- exterior density, electron fraction, and temperature populate
  `P_exterior` and `eps_exterior`.

No local return code is captured. GRHayLib table semantics, bounds behavior,
and result validity remain out of scope.

### Radius branch and visible writes

Each point first receives zero in all three velocity components. When
`r[index] > ConstantDensitySphere_sphere_radius`, the exterior branch writes
exterior `rho`, `Y_e`, `temperature`, `press`, and `eps`, then repeats all
three zero-velocity assignments. Otherwise, including equality, the interior
branch writes interior thermodynamics and overwrites velocity with
`ConstantDensitySphere_vx_interior`, `_vy_interior`, and `_vz_interior`.

The duplicate exterior zero writes and the equality-side choice are visible
source facts. Coordinate `r` meaning, grid geometry, boundary placement, and
resulting sphere properties are external or runtime concerns.

## Caveats

- Visible comparisons and macro calls do not establish termination.
- Two external EOS calls are visible; no checked-in local oracle validates
  their outputs.
- Schedule description calls this three-dimensional setup a "1D test"; see
  [GID-0011](../contradictions.md#gid-0011).
- ThornGuide and the function use different EOS literal capitalization; see
  [GID-0014](../contradictions.md#gid-0014).
- Radius comparison and assignments do not prove geometric or numerical
  correctness.

## Sources

- [ConstantDensitySphere implementation](../../../GRHayLID/src/ConstantDensitySphere.c)
- [Parameter declarations](../../../GRHayLID/param.ccl)
- [ThornGuide source](../../../GRHayLID/doc/documentation.tex)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)

## Related Pages

- [IsotropicGas Initial Data](isotropic-gas.md)
- [GRHayLib Contract](../integration/grhaylib-contract.md)
- [GID-0011](../contradictions.md#gid-0011)
- [GID-0014](../contradictions.md#gid-0014)
