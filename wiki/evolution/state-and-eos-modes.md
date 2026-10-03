# State and EOS Modes

> Status: confirmed · Last reconciled: 10-02-2026
> Up: [Evolution](index.md)

## Summary

IllinoisGRMHD declares one five-component hydrodynamic conservative group.
For locally handled `EOS_type` values `Hybrid`, `Simple`, and `Tabulated`, the
primary `EOS_type`/`evolve_entropy` selector tree chooses one of four source
families. `Simple` follows the Hybrid branch, so `Simple` with entropy
evolution selects `HybridEntropy`. No local catch-all establishes a family for
another shared-selector value. An active `ID_converter_ILGRMHD` branch may
independently add non-entropy Hybrid scheduling; details belong to
[Schedule Lifecycle](../architecture/schedule-lifecycle.md) and
[Migration and Backward Compatibility](../integration/migration-and-backward-compatibility.md).

Claim evidence:

- Claim: For locally handled selectors, the primary schedule maps `Simple` plus `evolve_entropy=true` to HybridEntropy; no catch-all covers another shared-selector value, and the active compatibility branch can add Hybrid scheduling independently.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/schedule.ccl`, `EOS_type` and `evolve_entropy` family tree
- Corroboration: registered four family `make.code.defn` `SRCS` roles and matching scheduled function names
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=Hybrid,Simple,Tabulated with evolve_entropy false/true; date=07-17-2026`

## Detail

### State roles

- `grmhd_conservatives` evolves `rho_star`, `tau`, and the three `Stilde`
  components with three timelevels. `grmhd_conservatives_rhs` and
  `grmhd_flux_temps` hold their RHS and face-flux data.
- `grmhd_velocities` holds IllinoisGRMHD's primitive `v^i=u^i/u^0`; `u0`
  supports later consumers. Density, pressure, internal energy, entropy,
  electron fraction, and temperature are inherited HydroBase fields used as
  applicable to each family. Hybrid/Simple entropy uses a separate
  `hybrid_entropy` proxy; HydroBase entropy is marked unavailable with NaN.
- `ent_star` and `ent_star_rhs` receive scheduled storage only when
  `evolve_entropy` is true; `ent_star_flux` is stored unconditionally and is
  repeated in that conditional storage block. `Ye_star` and `Ye_star_rhs`
  receive scheduled storage only when `EOS_type` is `Tabulated`;
  `Ye_star_flux` is stored unconditionally and is repeated in that
  conditional storage block.
- Reconstructed-state temporaries hold left/right velocities and staggered B
  data shared by hydro and induction work. Exact A/B placement belongs to
  [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md).

### Primary handoff for locally handled selectors

| `EOS_type` selector | `evolve_entropy` | Scheduled/build family | Carried extras |
| --- | --- | --- | --- |
| `Hybrid` or `Simple` | false | `src/Hybrid/`; `IllinoisGRMHD_hybrid_*` | no entropy or electron-fraction conservative |
| `Hybrid` or `Simple` | true | `src/HybridEntropy/`; `IllinoisGRMHD_hybrid_entropy_*` | thorn-owned `hybrid_entropy`; `ent_star`, entropy RHS, and entropy flux |
| `Tabulated` | false | `src/Tabulated/`; `IllinoisGRMHD_tabulated_*` | primitive `Y_e` and `temperature`; `Ye_star`, electron-fraction RHS, and flux |
| `Tabulated` | true | `src/TabulatedEntropy/`; `IllinoisGRMHD_tabulated_entropy_*` | compiled body only; startup rejects this entropy mode |

Each branch schedules its own Prim2Con, Con2Prim, matter-boundary, source,
flux, initial-perturbation, and conservative-perturbation function. All four
directories are also independently listed in the common build. These are
declared handoffs, not evidence that any branch was executed by a fixture.
The primary tree has no catch-all after its `Hybrid`/`Simple` and `Tabulated`
conditions. Separately, an active `ID_converter_ILGRMHD` condition schedules
the non-entropy Hybrid family without an `EOS_type` or `evolve_entropy` guard.
The shared selector defaults and behavior of any other selector value live
outside this thorn and are not inferred.

### Per-family dataflow differences

- Every Prim2Con function assembles rho, pressure, velocity, and centered B.
  Entropy families additionally assemble/write entropy; tabulated families
  additionally assemble/write `Y_e` and temperature. Their computed
  conservatives add `ent_star`, `Ye_star`, or both.
- Source initialization sets `rho_star_rhs` to zero. It also zeros each
  enabled advected-extra RHS (`ent_star_rhs`, `Ye_star_rhs`) before flux
  contributions. Primitive source input carries the same family extras.
- PPM reconstructs rho, pressure, transverse B, and velocities in every
  family. HybridEntropy also reconstructs entropy; Tabulated reconstructs
  `Y_e`; TabulatedEntropy reconstructs both. Tabulated code seeds face
  temperature from the current cell and locally calls the tabulated
  pressure-to-thermodynamics routine after bounds enforcement. Hybrid code
  supplies a locally calculated effective Gamma to density steepening;
  tabulated code supplies `1.0` per its source comment.
- Family-specific GRHayL characteristic-speed, HLLE-flux, limiting, EOS, and
  inversion calls are treated only as call boundaries here. Their external
  implementations are not inferred.

### Supported entropy boundary and Simple thermodynamics

At `CCTK_WRAGH`, `IllinoisGRMHD_check_eos_support` follows modern or legacy
GRHayL initialization. It rejects tabulated entropy evolution until the shared
packing/flux/recovery/atmosphere contract is corrected. Hybrid entropy requires
`neos=1` and `Gamma_th=Gamma_ppoly[0]`. Simple entropy uses the same internal
proxy machinery. It also rejects `Palenzuela1D` and `Palenzuela1D_entropy`
as the main or any backup recovery routine when `EOS_type` is `Simple`.
No claim is made that a rejected mode evolved successfully.

For Simple EOS, the installed callback computes epsilon, h, and sound speed
directly from pressure/density and Gamma. Other EOS calls delegate to the
original callback. This avoids auxiliary cold-curve subtraction in the enthalpy and
sound-speed evaluation only; it does not modify GRHayL recovery internals.

Claim evidence:
- Claim: Startup rejects unsupported entropy selections and both Palenzuela variants in Simple main/backup recovery, and installs direct Simple thermodynamics while preserving other EOS callbacks.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/check_eos_support.c, IllinoisGRMHD_check_eos_support and IllinoisGRMHD_compute_h_and_cs2`
- Corroboration: registered `IllinoisGRMHD/schedule.ccl, affected declarations`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`

## Sources

- `IllinoisGRMHD/src/check_eos_support.c` — startup guard and EOS callback.

- `IllinoisGRMHD/interface.ccl` — groups `grmhd_conservatives`, `ent_star`,
  `Ye_star`, reconstructed temporaries, RHS groups, and flux groups.
- `IllinoisGRMHD/schedule.ccl` — conditional `STORAGE` blocks and the
  `EOS_type`/`evolve_entropy` family schedule tree.
- `IllinoisGRMHD/src/make.code.defn`,
  `IllinoisGRMHD/src/Hybrid/make.code.defn`,
  `IllinoisGRMHD/src/HybridEntropy/make.code.defn`,
  `IllinoisGRMHD/src/Tabulated/make.code.defn`, and
  `IllinoisGRMHD/src/TabulatedEntropy/make.code.defn` — common `SUBDIRS` and
  independent family `SRCS` roles.
- Registered aggregates `illinoisgrmhd-hybrid`,
  `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and
  `illinoisgrmhd-tabulated-entropy` — four `*_prims_to_conservs`,
  `*_evaluate_sources_rhs`, and `*_calculate_flux_dir_rhs` function families.

## See Also

- Parent: [Evolution](index.md)
- Depends on: [Cactus Surface and Build](../architecture/cactus-surface-and-build.md)
- See also: [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
- See also: [Reconstruction, Fluxes, and Sources](reconstruction-fluxes-and-sources.md)
