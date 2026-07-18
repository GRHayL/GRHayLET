# Migration and Backward Compatibility

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Integration](index.md)

## Summary

ThornGuide documents migration from old IllinoisGRMHD parfiles to GRHayLib-
based controls. Current tree still compiles compatibility C files, retains
deprecated variables and parameters, and conditionally schedules compatibility
when `ID_converter_ILGRMHD` is active. This current-tree fact does not guarantee
support in any Einstein Toolkit release.

Claim status: stale; contradiction: CONTR-0001.
[Contradiction record](../contradictions.md#contr-0001)
documents ThornGuide's `ET_2024_11` sunset statement versus current compiled,
conditionally scheduled compatibility surfaces.

Claim evidence:

- Claim: Current tree retains compiled, conditionally scheduled compatibility,
  but this does not guarantee support in any external Einstein Toolkit release.
- Role: public/scientific contract
- Deciding authority: `IllinoisGRMHD/src/make.code.defn::SRCS` and
  `IllinoisGRMHD/schedule.ccl::Backward compatibility scheduling`
- Corroboration: `IllinoisGRMHD/doc/documentation.tex::Updating Old Parfiles`
  supplies conflicting sunset claim; `IllinoisGRMHD/interface.ccl::backward compatibility`
  retains old variables
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=static source/build/CCL/doc inspection; date=07-17-2026`

## Detail

### Old-to-new parfile map

Mappings below are ThornGuide migration guidance, not independently inferred
equivalence:

| Old setting | Documented replacement |
| --- | --- |
| `ID_converter_ILGRMHD::Gamma_Initial` | `GRHayLib::Gamma_ppoly_in[0]` |
| `ID_converter_ILGRMHD::random_seed` | `IllinoisGRMHD::random_seed` |
| `ID_converter_ILGRMHD::random_pert` | `IllinoisGRMHD::random_pert` |
| `ID_converter_ILGRMHD::K_Initial` | `GRHayLib::k_ppoly0` |
| `Convert_to_HydroBase::Convert_to_HydroBase_every` | `IllinoisGRMHD::Convert_to_HydroBase_every` |
| `IllinoisGRMHD::GAMMA_SPEED_LIMIT` | `GRHayLib::max_Lorentz_factor` |
| `IllinoisGRMHD::K_poly` | `GRHayLib::k_ppoly0` |
| `IllinoisGRMHD::rho_b_atm` | `GRHayLib::rho_b_atm` |
| `IllinoisGRMHD::rho_b_max` | `GRHayLib::rho_b_max` |
| `IllinoisGRMHD::Psi6threshold` | `GRHayLib::Psi6threshold` |
| `IllinoisGRMHD::neos` | `GRHayLib::neos` |
| `IllinoisGRMHD::gamma_th` | `GRHayLib::Gamma_th` |
| `IllinoisGRMHD::damp_lorenz` | `GRHayLib::Lorenz_damping_factor` |

Guide says remove thorns `ID_converter_ILGRMHD` and `Convert_to_HydroBase` and
their parameters from updated parfiles. It says
`ID_converter_ILGRMHD::pure_hydro_run` has no direct replacement because
GRHayLHD provides that role. For migrated random perturbation, it additionally
requires `IllinoisGRMHD::perturb_initial_data` to trigger operation.

Guide marks `IllinoisGRMHD::tau_atm` and
`IllinoisGRMHD::conserv_to_prims_debug` for removal: first is described as
automatically computed by GRHayLib; second feature as unavailable. It also
marks `verbose` options `essential` and `essential+iteration output` deprecated.
Current `param.ccl` still accepts all of these.

### Magnetic-definition migration

ThornGuide says GRHayL-based quantities use old magnetic quantities rescaled
by `(4*pi)^(-1/2)`. To support old initial-data thorns, it says default assumes
HydroBase B and A use old definition. Current `rescale_magnetics=yes` implements
that ingress factor and inverse egress factor. Setting `no` makes both factors
one. These are local conversion facts and documented migration intent; no
claim is made about definitions inside external thorns.

### Compatibility retained in current tree

`src/make.code.defn` includes both `backward_compatible_initialize.c` and
`backward_compatible_data.c` in current source list. `interface.ccl` retains
deprecated `psi6phi`, `rho_b`, `P`, and `Bx/By/Bz` grid functions. `param.ccl`
retains ten old controls listed on parameter owner page.

All compatibility scheduling is under
`CCTK_IsThornActive("ID_converter_ILGRMHD")`:

- storage for deprecated groups is enabled;
- `IllinoisGRMHD_backward_compatible_initialize` is declared at `CCTK_WRAGH`
  with `GLOBAL` option;
- `IllinoisGRMHD_backward_compatible_data` is declared at `CCTK_PostInitial`
  and in `MoL_PostRHS`;
- HydroBase egress is declared at initial conversion and `CCTK_ANALYSIS`;
- non-entropy Hybrid evolution functions are declared in variant groups.

Initializer allocates GRHayL parameter/EOS structures, rejects `neos>1`, and
sets a fixed compatibility configuration including Noble2D primary and Font1D
first backup, no entropy or temperature evolution, and hybrid EOS setup from
deprecated parameters. These are visible call inputs, not claims about GRHayL
algorithms.

Data copier assigns HydroBase `rho/press`, centered B, and `phitilde` into
deprecated `rho_b/P`, `Bx/By/Bz`, and `psi6phi`. Conversion routine separately
checks whether old `Convert_to_HydroBase` thorn is active and, if so, reads its
cadence dynamically.

Neither `configuration.ccl` nor current source list establishes availability
of old thorns. Schedule gate only describes behavior if active.

### Sunset statement versus release support

ThornGuide says backward compatibility “ends in the `ET_2024_11` release.”
Current code/build/CCL decide what this checkout contains and can declare for
scheduling; they do not prove what any released toolkit includes, supports,
or will remove. `CONTR-0001` remains stale documentation conflict until a
maintainer/release decision and source/documentation reconciliation resolves
it.

## Sources

- [`IllinoisGRMHD/doc/documentation.tex`](../../IllinoisGRMHD/doc/documentation.tex) —
  `Updating Old Parfiles` mappings, magnetic migration, and sunset statement.
- [`IllinoisGRMHD/interface.ccl`](../../IllinoisGRMHD/interface.ccl) —
  `backward compatibility` deprecated variable groups.
- [`IllinoisGRMHD/param.ccl`](../../IllinoisGRMHD/param.ccl) — deprecated
  parameter block and accepted verbose options.
- [`IllinoisGRMHD/schedule.ccl`](../../IllinoisGRMHD/schedule.ccl) —
  `Backward compatibility scheduling` active-thorn gate and schedule sites.
- [`IllinoisGRMHD/src/make.code.defn`](../../IllinoisGRMHD/src/make.code.defn) —
  current `SRCS` inclusion.
- [`IllinoisGRMHD/src/backward_compatible_initialize.c`](../../IllinoisGRMHD/src/backward_compatible_initialize.c) —
  `IllinoisGRMHD_backward_compatible_initialize`.
- [`IllinoisGRMHD/src/backward_compatible_data.c`](../../IllinoisGRMHD/src/backward_compatible_data.c) —
  `IllinoisGRMHD_backward_compatible_data`.
- [`IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c`](../../IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c) —
  retained old-thorn cadence lookup.

## See Also

- Parent: [Integration](index.md)
- Contrasts with: [Parameters and Runtime Controls](parameters-and-runtime-controls.md)
- Depends on: [Cactus Surface and Build](../architecture/cactus-surface-and-build.md)
- See also: [HydroBase, GRHayLib, and Tmunu](hydrobase-grhaylib-and-tmunu.md)
