# Parameters and Runtime Controls

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Integration](index.md)

## Summary

`param.ccl` declares 22 IllinoisGRMHD-owned parameters and consumes two shared
GRHayLib selectors. Tables below preserve declared type, legal range or
keywords, default, and explicit steerability. “Not declared” means this file
has no `STEERABLE` clause; no external default is invented.

Claim evidence:

- Claim: Local `param.ccl` owns 22 parameter declarations and consumes two
  shared selectors; tables preserve every locally visible type, range or
  keyword, default, and explicit steerability clause.
- Role: public/scientific contract
- Deciding authority: `IllinoisGRMHD/param.ccl` complete declaration surface
- Corroboration: `IllinoisGRMHD/schedule.ccl` selector and control use sites
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=static row-by-row inspection; date=07-17-2026`

## Detail

### Current IllinoisGRMHD controls

For `CCTK_BOOLEAN` rows, `param.ccl` has an empty explicit range block; table
uses “boolean” for standard CCL boolean values rather than inventing entries.

| Parameter | Type | Legal values/range | Default | Steerability | Local role |
| --- | --- | --- | --- | --- | --- |
| `Convert_to_HydroBase_every` | `INT` | `0:*` | `0` | `RECOVER` | Zero disables locally cadence-gated schedule sites; positive N converts on iterations divisible by N. |
| `update_Tmunu` | `CCTK_BOOLEAN` | boolean | `yes` | `RECOVER` | Gates `AddToTmunu` routine and MoL constrained-group registration. |
| `rescale_magnetics` | `CCTK_BOOLEAN` | boolean | `yes` | not declared | Selects inverse/same `(4*pi)^(1/2)` factors at HydroBase ingress/egress. |
| `Symmetry` | `KEYWORD` | `none` | `none` | not declared | Only public keyword is `none`; description says equatorial support is in progress. |
| `Sym_Bz` | `REAL` | `-1.0:1.0` | `1.0` | not declared | Supplies z parity factors in magnetic/symmetry code. Description asks for `+1` or `-1`, while syntactic range spans interval. |
| `Matter_BC` | `KEYWORD` | `copy`, `outflow`, `frozen` | `outflow` | not declared | Selects matter boundary behavior. |
| `EM_BC` | `KEYWORD` | `copy`, `frozen` | `copy` | not declared | Selects electromagnetic boundary behavior. |
| `verbose` | `KEYWORD` | `no`, `yes`, `essential`, `essential+iteration output` | `yes` | `ALWAYS` | Current recovery files test only `yes`; latter two keywords are declared deprecated. |
| `random_seed` | `INT` | `0:99999999` | `0` | `ALWAYS` | Seeds `srand()` in perturbation routines. |
| `random_pert` | `REAL` | `*:*` | `0` | `ALWAYS` | Multiplicative perturbation magnitude. |
| `perturb_initial_data` | `CCTK_BOOLEAN` | boolean | `no` | not declared | Gates primitive perturbation after HydroBase ingress and before Prim2Con. |
| `perturb_every_con2prim` | `CCTK_BOOLEAN` | boolean | `no` | `ALWAYS` | Gates conservative perturbation before every scheduled Con2Prim. |

Frozen matter and EM modes must be selected together; `IllinoisGRMHD_InitSymBound`
errors if exactly one is `frozen`. Detailed algorithms belong to
[Matter Boundaries and Perturbations](../evolution/matter-boundaries-and-perturbations.md)
and [Electromagnetic Boundaries and Symmetry](../magnetics/electromagnetic-boundaries-and-symmetry.md).

`Convert_to_HydroBase_every=0` removes two locally guarded conversion schedule
sites. Active NRPyLeakageET site has different gate, while conversion routine
still performs cadence modulo; see
[HydroBase, GRHayLib, and Tmunu](hydrobase-grhaylib-and-tmunu.md).

### Deprecated IllinoisGRMHD controls

These remain declared. Most are consumed only by retained
`ID_converter_ILGRMHD` compatibility initializer; `tau_atm` and
`conserv_to_prims_debug` have no non-declaration use in current tree.

| Parameter | Type | Legal values/range | Default | Steerability | Current local use |
| --- | --- | --- | --- | --- | --- |
| `GAMMA_SPEED_LIMIT` | `REAL` | `1:*` (description says positive `>1`) | `10.0` | not declared | Compatibility initializer argument. |
| `tau_atm` | `REAL` | `0:*` | `1e100` | `ALWAYS` | No current C use; ThornGuide says replacement is automatic in GRHayLib. |
| `rho_b_atm` | `REAL` | `*:*` | `1e200` | `ALWAYS` | Compatibility hybrid-EOS initialization. |
| `rho_b_max` | `REAL` | `0:*` | `1e300` | `ALWAYS` | Compatibility hybrid-EOS initialization. |
| `conserv_to_prims_debug` | `INT` | `0:1` | `0` | `ALWAYS` | No current C use; ThornGuide says feature is unavailable. |
| `Psi6threshold` | `REAL` | `*:*` | `1e100` | `ALWAYS` | Compatibility parameter initialization. |
| `neos` | `INT` | `1:10` | `1` | not declared | Compatibility initializer errors for values above one, then supplies hybrid-EOS initialization. |
| `gamma_th` | `REAL` | `0:*`; discrete `-1` sentinel described as forbidden | `-1` | not declared | Compatibility polytropic and thermal gamma input; sentinel default is intended to force explicit old-parfile setting. |
| `K_poly` | `REAL` | `0:*` | `1.0` | not declared | Compatibility hybrid-EOS initialization. |
| `damp_lorenz` | `REAL` | `*:*` | `0.0` | `ALWAYS` | Compatibility parameter initialization. |

Migration replacements and retained compatibility activation belong to
[Migration and Backward Compatibility](migration-and-backward-compatibility.md).

### Shared selectors consumed here

| Shared parameter | Local declaration | Local behavior | Unknown at this boundary |
| --- | --- | --- | --- |
| `GRHayLib::EOS_type` | `USES KEYWORD EOS_type` | `Simple` or `Hybrid` selects hybrid family; `Tabulated` selects tabulated family and enables `Ye_star` storage. | Owner's legal keyword set, default, and steerability are not declared here. |
| `GRHayLib::evolve_entropy` | `USES CCTK_BOOLEAN evolve_entropy` | True enables entropy storage and entropy variant within selected EOS family. | Owner's default and steerability are not declared here. |

IllinoisGRMHD `schedule.ccl` contains no final catch-all family for an
`EOS_type` outside its tested strings. This describes local schedule selection,
not GRHayLib's public parameter contract. Case files explicitly set only
`Simple` or `Hybrid`; none explicitly sets `evolve_entropy`.

### Interaction ownership

- `update_Tmunu`, rescaling, velocity conversion, and cadence:
  [HydroBase, GRHayLib, and Tmunu](hydrobase-grhaylib-and-tmunu.md).
- EOS/entropy family schedule selection:
  [State and EOS Modes](../evolution/state-and-eos-modes.md).
- Matter boundary and perturbation behavior:
  [Matter Boundaries and Perturbations](../evolution/matter-boundaries-and-perturbations.md).
- EM boundary and symmetry behavior:
  [Electromagnetic Boundaries and Symmetry](../magnetics/electromagnetic-boundaries-and-symmetry.md).

## Sources

- [`IllinoisGRMHD/param.ccl`](../../IllinoisGRMHD/param.ccl) — complete
  Illinois-owned declaration surface and two `USES` declarations; tables were
  checked row-by-row against this file.
- [`IllinoisGRMHD/schedule.ccl`](../../IllinoisGRMHD/schedule.ccl) — cadence,
  Tmunu, perturbation, EOS, entropy, and compatibility conditions.
- [`IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c`](../../IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c) and
  [`convert_IllinoisGRMHD_to_HydroBase.c`](../../IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c) — rescaling and cadence use sites.
- [`IllinoisGRMHD/src/InitSymBound.c`](../../IllinoisGRMHD/src/InitSymBound.c) —
  boundary coupling and symmetry controls.
- [`IllinoisGRMHD/src/MoL_registration.c`](../../IllinoisGRMHD/src/MoL_registration.c) —
  `update_Tmunu` and shared-selector effects.
- [`IllinoisGRMHD/src/backward_compatible_initialize.c`](../../IllinoisGRMHD/src/backward_compatible_initialize.c) —
  deprecated parameter consumers.
- [`IllinoisGRMHD/src/Hybrid/perturb_primitives.c`](../../IllinoisGRMHD/src/Hybrid/perturb_primitives.c) and
  [`perturb_conservatives.c`](../../IllinoisGRMHD/src/Hybrid/perturb_conservatives.c) — representative seed/magnitude use; parallel files exist in all four families.
- [`IllinoisGRMHD/doc/documentation.tex`](../../IllinoisGRMHD/doc/documentation.tex) —
  `Parameters` and `Updating Old Parfiles` intent.

## See Also

- Parent: [Integration](index.md)
- Contrasts with: [Migration and Backward Compatibility](migration-and-backward-compatibility.md)
- Example: [Balsara and TOV Cases](../validation/balsara-and-tov-cases.md)
- Depends on: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
