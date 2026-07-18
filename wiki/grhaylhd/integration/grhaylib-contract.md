# GRHayLib Contract

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page inventories GRHayLib handles and API-call families visibly used by
GRHayLHD. Calls prove local delegation only. GRHayLib initialization,
algorithms, parameter semantics, error contracts, and numerical properties
remain externally unverified.

## Summary

Local routines read `ghl_params` and `ghl_eos`, then delegate metric setup,
primitive limits, EOS work, conversion/recovery, reconstruction, wave speeds,
HLLE fluxes, sources, atmosphere handling, diagnostics, and stress-energy
calculation to GRHayLib APIs. Local CCL dispatches on shared `EOS_type` and
`evolve_entropy`; equivalence between those CCL values and library objects is
not established locally.

## Variant Applicability

| Applicability | Visible library boundary |
| --- | --- |
| Common | Metric/auxiliary setup, primitive limits, conservative conversion, multi-method recovery, PPM, characteristic speeds, sources, error handling, and stress-energy helpers. |
| Hybrid/Simple | Hybrid HLLE calls, hybrid EOS helper, and explicit `ghl_hybrid_Font1D` recovery fallback. |
| Hybrid/Simple+Entropy | Hybrid-entropy HLLE calls and entropy-bearing primitive/conservative paths. |
| Tabulated | Tabulated bounds/EOS calls, tabulated HLLE calls, and electron-fraction paths. |
| Tabulated+Entropy | Tabulated-entropy HLLE calls plus both electron-fraction and entropy paths. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-GHL-01` | Local registration reads `ghl_params->evolve_entropy` and `ghl_eos->eos_type`. | visible-implementation | MoL registration branches | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` |
| `INT-GHL-02` | Hybrid recovery visibly delegates primary recovery to `ghl_con2prim_multi_method` and has an explicit Font1D fallback. | visible-implementation | Hybrid recovery function | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` |
| `INT-GHL-03` | Tabulated-entropy flux code visibly delegates reconstruction, EOS bounds, directional speeds, and HLLE fluxes. | visible-implementation | TabulatedEntropy flux function | `c:GRHayLHD/src/TabulatedEntropy/evaluate_fluxes_rhs.c#symbol=GRHayLHD_tabulated_entropy_evaluate_fluxes_rhs` |
| `INT-GHL-04` | Stress-energy path visibly calls metric, auxiliary, and Tmunu helpers. | visible-implementation | Tmunu function | `c:GRHayLHD/src/compute_Tmunu.c#symbol=GRHayLHD_compute_Tmunu` |
| `INT-GHL-05` | README states support including all Con2Prim routines. | declared | Purpose section | `doc:GRHayLHD/README#section=1. Purpose` |
| `INT-GHL-06` | Internal behavior of called GRHayLib APIs is delegated and unverified locally. | out-of-scope | Local calls expose no library implementation | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` |
| `INT-GHL-07` | Parameter CCL uses shared `EOS_type`. | declared | Shared parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=EOS_type` |
| `INT-GHL-08` | Parameter CCL uses shared `evolve_entropy`. | declared | Shared parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=evolve_entropy` |
| `INT-GHL-09` | Local header visibly includes `GRHayLib.h`. | visible-implementation | Include directive | `macro:GRHayLHD/src/GRHayLHD.h#include=GRHayLib.h` |

## Details

### Shared handles and selection boundary

`GRHayLHD.h` includes `GRHayLib.h`. Local C uses global `ghl_params` for
entropy selection, primitive/reconstruction limits, atmosphere thresholds,
and recovery configuration; it uses `ghl_eos` for EOS-family selection and
EOS-dependent calls. Parameter CCL only `USES` shared `EOS_type` and
`evolve_entropy`. Schedule CCL selects four local source families from those
values, while registration checks library-object fields. Their initialization
and agreement are external boundaries.

### Visible API-call families

- Metric and geometry: `ghl_initialize_metric`,
  `ghl_enforce_detgtij_and_initialize_ADM_metric`,
  `ghl_compute_ADM_auxiliaries`, and
  `ghl_initialize_extrinsic_curvature`.
- Limits and errors: `ghl_limit_v_and_compute_u0`,
  `ghl_enforce_primitive_limits_and_compute_u0`,
  `ghl_apply_conservative_limits`, and `ghl_abort_if_error`.
- EOS work: `ghl_hybrid_get_K_and_Gamma`,
  `ghl_tabulated_enforce_bounds_rho_Ye_P`, and
  `ghl_tabulated_compute_eps_T_from_P`.
- Conversion and recovery: `ghl_compute_conservs`,
  `ghl_undensitize_conservatives`, `ghl_con2prim_multi_method`,
  `ghl_hybrid_Font1D`, `ghl_set_prims_to_constant_atm`, and diagnostic
  initialization.
- Reconstruction and fluxes: `ghl_compute_ftilde`, both PPM entry points,
  three directional characteristic-speed entry points, and directional HLLE
  entry points specialized for each source family.
- Sources and stress energy: `ghl_calculate_source_terms` and
  `ghl_compute_TDNmunu`.

These names and local call order are visible implementation. This KB does not
reconstruct called algorithms from sibling sources or infer accuracy,
convergence, supported parameter sets, or failure semantics.

### Documentation boundary

README says thorn supports most library features, including all Con2Prim
routines. Local recovery code proves calls to multi-method delegation and an
explicit Hybrid-family Font1D fallback. It does not enumerate or verify every
runtime-selectable external Con2Prim implementation, so broader statement
remains documentation intent.

## Caveats

- `EOS_type` and `evolve_entropy` are shared selections; local CCL declares no
  defaults for them.
- Type and function names do not establish ABI compatibility, successful
  linking, or runtime behavior.
- Absence of an explicit Tabulated Font1D call is not proof about recovery
  inside `ghl_con2prim_multi_method`.

## Sources

- [Purpose statement](../../../GRHayLHD/README)
- [Shared parameter declarations](../../../GRHayLHD/param.ccl)
- [Shared header](../../../GRHayLHD/src/GRHayLHD.h)
- [MoL registration](../../../GRHayLHD/src/MoL_registration.c)
- [Hybrid recovery](../../../GRHayLHD/src/Hybrid/conservs_to_prims.c)
- [TabulatedEntropy fluxes](../../../GRHayLHD/src/TabulatedEntropy/evaluate_fluxes_rhs.c)
- [Stress-energy computation](../../../GRHayLHD/src/compute_Tmunu.c)

## Related Pages

- [EOS and Entropy Variants](../evolution/eos-entropy-variants.md)
- [Conservative Recovery](../evolution/conservative-recovery.md)
- [RHS Fluxes and Sources](../evolution/rhs-fluxes-and-sources.md)
