# EOS and Entropy Variants

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page owns canonical mapping from local `EOS_type` and `evolve_entropy`
dispatch to storage, registration, state, seven operation families, recovery
fallback, and checked-in test selection. GRHayLib initialization and EOS
semantics remain external.

## Summary

Schedule CCL declares four branches: Hybrid or Simple with entropy off/on, and
Tabulated with entropy off/on. Entropy adds `ent_star`; tabulated mode adds
`Ye_star`; both axes combine independently. Authored Balsara0 and TOV test
inputs set Simple and Hybrid respectively but do not visibly set shared
`evolve_entropy`, whose default is not declared in local `param.ccl`.

## Variant Applicability

This is canonical four-mode matrix. Each operation cell gives exact scheduled
symbol suffixes following row prefix; full symbols are `GRHayLHD_` + prefix +
suffix.

| Applicability | Local dispatch / source family | Conditional storage and MoL registration | Primitive, evolved, RHS, and flux state | Seven scheduled operation symbols | Recovery fallback | Local test status |
| --- | --- | --- | --- | --- | --- | --- |
| Hybrid/Simple | `EOS_type` Hybrid or Simple; entropy false; `Hybrid` | Core storage; calls evolved-group API for core/RHS and constrained-group API for `rho`, `press`, `eps` | Base thermodynamics and velocity; five core conservatives/RHS/flux | Prefix `hybrid_`: `prims_to_conservs`, `conservs_to_prims`, `evaluate_sources_rhs`, `evaluate_fluxes_rhs`, `outer_boundaries`, `perturb_primitives`, `perturb_conservatives` | Visible explicit `ghl_hybrid_Font1D` after weighted retries | Balsara0 sets Simple; TOV sets Hybrid; entropy selection unresolved locally |
| Hybrid/Simple+Entropy | Same EOS choices; entropy true; `HybridEntropy` | Add `ent_star[3]`, RHS, flux; calls APIs with entropy evolved/RHS and HydroBase entropy indices | Add entropy primitive, evolved, RHS, flux, reconstruction, and advection state | Prefix `hybrid_entropy_` with same seven exact suffixes | Visible explicit `ghl_hybrid_Font1D` after weighted retries | No authored test explicitly selects entropy locally |
| Tabulated | `EOS_type` Tabulated; entropy false; `Tabulated` | Add `Ye_star[3]`, RHS, flux; calls APIs with electron-fraction evolved/RHS and `Y_e`/temperature indices | Add `Y_e`, temperature, `Ye_star`, RHS, flux, reconstruction, and advection state | Prefix `tabulated_` with same seven exact suffixes | No equivalent explicit Font1D call visible | No authored test explicitly selects Tabulated locally |
| Tabulated+Entropy | Tabulated; entropy true; `TabulatedEntropy` | Adds both optional storage sets and both optional group-index call sets | Adds entropy and electron-fraction primitive/evolved/RHS/flux/reconstruction/advection state | Prefix `tabulated_entropy_` with same seven exact suffixes | No equivalent explicit Font1D call visible | No authored test explicitly selects this combination locally |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-MODE-01` | Parameter CCL uses shared `EOS_type`. | declared | Shared parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=EOS_type` |
| `EV-MODE-02` | Parameter CCL uses shared `evolve_entropy` without a local default. | declared | Shared parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=evolve_entropy` |
| `EV-MODE-03` | Schedule declares Hybrid/Simple entropy dispatch and concrete operation symbols. | declared | Concrete Prim2Con schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_entropy_prims_to_conservs` |
| `EV-MODE-04` | Schedule declares Tabulated entropy dispatch and concrete operation symbols. | declared | Concrete Prim2Con schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_tabulated_entropy_prims_to_conservs` |
| `EV-MODE-05` | Function visibly calls MoL APIs with optional group indices on independent entropy and tabulated conditions. | visible-implementation | Registration branches | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` |
| `EV-MODE-06` | Balsara0 authored test input sets Simple EOS. | checked-in-observation | Exact assignment | `par:GRHayLHD/test/Balsara0.par#parameter=GRHayLib::EOS_type` |
| `EV-MODE-07` | TOV authored test input sets Hybrid EOS. | checked-in-observation | Exact assignment | `par:GRHayLHD/test/TOV.par#parameter=GRHayLib::EOS_type` |
| `EV-MODE-08` | Balsara0 authored input contains no explicit `GRHayLib::evolve_entropy` assignment. | checked-in-observation | Full authored-input inspection | `par:GRHayLHD/test/Balsara0.par#file` |
| `EV-MODE-09` | TOV authored input contains no explicit `GRHayLib::evolve_entropy` assignment. | checked-in-observation | Full authored-input inspection | `par:GRHayLHD/test/TOV.par#file` |
| `EV-MODE-10` | Shared entropy selection remains locally unresolved because authored-input omissions do not establish the external shared parameter value. | unresolved | Shared selection relationship | `ccl:GRHayLHD/param.ccl#parameter=evolve_entropy` |

## Details

### Exact concrete symbols

For each prefix in canonical matrix, seven suffixes expand to exact CCL and C
symbols. Example Hybrid/Simple row expands to
`GRHayLHD_hybrid_prims_to_conservs`,
`GRHayLHD_hybrid_conservs_to_prims`,
`GRHayLHD_hybrid_evaluate_sources_rhs`,
`GRHayLHD_hybrid_evaluate_fluxes_rhs`,
`GRHayLHD_hybrid_outer_boundaries`,
`GRHayLHD_hybrid_perturb_primitives`, and
`GRHayLHD_hybrid_perturb_conservatives`. Other rows substitute their exact
prefix from matrix. All 28 concrete schedules were inspected in local CCL;
all 28 matching C units were inspected independently rather than inferred from
filename symmetry.

### Storage and registration

Core five conservatives and RHS always have storage. `EOS_type == Tabulated`
adds `Ye_star[3]`, `Ye_star_rhs`, and `Ye_star_flux`; `evolve_entropy` adds
`ent_star[3]`, `ent_star_rhs`, and `ent_star_flux`. Function always passes core
conservative/RHS indices to an evolved-group API and base HydroBase indices to
a constrained-group API. Entropy condition adds entropy group indices.
Tabulated condition adds electron-fraction evolved/RHS plus HydroBase electron-
fraction and temperature indices. External registration effects remain
unverified.

### Reconstruction and advection extras

HybridEntropy flux code reconstructs entropy and stores entropy flux.
Tabulated reconstructs electron fraction, applies tabulated bound/EOS helper
calls, and stores electron-fraction flux. TabulatedEntropy does both. Source
routines initialize matching optional RHS fields to zero; flux routines add
directional differences for each matching optional flux.

### Test-selection boundary

`GRHayLHD/test/Balsara0.par` sets `GRHayLib::EOS_type = "Simple"` and
`GRHayLHD/test/TOV.par` sets it to `"Hybrid"`. Neither authored file visibly
assigns `GRHayLib::evolve_entropy`. Because GRHayLHD `USES` rather than locally
declares that boolean's default, no local evidence selects entropy-on or
entropy-off for either test. Companion parfiles repeat EOS assignments but do
not establish runtime provenance or replace this authored-input observation.

## Caveats

- Hybrid and Simple share one local dispatch family; their external EOS
  differences are not reconstructed here.
- CCL conditions and GRHayLib runtime-object fields are locally parallel but
  their initialization equivalence is delegated.
- "No explicit test selection" is an absence observation, not proof of
  entropy-off behavior.
- Recovery helper semantics remain external; matrix records visible calls and
  absence of an equivalent explicit Font1D call only.

## Sources

- [Parameter declarations](../../../GRHayLHD/param.ccl)
- [Schedule dispatch and storage](../../../GRHayLHD/schedule.ccl)
- [MoL registration](../../../GRHayLHD/src/MoL_registration.c)
- [Balsara0 authored input](../../../GRHayLHD/test/Balsara0.par)
- [TOV authored input](../../../GRHayLHD/test/TOV.par)
- Four variant implementation families under [GRHayLHD/src](../../../GRHayLHD/src)

## Related Pages

- [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
- [Architecture Variables and Storage](../architecture/variables-and-storage.md)
- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
