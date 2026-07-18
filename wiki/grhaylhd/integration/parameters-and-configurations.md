# Parameters and Configurations

> Page status: reviewed · Last reviewed: 07-18-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page maps exactly nine locally declared GRHayLHD parameters, two shared
GRHayLib selections used locally, and five checked-in parfiles. It records
declarations, visible consumers, and checked-in assignments. External
parameter semantics, active-thorn internals, execution, and configuration-to-
oracle provenance remain out of scope.

## Summary

Local CCL declares cadence, Tmunu, symmetry, matter-boundary, verbosity, and
four perturbation controls. It uses shared `EOS_type` and `evolve_entropy`
without redeclaring their defaults. Five parfiles have three distinct roles:
one shipped example, two authored test inputs, and two checked-in companion
configurations whose headers assert Cactus generation. Only EOS, matter
boundary, and Tmunu update receive explicit assignments among these eleven
selections.

## Variant Applicability

| Applicability | Parameter interaction |
| --- | --- |
| Common | Nine local controls apply without an EOS-family requirement; cadence, Tmunu, boundary, verbosity, and perturbation guards own common lifecycle choices. |
| Hybrid/Simple | `EOS_type` selects Hybrid/Simple source family; no optional state unless entropy is selected. |
| Hybrid/Simple+Entropy | `evolve_entropy` selects entropy storage, registration, and concrete HybridEntropy schedules. |
| Tabulated | `EOS_type` selects Tabulated storage, registration, and schedules. |
| Tabulated+Entropy | Both shared selections enable both optional state families. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-PAR-01` | Parameter CCL declares nine local controls. | declared | Local cadence declaration and complete file inspection | `ccl:GRHayLHD/param.ccl#parameter=Convert_to_HydroBase_every` |
| `INT-PAR-02` | Parameter CCL uses two shared GRHayLib selections without local defaults. | declared | Shared selection declaration | `ccl:GRHayLHD/param.ccl#parameter=evolve_entropy` |
| `INT-PAR-03` | Shipped example sets Simple EOS and copy matter boundary. | checked-in-observation | Example assignments | `par:GRHayLHD/par/Balsara0.par#parameter=GRHayLib::EOS_type` |
| `INT-PAR-04` | Authored Balsara0 input sets `update_Tmunu=no`, Simple EOS, and copy matter boundary. | checked-in-observation | Authored test input | `par:GRHayLHD/test/Balsara0.par#parameter=GRHayLHD::update_Tmunu` |
| `INT-PAR-05` | Authored TOV input sets Hybrid EOS. | checked-in-observation | Exact authored assignment | `par:GRHayLHD/test/TOV.par#parameter=GRHayLib::EOS_type` |
| `INT-PAR-06` | Balsara0 companion contains checked-in Cactus-generation and original-path header assertions. | checked-in-observation | Whole companion file | `par:GRHayLHD/test/Balsara0/Balsara0.par#file` |
| `INT-PAR-07` | Authored TOV input omits explicit assignments for all nine local GRHayLHD controls. | checked-in-observation | Full authored-input inspection | `par:GRHayLHD/test/TOV.par#file` |
| `INT-PAR-08` | TOV companion contains checked-in Cactus-generation and original-path header assertions. | checked-in-observation | Whole companion file | `par:GRHayLHD/test/TOV/TOV.par#file` |

## Details

### Exact parameter inventory

“Not marked” below means local declaration contains no `STEERABLE` attribute;
no external default semantics are inferred.

| Parameter | Local declaration, domain, default | Steerability | Schedule guards and visible C consumers | Mode interaction | Assignments across five parfiles / caveat |
| --- | --- | --- | --- | --- | --- |
| `Convert_to_HydroBase_every` | `INT`; `0:*`; default `0` | Not marked | Guards reverse conversion at initial and analysis phases; reverse converter uses it as modulo divisor; leakage call is outside positive guard | Common | No explicit assignment; local default disables guarded contexts, but leakage precondition remains [GRH-0006](../contradictions.md#grh-0006). |
| `update_Tmunu` | `CCTK_BOOLEAN`; default `yes` | `ALWAYS` | Guards `AddToTmunu` schedule and Tmunu constrained registration | Common | Authored and companion Balsara0 set `no`; other three omit it. Steering lifecycle is [GRH-0007](../contradictions.md#grh-0007). |
| `Symmetry` | `KEYWORD`; permitted `none`; default `none` | Not marked | Read by symmetry initialization, all recovery symmetry branches, and all outer-boundary routines | Common | No explicit assignment. Equatorial source branches are not locally selectable; see [GRH-0001](../contradictions.md#grh-0001). |
| `Matter_BC` | `KEYWORD`; `copy`, `outflow`, `frozen`; default `outflow` | Not marked | Read by all four outer-boundary routines; outer-boundary group itself is not parameter-guarded | Common | Example, authored Balsara0, and Balsara0 companion set `copy`; TOV files omit it. |
| `verbose` | `KEYWORD`; `no`, `yes`; default `yes` | `ALWAYS` | Read by all four recovery routines for aggregate messages | Common | No explicit assignment; logging and counters do not prove solver success. |
| `random_seed` | `INT`; `0:99999999`; default `0` | Not marked | Passed to `srand` in all eight perturbation routines when their containing schedule group runs | Common | No explicit assignment; concurrency behavior remains unverified. |
| `random_pert` | `REAL`; `*:*`; default `0` | Not marked | Passed to `one_plus_pert` in all eight perturbation routines | Common | No explicit assignment; unrestricted local range does not establish numerically safe magnitudes. |
| `perturb_initial_data` | `CCTK_BOOLEAN`; default `no` | Not marked | CCL condition schedules primitive perturbation between forward conversion and Prim2Con; no direct C read | Dispatches one of four primitive perturbation functions | No explicit assignment. |
| `perturb_every_con2prim` | `CCTK_BOOLEAN`; default `no` | Not marked | CCL condition schedules conservative perturbation after sync and before recovery; no direct C read | Dispatches one of four conservative perturbation functions | No explicit assignment. |
| `EOS_type` | `USES KEYWORD`; domain/default not locally declared; schedule compares `Hybrid`, `Simple`, `Tabulated` | External declaration not inspected | Controls conditional storage and four schedule families; C registration checks `ghl_eos->eos_type` | Primary EOS-family axis | Example/authored/companion Balsara0 set `Simple`; authored/companion TOV set `Hybrid`. External EOS semantics remain unverified. |
| `evolve_entropy` | `USES CCTK_BOOLEAN`; default not locally declared | External declaration not inspected | Controls entropy storage and schedule branches; C registration checks `ghl_params->evolve_entropy` | Entropy axis independent of EOS family | No explicit assignment in any five parfiles; selected value remains locally unresolved. |

### Five configuration roles

| File | Evidence role | Checked-in local observations |
| --- | --- | --- |
| `GRHayLHD/par/Balsara0.par` | Shipped example | Time-terminated Balsara-labeled configuration; activates GRHayLib, GRHayLHD, and GRHayLID; sets Simple EOS and copy matter boundary; requests generated parfile output. |
| `GRHayLHD/test/Balsara0.par` | Authored test input | Iteration-10 configuration; sets `update_Tmunu=no`, Simple EOS, copy matter boundary; requests five x-direction output variables and generated parfile output. |
| `GRHayLHD/test/TOV.par` | Authored test input | Iteration-2 TOV configuration; sets Hybrid EOS, names external Font1D as first Con2Prim backup, and requests HydroBase density/pressure plus generated parfile output. |
| `GRHayLHD/test/Balsara0/Balsara0.par` | Checked-in companion | Header asserts automatic generation by Cactus 4.15.0 and names an original Balsara0 input path; captures Balsara assignments and one `ActiveThorns` list. |
| `GRHayLHD/test/TOV/TOV.par` | Checked-in companion | Header asserts automatic generation by Cactus 4.15.0 and names an original TOV input path; captures TOV assignments and one `ActiveThorns` list. |

Active-thorn names show declared configuration only. For example, Balsara
files activate `GRHayLID`; this page does not import its implementation or use
it as GRHayLHD authority. Likewise, TOV's Font1D assignment does not establish
external recovery selection or behavior.

### Provenance boundary

Companion comments are checked-in assertions about Cactus generation and an
original path. They do not prove byte/content identity with current authored
inputs, establish that companions produced checked-in numeric oracles, order
all artifacts chronologically, or establish a shared provenance chain. See
[GRH-0009](../contradictions.md#grh-0009) and
[GRH-0010](../contradictions.md#grh-0010).

## Caveats

- Local defaults are declaration facts, not observed runtime selections.
- Omitted assignments do not reveal defaults owned by external shared
  parameters.
- Schedule guards establish declared dispatch, not successful execution.
- No test was run and no oracle was regenerated during this review.

## Sources

- [Parameter declarations](../../../GRHayLHD/param.ccl)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Shipped Balsara0 example](../../../GRHayLHD/par/Balsara0.par)
- [Authored Balsara0 input](../../../GRHayLHD/test/Balsara0.par)
- [Authored TOV input](../../../GRHayLHD/test/TOV.par)
- [Balsara0 companion](../../../GRHayLHD/test/Balsara0/Balsara0.par)
- [TOV companion](../../../GRHayLHD/test/TOV/TOV.par)

## Related Pages

- [EOS and Entropy Variants](../evolution/eos-entropy-variants.md)
- [Perturbations and Diagnostics](../evolution/perturbations-and-diagnostics.md)
- [Matter Boundaries and Symmetry](../evolution/matter-boundaries-and-symmetry.md)
- [HydroBase Velocity Conversion](hydrobase-velocity-conversion.md)
- [ADM, MoL, and Tmunu Contracts](adm-mol-tmunu-contracts.md)
