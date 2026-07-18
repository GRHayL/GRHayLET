# Coverage Gaps

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Validation](index.md)

## Scope and Non-Scope

This page classifies dimensions not visibly owned by local test declarations,
authored inputs, or checked-in oracles, then ranks proposed future evidence.
Absence means “not visibly covered,” never broken or unsupported. Proposals do
not change tests, source, configurations, or oracles.

## Summary

Checked-in cases select Simple and Hybrid EOS only, with entropy choice
unresolved. Oracles cover Balsara density/pressure/native velocity in x and TOV
density/pressure in x/y/z. No checked-in oracle owns optional state, internal
conservatives/RHS/fluxes, recovery diagnostics, several lifecycle boundaries,
or scientific assurance dimensions. Highest-risk proposals target four-mode
dispatch, cadence/leakage, recovery diagnostics, and Tmunu lifecycle.

## Variant Applicability

| Applicability | Visible gap |
| --- | --- |
| Common | No direct oracle ownership for conservatives, RHS, flux temporaries, `failure_checker`, cadence/leakage, Tmunu contribution, restart/timelevels, or long-horizon behavior. |
| Hybrid/Simple | Existing inputs select Simple/Hybrid, but recovery-failure, fallback, boundary-choice, perturbation, and directional-helper distinctions lack owned observations. |
| Hybrid/Simple+Entropy | No explicit authored selection or entropy-state oracle. |
| Tabulated | No explicit authored selection or `Y_e`/temperature/`Ye_star` oracle. |
| Tabulated+Entropy | No explicit authored selection or combined optional-state oracle. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `VAL-GAP-01` | Authored inputs visibly select Simple and Hybrid but do not explicitly select Tabulated or entropy evolution. | coverage-gap | Both authored EOS assignments plus shared declaration | `par:GRHayLHD/test/Balsara0.par#parameter=GRHayLib::EOS_type` |
| `VAL-GAP-02` | Eleven checked-in oracles omit direct ownership for optional evolved state and internal conservative/RHS/flux/diagnostic fields. | coverage-gap | Complete oracle inventory | `oracle:GRHayLHD/test/Balsara0/rho.x.asc#file` |
| `VAL-GAP-03` | Copy matter boundary is explicit in Balsara inputs; no authored input explicitly selects frozen boundary. | coverage-gap | Balsara boundary assignment | `par:GRHayLHD/test/Balsara0.par#parameter=GRHayLHD::Matter_BC` |
| `VAL-GAP-04` | Local checked-in evidence does not own positive cadence/leakage precondition behavior. | coverage-gap | Leakage schedule boundary | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=GRHayLHD_RHS` |
| `VAL-GAP-05` | Local evidence does not establish serial/parallel RNG behavior or reproducibility. | coverage-gap | RNG macro and perturbation call pattern | `macro:GRHayLHD/src/GRHayLHD.h#name=one_plus_pert` |
| `VAL-GAP-06` | Two single-tolerance case declarations provide no local multi-resolution convergence or conservation evidence. | coverage-gap | Balsara0 declaration | `test:GRHayLHD/test/test.ccl#case=Balsara0` |
| `VAL-GAP-07` | Companion and oracle assertions do not establish current-input identity, production, chronology, or a shared provenance chain. | unresolved | Balsara companion/oracle issue | `par:GRHayLHD/test/Balsara0/Balsara0.par#file` |

## Details

### Visible gap inventory

| Area | What checked-in evidence owns | Not visibly covered |
| --- | --- | --- |
| Mode dispatch | Simple Balsara0 and Hybrid TOV assignments; entropy omitted | Explicit HybridEntropy, Tabulated, and TabulatedEntropy selection; `ent_star`, `Ye_star`, `Y_e`, temperature ownership |
| Internal hydrodynamic state | HydroBase density/pressure and native velocity oracle files | Core conservatives, optional conservatives, every RHS group, flux temporaries, and `failure_checker` |
| Matter boundaries and symmetry | Copy explicitly selected in Balsara; only `none` locally selectable for symmetry | Explicit frozen/outflow distinction, inflow clipping, six-face behavior, AMR/coarsest-level gates, dormant equatorial branches |
| Perturbations and RNG | Disabled local defaults and visible perturbation source | Initial/per-recovery perturbation outputs, seed repeatability, serial/parallel RNG behavior, field-by-field magnitude evidence |
| HydroBase conversion and leakage | Conversion formulas and schedule declarations | Positive cadence behavior, zero-cadence leakage path, `w_lorentz` warning path, leakage integration |
| Tmunu | TOV configuration supplies TmunuBase settings; local function visibly adds components | Direct stress-energy observation, contribution ordering, steering transition, registration-error handling |
| Recovery | TOV config names external Font1D backup; recovery code has diagnostic paths | Forced primary failure/fallback, averaging retries, atmosphere terminal path, NaN screen, `failure_checker` legend/writeback distinction |
| Flux/source helpers | Visible interpolation/derivative formulas and directional calls | Face-helper focused observation, x/y/z directional consistency, source-before-flux runtime ordering |
| Lifecycle/storage | Interface tags and declared schedules; TOV contains multiple refinement levels | Restart/checkpoint behavior, three-timelevel evolution, synchronization, AMR boundary gating, prolongation behavior |
| Build/dependency | Local configuration requires HDF5 | Requirement/link behavior under build configurations |
| Scientific assurance | Two short declared cases and checked-in numeric content | Resolution series, convergence order, conservation budgets, long-term stability, and broader EOS/state ranges |

TOV's multiple configured refinement levels do not by themselves own AMR gate
behavior: checked-in files do not isolate boundary function entry/skip paths or
attribute output to those gates.

### Ranked proposals

Every row is a proposal. “Expected evidence” names an artifact needed to close
a gap; it does not predict a successful result.

| Rank / operational risk | Proposed smallest configuration | Expected evidence | Dependencies |
| --- | --- | --- | --- |
| P0 / state-layout divergence | Add one minimal input for each missing mode: Hybrid+entropy, Tabulated, Tabulated+entropy; request active optional primitives/conservatives and a core baseline | Authored assignments plus mode-specific ASCII fields showing field presence and iteration blocks | EOS tables for tabulated modes; external GRHayLib initialization; output support for internal fields |
| P0 / zero-divisor precondition | Minimal existing case with NRPyLeakageET active at default cadence, paired with explicit positive cadence | Scheduler/runtime record and HydroBase velocity/Lorentz outputs sufficient to distinguish each call path | NRPyLeakageET availability; diagnostic capture; source fix may be prerequisite if zero path is unsafe |
| P0 / recovery fallback and diagnostics | Small grid with controlled conservative states that trigger primary failure, weighted retry, explicit Hybrid fallback, and atmosphere path separately | `failure_checker`, corrected conservatives/primitives, and bounded diagnostic log for each targeted path | Reproducible fault injection or crafted state; external recovery semantics |
| P0 / Tmunu lifecycle | Minimal matter configuration with `update_Tmunu=no/yes`, plus a supported steering-transition probe | TmunuBase component output, registration status, and contribution timing evidence | MoL/TmunuBase lifecycle semantics; a defined supported steering contract |
| P1 / boundary distinctions | Small uniform-state grids selecting copy, outflow, and frozen separately; inject normal inflow on one face at a time | Six-face primitive/conservative observations with signed normal-velocity cases | Boundary output including ghost zones; one coarse level first, then AMR gate case |
| P1 / RNG and perturbations | Tiny serial and multi-thread grids with fixed seed, nonzero magnitude, and one perturbation guard at a time | Field-by-field outputs and repeat-run comparison for serial/parallel call patterns | Platform C/OpenMP RNG contract or instrumentation; controlled thread count |
| P1 / directional helpers | Symmetric smooth-data inputs rotated across x/y/z and focused face/helper output | Direction-paired flux/RHS fields and face values | Instrumentation or output registration for temporaries; external PPM/HLLE semantics remain separate |
| P1 / restart and AMR lifecycle | Short checkpoint/restart pair with multiple levels and boundary gating, compared at same logical iteration | Restarted evolved/primitive state, timelevel ownership, and gate-specific diagnostics | Checkpoint-capable driver; exact restart protocol; internal-field output |
| P2 / dependency surface | Minimal configuration/build matrix varying HDF5 availability within supported build process | Build-system evidence identifying declared requirement and link boundary | Authorized build environment; external HDF5/Cactus configuration |
| P2 / scientific assurance | Resolution series for shock and star cases, then longer-duration variants | Convergence measures, conservation budgets, and time-series stability artifacts | Agreed norms/reference quantities, compute budget, external numerical-method expectations |

### Provenance and RNG classification

Balsara oracle headers name one parameter path; Balsara companion names another
original path and asserts Cactus generation. TOV companion makes analogous
generation/original-path assertions, while TOV oracle headers contain only a
generic CarpetIOASCII label. None proves identity with current authored input,
companion-to-oracle production, artifact chronology, or a shared chain. These
remain provenance ambiguities [GRH-0009](../contradictions.md#grh-0009) and
[GRH-0010](../contradictions.md#grh-0010), not inferred lineage.

RNG source comment asserts `rand()` thread safety, while visible routines seed
once and call `rand()` inside OpenMP loops. Local evidence neither verifies nor
refutes runtime semantics. This remains caveat and focused coverage gap, not a
contradiction.

## Caveats

- Gap ranking reflects operational impact of unobserved boundaries, not proof
  of defects.
- Proposed configurations may require external thorns, EOS data, framework
  semantics, instrumentation, or source changes outside this KB task.
- Checked-in absence can change when files are added; registry and source-map
  reconciliation must precede status changes.
- No proposal was implemented, built, or run during KB review.

## Sources

- [Test declarations](../../../GRHayLHD/test/test.ccl)
- [Authored Balsara0 input](../../../GRHayLHD/test/Balsara0.par)
- [Authored TOV input](../../../GRHayLHD/test/TOV.par)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Perturbation macro](../../../GRHayLHD/src/GRHayLHD.h)
- [Balsara0 checked-in artifacts](../../../GRHayLHD/test/Balsara0)
- [TOV checked-in artifacts](../../../GRHayLHD/test/TOV)

## Related Pages

- [Test Inventory and Oracles](test-inventory-and-oracles.md)
- [EOS and Entropy Variants](../evolution/eos-entropy-variants.md)
- [Perturbations and Diagnostics](../evolution/perturbations-and-diagnostics.md)
- [HydroBase Velocity Conversion](../integration/hydrobase-velocity-conversion.md)
- [ADM, MoL, and Tmunu Contracts](../integration/adm-mol-tmunu-contracts.md)
