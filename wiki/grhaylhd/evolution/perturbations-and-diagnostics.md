# Perturbations and Diagnostics

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page records primitive and conservative perturbation schedule positions,
guards, mode field sets, and visible RNG use. It does not establish RNG thread
safety, deterministic sequences, perturbation independence, statistical
quality, or diagnostic validity.

## Summary

`perturb_initial_data` schedules primitive perturbation after HydroBase velocity
conversion and before Prim2Con. `perturb_every_con2prim` schedules conservative
perturbation after conservative sync and before recovery. Every variant seeds
`rand()` once before an OpenMP loop and multiplies each selected field by an
independently invoked `one_plus_pert(random_pert)` expression.

## Variant Applicability

| Applicability | Primitive fields multiplied | Conservative fields multiplied |
| --- | --- | --- |
| Common | `rho`, `press`, `vx`, `vy`, `vz` | `rho_star`, `tau`, `Stildex`, `Stildey`, `Stildez` |
| Hybrid/Simple | Common set | Common set |
| Hybrid/Simple+Entropy | Common plus `entropy` | Common plus `ent_star` |
| Tabulated | Common plus `Y_e`, `temperature` | Common plus `Ye_star` |
| Tabulated+Entropy | Common plus `entropy`, `Y_e`, `temperature` | Common plus `ent_star`, `Ye_star` |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-PERT-01` | Initial primitive perturbation is disabled by default. | declared | Parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=perturb_initial_data` |
| `EV-PERT-02` | Per-recovery conservative perturbation is disabled by default. | declared | Parameter declaration | `ccl:GRHayLHD/param.ccl#parameter=perturb_every_con2prim` |
| `EV-PERT-03` | Hybrid primitive routine multiplies core primitive fields. | visible-implementation | Primitive function | `c:GRHayLHD/src/Hybrid/perturb_primitives.c#symbol=GRHayLHD_hybrid_perturb_primitives` |
| `EV-PERT-04` | Hybrid conservative routine multiplies five core conservatives. | visible-implementation | Conservative function | `c:GRHayLHD/src/Hybrid/perturb_conservatives.c#symbol=GRHayLHD_hybrid_perturb_conservatives` |
| `EV-PERT-05` | HybridEntropy primitive routine adds entropy. | visible-implementation | Primitive function | `c:GRHayLHD/src/HybridEntropy/perturb_primitives.c#symbol=GRHayLHD_hybrid_entropy_perturb_primitives` |
| `EV-PERT-06` | HybridEntropy conservative routine adds `ent_star`. | visible-implementation | Conservative function | `c:GRHayLHD/src/HybridEntropy/perturb_conservatives.c#symbol=GRHayLHD_hybrid_entropy_perturb_conservatives` |
| `EV-PERT-07` | Tabulated primitive routine adds electron fraction and temperature. | visible-implementation | Primitive function | `c:GRHayLHD/src/Tabulated/perturb_primitives.c#symbol=GRHayLHD_tabulated_perturb_primitives` |
| `EV-PERT-08` | Tabulated conservative routine adds `Ye_star`. | visible-implementation | Conservative function | `c:GRHayLHD/src/Tabulated/perturb_conservatives.c#symbol=GRHayLHD_tabulated_perturb_conservatives` |
| `EV-PERT-09` | TabulatedEntropy primitive routine includes entropy, electron fraction, and temperature. | visible-implementation | Primitive function | `c:GRHayLHD/src/TabulatedEntropy/perturb_primitives.c#symbol=GRHayLHD_tabulated_entropy_perturb_primitives` |
| `EV-PERT-10` | TabulatedEntropy conservative routine includes both optional evolved fields. | visible-implementation | Conservative function | `c:GRHayLHD/src/TabulatedEntropy/perturb_conservatives.c#symbol=GRHayLHD_tabulated_entropy_perturb_conservatives` |
| `EV-PERT-11` | Macro visibly calls `rand()`, but local evidence does not verify behavior under parallel calls. | coverage-gap | Perturbation macro | `macro:GRHayLHD/src/GRHayLHD.h#name=one_plus_pert` |
| `EV-PERT-12` | Initial perturbation guard schedules primitive perturbation after native-velocity conversion and before Prim2Con. | declared | Initial perturbation schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_perturb_primitives` |
| `EV-PERT-13` | Per-recovery guard schedules conservative perturbation after conservative sync and before recovery. | declared | Recovery perturbation schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_perturb_conservatives` |
| `EV-PERT-14` | Hybrid primitive-perturbation schedule declares `HydroBase::eps` in `READS` and `WRITES`. | declared | Hybrid schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_perturb_primitives` |
| `EV-PERT-15` | ThornGuide describes perturbation controls as debugging/testing aids. | declared | Parameters section | `doc:GRHayLHD/doc/documentation.tex#section=Parameters` |
| `EV-PERT-16` | HybridEntropy primitive-perturbation schedule declares `HydroBase::eps` in `READS` and `WRITES`. | declared | HybridEntropy schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_entropy_perturb_primitives` |
| `EV-PERT-17` | Tabulated primitive-perturbation schedule declares `HydroBase::eps` in `READS` and `WRITES`. | declared | Tabulated schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_tabulated_perturb_primitives` |
| `EV-PERT-18` | TabulatedEntropy primitive-perturbation schedule declares `HydroBase::eps` in `READS` and `WRITES`. | declared | TabulatedEntropy schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_tabulated_entropy_perturb_primitives` |
| `EV-PERT-19` | Hybrid primitive body does not reference `eps`. | visible-implementation | Hybrid primitive body | `c:GRHayLHD/src/Hybrid/perturb_primitives.c#symbol=GRHayLHD_hybrid_perturb_primitives` |
| `EV-PERT-20` | HybridEntropy primitive body does not reference `eps`. | visible-implementation | HybridEntropy primitive body | `c:GRHayLHD/src/HybridEntropy/perturb_primitives.c#symbol=GRHayLHD_hybrid_entropy_perturb_primitives` |
| `EV-PERT-21` | Tabulated primitive body does not reference `eps`. | visible-implementation | Tabulated primitive body | `c:GRHayLHD/src/Tabulated/perturb_primitives.c#symbol=GRHayLHD_tabulated_perturb_primitives` |
| `EV-PERT-22` | TabulatedEntropy primitive body does not reference `eps`. | visible-implementation | TabulatedEntropy primitive body | `c:GRHayLHD/src/TabulatedEntropy/perturb_primitives.c#symbol=GRHayLHD_tabulated_entropy_perturb_primitives` |

## Details

### Parameters and positions

`random_seed` locally permits 0 through 99999999 and defaults to 0.
`random_pert` permits any real and defaults to 0. Both perturbation booleans
default to no. ThornGuide describes these features as debugging or other
testing controls for initial or evolution data; this is documented intent, not
validation.

When initial perturbation is enabled, schedule group follows
`convert_HydroBase_to_GRHayLHD` and precedes `GRHayLHD_prims_to_conservs`.
Thus native velocity and HydroBase thermodynamic fields are perturbed before
conservative construction. When per-Con2Prim perturbation is enabled,
conservative group follows `GRHayLHD_sync_conservatives` and precedes recovery.

### Exact four-family field matrix

All eight routines were inspected independently. Primitive schedules declare
`HydroBase::eps` in both `READS` and `WRITES`, but none of four visible
primitive bodies references `eps`; declared and visible read/write sets differ
under [GRH-0012](../contradictions.md#grh-0012). Visible C bodies determine
the actual multiplication list above.
Conservative routines multiply core evolved fields and active optional evolved
fields only. No RHS, flux temporary, or `failure_checker` field is perturbed.

Each assignment is multiplicative, with a fresh expansion of
`one_plus_pert(random_pert)`. Macro is `1 + perturb * rand() / RAND_MAX` in
visible header. No distribution, sign range beyond expression, or independence
claim is inferred.

### RNG and OpenMP review boundary

Each routine calls `srand(random_seed)` before `#pragma omp parallel for`; macro
then calls `rand()` inside loop. Source comment states `rand()` is thread-safe,
but local tree contains no cited test, platform contract, synchronization
analysis, or reproducibility oracle for concurrent RNG use. This is retained as
coverage gap/caveat, not proof of safety, race, failure, or contradiction.

### Recovery diagnostics relationship

Conservative perturbations occur before recovery and can therefore feed visible
recovery diagnostics. Recovery writes `failure_checker` each point and later
recomputes conservatives. Perturbation code itself does not write diagnostic
gridfunction. Existing terminal failure-code overwrite ambiguity is owned by
[GRH-0004](../contradictions.md#grh-0004).

## Caveats

- Schedule positions are declarations, not observed executions.
- C library/OpenMP RNG behavior is external and unverified locally.
- Re-seeding on every routine call is visible; resulting repeatability under
  serial or parallel execution is not established.
- Primitive declaration/body read/write-set mismatch remains open
  [GRH-0012](../contradictions.md#grh-0012).
- Zero default magnitude makes multiplier expression one in ordinary arithmetic,
  but no runtime/no-effect claim is made.

## Sources

- [ThornGuide parameters](../../../GRHayLHD/doc/documentation.tex)
- [Parameter declarations](../../../GRHayLHD/param.ccl)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Perturbation macro](../../../GRHayLHD/src/GRHayLHD.h)
- Eight perturbation implementations under [variant source tree](../../../GRHayLHD/src)

## Related Pages

- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Conservative Recovery](conservative-recovery.md)
- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
