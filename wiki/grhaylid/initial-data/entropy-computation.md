# Entropy Computation

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns local entropy-selection declarations and visible implementation
for Hybrid and Tabulated EOS branches. GRHayLib EOS semantics, bounds behavior,
computed values, successful scheduling, and behavior for externally defined
EOS selections remain out of scope.

## Summary

`HydroBase::initial_entropy="GRHayLID"` activates a schedule-level dispatch.
`EOS_type="Hybrid"` declares the hybrid entropy routine; `EOS_type="Tabulated"`
declares the tabulated routine. No third branch or error arm is declared. The
hybrid body writes entropy from density and pressure. The tabulated body first
passes density, electron fraction, and temperature by address to a bounds
helper, then calls a pressure/energy/entropy helper that writes pressure,
specific internal energy, and entropy.

## Mode Applicability

| Applicability | Visible behavior |
| --- | --- |
| Entropy/Hybrid | Schedule selects `GRHayLID_compute_entropy_hybrid`; body reads density and pressure and writes entropy. |
| Entropy/Tabulated | Schedule selects `GRHayLID_compute_entropy_tabulated`; body may receive in-place changes to density, electron fraction, and temperature and writes pressure, energy, and entropy through an external call. |
| Common | README describes entropy computation as usable with any initial-data thorn; only Hybrid and Tabulated schedule arms are locally declared. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-ENT-01` | HydroBase's `initial_entropy` keyword is extended with the value `GRHayLID`. | declared | Shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` |
| `ID-ENT-02` | Hybrid implementation visibly writes entropy from density and pressure through a GRHayLib call. | visible-implementation | Hybrid entropy function | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_hybrid` |
| `ID-ENT-03` | Tabulated implementation visibly passes density, electron fraction, and temperature by address to a bounds helper, then passes pressure, energy, and entropy output addresses to an EOS helper. | visible-implementation | Tabulated entropy function | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` |
| `ID-ENT-04` | Schedule CCL declares the Hybrid entropy arm after initial-data and beta-equilibrium aliases and before primitive-to-conservative initialization. | declared | Hybrid schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` |
| `ID-ENT-05` | Schedule CCL declares the Tabulated entropy arm with the same ordering and broader writes. | declared | Tabulated schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` |
| `ID-ENT-06` | README's “any EOS” wording is broader than the two locally declared EOS dispatch arms. | unresolved | Purpose wording compared with schedule declarations | `doc:GRHayLID/README#section=1. Purpose` |

## Details

### Selection contract

The outer condition compares `initial_entropy` with `GRHayLID`. Inside it, one
arm compares `EOS_type` with `Hybrid` and one `else if` compares it with
`Tabulated`. The CCL contains no final `else` or diagnostic. Therefore, for a
selected entropy keyword and any other EOS value, this local schedule block
declares no entropy routine. “Schedules nothing” here describes visible CCL
only; external parameter domains and runtime behavior are not inferred. See
[GID-0010](../contradictions.md#gid-0010).

### Hybrid visible dataflow

For each grid point, `GRHayLID_compute_entropy_hybrid` passes `ghl_eos`,
`rho[index]`, and `press[index]` to
`ghl_hybrid_compute_entropy_function`, assigning the returned value to
`entropy[index]`. Schedule CCL correspondingly declares reads of density and
pressure and a write to entropy.

### Tabulated visible dataflow

For each point, `GRHayLID_compute_entropy_tabulated` calls
`ghl_tabulated_enforce_bounds_rho_Ye_T` with addresses of `rho`, `Y_e`, and
`temperature`, making in-place mutation visible at the call boundary. It then
calls `ghl_tabulated_compute_P_eps_S_from_T` with their current values and
addresses for `press`, `eps`, and `entropy`. Schedule CCL declares all six as
writes except that the first three are also inputs. Result semantics remain
external.

### Documentation disagreement

README says entropy can be computed “for any EOS,” while CCL names only Hybrid
and Tabulated. ThornGuide says entropy is controlled by a `compute_entropy`
parameter, but no such local parameter is declared; the visible control is the
`initial_entropy` extension plus `EOS_type`. These are tracked as
[GID-0009](../contradictions.md#gid-0009) and
[GID-0001](../contradictions.md#gid-0001).

## Caveats

- Local calls do not establish GRHayLib formula, bounds, or error behavior.
- OpenMP loop pragmas do not prove parallel execution or thread safety.
- Schedule declarations do not prove that external bins/groups exist or run.
- No checked-in test or oracle observes either entropy branch.

## Sources

- [Purpose statement](../../../GRHayLID/README)
- [Keyword extensions](../../../GRHayLID/param.ccl)
- [Entropy schedules](../../../GRHayLID/schedule.ccl)
- [Entropy implementations](../../../GRHayLID/src/ComputeEntropy.c)

## Related Pages

- [Beta Equilibrium](beta-equilibrium.md)
- [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [HydroBase Keyword Extensions](../integration/hydrobase-keyword-extensions.md)
- [GID-0001: documented control mismatch](../contradictions.md#gid-0001)
- [GID-0009: any-EOS wording](../contradictions.md#gid-0009)
- [GID-0010: undeclared fallback arm](../contradictions.md#gid-0010)
