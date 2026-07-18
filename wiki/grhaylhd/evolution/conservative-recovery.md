# Conservative Recovery

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page traces visible conservative-to-primitive control flow and final
writes in all four variants. GRHayLib solver math, success criteria, atmosphere
values, and physical validity remain external.

## Summary

All variants route non-positive conserved density to atmosphere, otherwise
undensitize and call `ghl_con2prim_multi_method`, screen recovered fields for
NaN, retry failures with bounded-neighborhood weighted conservative inputs,
and use atmosphere after terminal failure. Hybrid families alone visibly call
an explicit Font1D fallback. After recovery, all variants limit primitives,
write them, then use a separate loop to recompute and overwrite conservatives
while accumulating change diagnostics.

## Variant Applicability

| Applicability | Recovered extras | Visible pre-solver conservative limits | Explicit post-retry fallback | Recomputed extras |
| --- | --- | --- | --- | --- |
| Common | Base thermodynamics and velocity | Mode-dependent | Atmosphere after terminal error | Five core conservatives |
| Hybrid/Simple | None | Yes | `ghl_hybrid_Font1D`, then atmosphere | None |
| Hybrid/Simple+Entropy | Entropy | Yes | `ghl_hybrid_Font1D`, then atmosphere | `ent_star` |
| Tabulated | `Y_e`, temperature | No local `ghl_apply_conservative_limits` call | Atmosphere; no explicit Font1D call | `Ye_star` |
| Tabulated+Entropy | Entropy, `Y_e`, temperature | No local `ghl_apply_conservative_limits` call | Atmosphere; no explicit Font1D call | `ent_star`, `Ye_star` |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-C2P-01` | Hybrid recovery visibly applies conservative limits, retries weighted neighborhoods, calls Font1D, writes primitives, and recomputes core conservatives. | visible-implementation | Full variant function | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` |
| `EV-C2P-02` | HybridEntropy follows Hybrid fallback structure and includes entropy in recovery and recomputation. | visible-implementation | Full variant function | `c:GRHayLHD/src/HybridEntropy/conservs_to_prims.c#symbol=GRHayLHD_hybrid_entropy_conservs_to_prims` |
| `EV-C2P-03` | Tabulated recovery includes electron fraction and temperature, but no explicit Font1D call. | visible-implementation | Full variant function | `c:GRHayLHD/src/Tabulated/conservs_to_prims.c#symbol=GRHayLHD_tabulated_conservs_to_prims` |
| `EV-C2P-04` | TabulatedEntropy includes both optional state sets and no explicit Font1D call. | visible-implementation | Full variant function | `c:GRHayLHD/src/TabulatedEntropy/conservs_to_prims.c#symbol=GRHayLHD_tabulated_entropy_conservs_to_prims` |
| `EV-C2P-05` | Interface declares `failure_checker` as per-substep-overwritten diagnostic. | declared | Diagnostic group | `ccl:GRHayLHD/interface.ccl#group=failure_checker` |

## Details

### Primary path and atmosphere path

Each function initializes diagnostics, metric, auxiliaries, zero magnetic
primitive components, and mode conservatives. Condition `cons.rho > 0.0`
enters solver path; its complement sets constant atmosphere, adds one to local
failure code, increments density-fix count, and marks success. Thus visible
atmosphere gate is non-positive density, although nearby decoder comment says
`rho_star < 0`.

Positive-density Hybrid paths call `ghl_apply_conservative_limits` before
undensitization. Tabulated paths visibly skip that helper. All call
`ghl_undensitize_conservatives` and `ghl_con2prim_multi_method`. A product of
required recovered primitive fields is tested with `isnan`; optional entropy,
electron fraction, and temperature join product in relevant modes. NaN forces
local singular-error code.

### Bounded-neighborhood retries

On error, code bounds each coordinate to local grid and scans clipped
`i-1..i+1`, `j-1..j+1`, `k-1..k+1` neighborhood, excluding center. It sums
every active conservative field. A loop with `avg_weight` values 1 through 4
constructs weighted neighbor/center combinations, divides by adjusted
`n_avg`, undensitizes, retries multi-method recovery, and repeats NaN screen.
This is visible algorithmic structure; no numerical-quality or race-free claim
is inferred from comments.

Hybrid and HybridEntropy then reload/apply conservative limits, undensitize,
and call `ghl_hybrid_Font1D`. If that or its NaN screen fails, they set
atmosphere. Tabulated and TabulatedEntropy proceed directly from exhausted
weighted retries to atmosphere. All terminal paths update aggregate failure
and horizon counters based on visible conditions.

### Post-recovery writes and recomputation

All paths call `ghl_enforce_primitive_limits_and_compute_u0` and abort through
helper on returned error. They write `rho`, `press`, `eps`, `u0`, and native
velocity. Entropy and tabulated families additionally write their mode
primitives.

A second OpenMP loop reconstructs metric and primitives, again zeros `BU`,
calls `ghl_compute_conservs`, and overwrites five core conservative fields plus
active `ent_star` and/or `Ye_star`. Before overwrite it snapshots original
conservatives; after computation it accumulates absolute differences and
denominators. Under `verbose == yes`, each variant prints aggregate backup,
limit, failure, iteration, and conservative-difference summaries; exact fields
follow mode.

### `failure_checker` legend versus write order

Source comments assign 1 to a density-atmosphere reset, 10 to speed limiting,
100 to "Both C2P and Font Fix failed", 1000 to backup use, 10000 to a tau
fix, and 100000 to a momentum fix. Tabulated and TabulatedEntropy repeat the
Font-Fix wording although neither file contains an explicit local Font1D call;
only Hybrid families visibly call `ghl_hybrid_Font1D`. This variant/legend
mismatch is tracked under [GRH-0014](../contradictions.md#grh-0014).

Separately, all four terminal branches execute
`failure_checker[index] += 100`, but later each point assigns
`failure_checker[index] = local_failure_checker + ...` without that terminal
100 contribution; this overwrite mismatch remains tracked under
[GRH-0004](../contradictions.md#grh-0004). No runtime diagnostic value is
asserted for either mismatch.

### Dormant symmetry-name mismatch

Hybrid uses `GRHayLHD::grmhd_conservatives` in its equatorial block. The other
three variants visibly request `GRHayLHD::grhd_conservatives`, which differs
from interface declaration. Only `Symmetry=none` is locally selectable, so
this remains dormant mismatch rather than supported symmetry behavior.

## Caveats

- Helper calls establish visible order, not solver semantics or recovery
  correctness.
- Absence of explicit Font1D in tabulated files does not exclude fallback
  inside external multi-method implementation.
- Comments calling second loop deterministic do not prove thread safety.
- Group-name mismatch: [GRH-0002](../contradictions.md#grh-0002).
- Failure-code overwrite mismatch: [GRH-0004](../contradictions.md#grh-0004).
- Code-100 legend mismatch: [GRH-0014](../contradictions.md#grh-0014).

## Sources

- [Interface diagnostic declaration](../../../GRHayLHD/interface.ccl)
- [Hybrid recovery](../../../GRHayLHD/src/Hybrid/conservs_to_prims.c)
- [HybridEntropy recovery](../../../GRHayLHD/src/HybridEntropy/conservs_to_prims.c)
- [Tabulated recovery](../../../GRHayLHD/src/Tabulated/conservs_to_prims.c)
- [TabulatedEntropy recovery](../../../GRHayLHD/src/TabulatedEntropy/conservs_to_prims.c)

## Related Pages

- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
- [Architecture Variables and Storage](../architecture/variables-and-storage.md)
