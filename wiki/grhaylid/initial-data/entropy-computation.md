# Entropy Computation

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

initial_entropy=GRHayLID selects a native-proxy routine for Simple/Hybrid
or a table routine for Tabulated. Native-proxy use requires explicit
allow_native_entropy_proxy=yes. Tabulated entropy stages doubles, checks
finite initialized inputs and EOS outputs, and commits the complete tuple
only after success. With beta enabled it preserves rho/Ye/T.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Entropy/Hybrid | Simple/Hybrid native proxy requires compatible-consumer opt-in. |
| Entropy/Tabulated | Physical table entropy with deliberate primitive replacement. |
| Common | A producer must initialize required fields before these routines. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ENTROPY-COMPUTATION-01` | The entropy selector is a HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` |
| `ENTROPY-COMPUTATION-02` | The native body requires opt-in, positive finite rho, nonnegative finite pressure, and finite helper output. | visible-implementation | Named local source | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_hybrid` |
| `ENTROPY-COMPUTATION-03` | The Tabulated body checks storage, stages input/output doubles, and checks EOS status before publication. | visible-implementation | Named local source | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` |
| `ENTROPY-COMPUTATION-04` | The native schedule arm permits Simple and Hybrid. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` |
| `ENTROPY-COMPUTATION-05` | README declares the native proxy representation and consumer limitation. | declared | Named local source | `doc:GRHayLID/README#section=1. Purpose` |
| `ENTROPY-COMPUTATION-06` | allow_native_entropy_proxy defaults to no. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=allow_native_entropy_proxy` |

## Details

Entropy runs after HydroBase_Initial and the beta alias, before
HydroBase_Prim2ConInitial. ParamCheck diagnoses unsupported EOS selections,
missing hydro selection, missing Tabulated Ye/T producers without beta, and
native-proxy selection without opt-in. Storage queries separately check the
actual entropy destination and Tabulated optional groups.

README declares Simple/Hybrid output as P/rho^(Gamma_piece-1), intended for
native GRHayL recovery rather than physical k_b/baryon entropy. The call
semantics remain delegated. Explicit opt-in records a compatible-consumer
contract; it does not convert the proxy into physical entropy.

Tabulated entropy rejects nonfinite input triples. Without beta, it passes
local double temporaries to the bounds helper, validates effective bounds,
and obtains pressure/energy/entropy into separate doubles. All six fields
are published after successful status and finite outputs. With beta, it
validates effective bounds without clamping, preserving density, Ye, and T
while recomputing thermodynamics. Allocation does not prove initialized input
ownership; external producers must be correctly ordered.

The guide now names the actual selector and documents side effects;
GID-0001/0009/0010 record resolved control and dispatch discrepancies.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [param.ccl](../../../GRHayLID/param.ccl)
- [ComputeEntropy.c](../../../GRHayLID/src/ComputeEntropy.c)
- [schedule.ccl](../../../GRHayLID/schedule.ccl)
- [README](../../../GRHayLID/README)

## Related Pages

- [Beta Equilibrium](beta-equilibrium.md)
- [Keyword Extensions](../integration/hydrobase-keyword-extensions.md)
- [GID-0001](../contradictions.md#gid-0001)
- [GID-0009](../contradictions.md#gid-0009)
- [GID-0010](../contradictions.md#gid-0010)
