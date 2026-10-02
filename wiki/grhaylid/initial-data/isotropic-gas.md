# Isotropic Gas Initial Data

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

IsotropicGas requires Tabulated EOS and both GRHayLID Ye/T selectors.
The body checks active Ye/T/metric storage and finite effective-bounds input
triples. EOS outputs use double temporaries, return codes and output finiteness
are checked before the grid loop, and each point's metric is checked.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| IsotropicGas | Tabulated gas/sphere initialization with checked inputs and EOS results. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ISOTROPIC-GAS-01` | IsotropicGas checks selectors, storage, finite effective bounds, EOS return codes, and outputs. | visible-implementation | Named local source | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` |
| `ISOTROPIC-GAS-02` | ThornGuide describes the supported three-dimensional family and Tabulated selector. | declared | Named local source | `doc:GRHayLID/doc/documentation.tex#section=IsotropicGas` |
| `ISOTROPIC-GAS-03` | Schedule description names three-dimensional data and declares metric reads and primitive writes. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_IsotropicGas` |

## Details

Gas evaluates one primitive triple; sphere evaluates interior and exterior
triples independently before entering the loop. Density and temperature must
be strictly positive, all inputs finite, and each triple within the configured
effective rho/Ye/T bounds. A failed lookup is diagnosed with its region and
requested triple; neither output is published after a failed EOS status.

The loop publishes the uniform requested triple and computed pressure/energy, with zero velocity.

Sentinel defaults still require explicit setup parameters. Header guards add
runtime input checks beyond the CCL ranges. Library output meanings, coordinate
r semantics, and physical table validity remain external. The revised guide
uses the same Tabulated capitalization as the bodies; GID-0011 and GID-0014
record the resolved schedule-description and guide discrepancies.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [IsotropicGas.c](../../../GRHayLID/src/IsotropicGas.c)
- [documentation.tex](../../../GRHayLID/doc/documentation.tex)
- [schedule.ccl](../../../GRHayLID/schedule.ccl)

## Related Pages

- [Other Tabulated Family](constant-density-sphere.md)
- [GRHayLib Contract](../integration/grhaylib-contract.md)
- [GID-0011](../contradictions.md#gid-0011)
- [GID-0014](../contradictions.md#gid-0014)
