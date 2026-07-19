# GRHayLib Contract

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page inventories the GRHayLib shared selection, header boundary, handle
fields, and API-call families visibly used by GRHayLID. It does not establish
GRHayLib initialization, ABI compatibility, handle contents, EOS-table
discovery, algorithms, error semantics, or numerical results.

## Summary

Parameter CCL shares GRHayLib and uses its `EOS_type` keyword. The local header
includes `GRHayLib.h`; interface CCL also declares that include. The claim
table anchors representative handle reads and `ghl_*` calls in beta-equilibrium
and tabulated-entropy functions; no complete per-unit usage count is asserted.
Local build configuration declares HDF5 as a requirement. All called-library
behavior remains external.

## Mode Applicability

| Applicability | Visible GRHayLib boundary |
| --- | --- |
| Common | Shared `EOS_type`, `GRHayLib.h`, `ghl_eos`, error-code names, and declared HDF5 requirement. |
| HydroTest1D | Reads `Gamma_ppoly` and calls a Hybrid polytropic-index helper. |
| IsotropicGas | Calls a Tabulated pressure/energy-from-temperature helper once before grid writes. |
| ConstantDensitySphere | Calls that helper for interior and exterior state. |
| BetaEquilibrium | Reads EOS kind, table bounds, atmosphere fields, and calls beta-equilibrium/Tabulated/error helpers. |
| Entropy/Hybrid | Calls a Hybrid entropy helper. |
| Entropy/Tabulated | Calls Tabulated bounds and pressure/energy/entropy helpers. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-GHL-01` | Parameter CCL shares GRHayLib and uses the shared `EOS_type` keyword without declaring its domain or default locally. | declared | Shared-keyword declaration | `ccl:GRHayLID/param.ccl#parameter=EOS_type` |
| `INT-GHL-02` | Local header visibly includes `GRHayLib.h`. | visible-implementation | Include directive | `macro:GRHayLID/src/GRHayLID.h#include=GRHayLib.h` |
| `INT-GHL-03` | Build configuration declares HDF5 as a requirement. | declared | Local configuration input | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` |
| `INT-GHL-04` | Beta-equilibrium implementation visibly reads `ghl_eos` fields and calls several `ghl_*` APIs with local error-code checks. | visible-implementation | Beta-equilibrium function | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` |
| `INT-GHL-05` | Semantics of the BetaEquilibrium calls, handle fields, and error helpers are delegated and unverified locally. | out-of-scope | External boundary at local function | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` |
| `INT-GHL-06` | Tabulated entropy implementation visibly invokes bounds and pressure/energy/entropy call families. | visible-implementation | Tabulated entropy function | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` |
| `INT-GHL-07` | Semantics and mutation guarantees of those Tabulated calls are delegated and unverified locally. | out-of-scope | External boundary at local function | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` |

## Details

### Selection and include boundary

`param.ccl` uses rather than redeclares `EOS_type`. Local schedules compare it
with `Hybrid` and `Tabulated`; setup bodies also compare it with `Tabulated`.
The domain, default, normalization, and relation to `ghl_eos->eos_type` are not
owned locally. `src/GRHayLID.h` includes `GRHayLib.h`, while interface CCL
declares `USES INCLUDE: GRHayLib.h`; neither declaration proves include
discovery or successful compilation.

### Visible handle fields

Local code reads these fields and no semantics beyond the reads are inferred:

- `Gamma_ppoly` in one-dimensional hydro setup;
- `eos_type`, `table_T_min`, and `table_T_max` in beta equilibrium; and
- `rho_atm`, `press_atm`, `eps_atm`, `Y_e_atm`, and `T_atm` in its atmosphere
  branch.

### Visible call families

| Family | Local names |
| --- | --- |
| Hybrid lookup/entropy | `ghl_hybrid_find_polytropic_index`, `ghl_hybrid_compute_entropy_function` |
| Tabulated pressure/energy/entropy | `ghl_tabulated_compute_P_eps_from_T`, `ghl_tabulated_compute_P_eps_S_from_T` |
| Tabulated bounds and beta equilibrium | `ghl_tabulated_enforce_bounds_rho_Ye_T`, `ghl_tabulated_compute_Ye_of_rho_beq_constant_T`, `ghl_tabulated_compute_Ye_from_rho` |
| Errors | `ghl_abort_if_error`, `ghl_error_codes_t`, `ghl_success` |

Names prove only visible local calls or comparisons. Return contracts,
mutation rules, supported values, table contents, formulas, and accuracy are
not present in this source tree.

### HDF5 boundary

`configuration.ccl` consists of `requires HDF5` and lacks a trailing newline.
The declaration is still locatable. It proves a checked-in build requirement,
not HDF5 discovery, linkage, table availability, or a successful build.

## Caveats

- C symbol/field names do not establish library implementation or ABI.
- `CCTK_ERROR`, `CCTK_VERROR`, and `ghl_abort_if_error` termination semantics
  are external.
- A build requirement does not establish dependency availability.
- No external GRHayLib or HDF5 source is admitted as GRHayLID domain evidence.

## Sources

- [Shared parameter declaration](../../../GRHayLID/param.ccl)
- [Local common header](../../../GRHayLID/src/GRHayLID.h)
- [Build requirement](../../../GRHayLID/configuration.ccl)
- [One-dimensional hydro calls](../../../GRHayLID/src/1D_tests_hydro_data.c)
- [IsotropicGas calls](../../../GRHayLID/src/IsotropicGas.c)
- [ConstantDensitySphere calls](../../../GRHayLID/src/ConstantDensitySphere.c)
- [Beta-equilibrium calls](../../../GRHayLID/src/BetaEquilibrium.c)
- [Entropy calls](../../../GRHayLID/src/ComputeEntropy.c)

## Related Pages

- [HydroBase Keyword Extensions](hydrobase-keyword-extensions.md)
- [Parameters and Configurations](parameters-and-configurations.md)
- [Beta Equilibrium](../initial-data/beta-equilibrium.md)
- [Entropy Computation](../initial-data/entropy-computation.md)
