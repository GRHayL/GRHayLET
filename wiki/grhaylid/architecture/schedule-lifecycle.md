# Declared Schedule Lifecycle

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page maps all seven local schedule declarations, their textual guards,
bins, groups, ordering clauses, aliases, and declared read/write sets. It does
not establish observed scheduling, bin or alias semantics, external group
existence, error behavior, or successful execution.

## Summary

Three `initial_hydro` choices declare hydrodynamic setup in
`HydroBase_Initial`. HydroTest1D can additionally declare magnetic setup there
under `initialize_magnetic_quantities`. Beta equilibrium is declared at
`CCTK_INITIAL` between `HydroBase_Initial` and
`HydroBase_Prim2ConInitial`, with alias `impose_beta_equilibrium`. Entropy
declarations select Hybrid or Tabulated routines using `initial_entropy` and
`EOS_type`, ordered after both `HydroBase_Initial` and that alias and before
`HydroBase_Prim2ConInitial`.

## Mode Applicability

| Applicability | Declared schedule selection |
| --- | --- |
| Common | `initial_hydro`, standalone beta-equilibrium, and standalone entropy guards form the common dispatch surface. |
| HydroTest1D | Declares hydrodynamic setup in `HydroBase_Initial`. |
| HydroTest1D+Magnetic | Adds magnetic setup after the one-dimensional hydrodynamic routine. |
| IsotropicGas | Declares IsotropicGas setup in `HydroBase_Initial`. |
| ConstantDensitySphere | Declares ConstantDensitySphere setup in `HydroBase_Initial`. |
| BetaEquilibrium | Declares beta-equilibrium update at `CCTK_INITIAL`. |
| Entropy/Hybrid | Declares Hybrid entropy computation at `CCTK_INITIAL`. |
| Entropy/Tabulated | Declares Tabulated entropy computation at `CCTK_INITIAL`. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ARCH-LIFE-01` | `initial_hydro="HydroTest1D"` declares hydrodynamic setup in `HydroBase_Initial`. | declared | Hydrodynamic schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_hydro_data` |
| `ARCH-LIFE-02` | Nested magnetic guard declares magnetic setup in `HydroBase_Initial` after the hydrodynamic routine. | declared | Magnetic schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_magnetic_data` |
| `ARCH-LIFE-03` | `initial_hydro="IsotropicGas"` declares IsotropicGas setup in `HydroBase_Initial`. | declared | IsotropicGas schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_IsotropicGas` |
| `ARCH-LIFE-04` | `initial_hydro="ConstantDensitySphere"` declares sphere setup in `HydroBase_Initial`. | declared | Sphere schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_ConstantDensitySphere` |
| `ARCH-LIFE-05` | Boolean guard declares beta equilibrium at `CCTK_INITIAL` with ordering and alias clauses. | declared | Beta-equilibrium schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |
| `ARCH-LIFE-06` | Entropy guard plus `EOS_type="Hybrid"` declares Hybrid entropy computation. | declared | Hybrid entropy schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` |
| `ARCH-LIFE-07` | Entropy guard plus `EOS_type="Tabulated"` declares Tabulated entropy computation. | declared | Tabulated entropy schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` |
| `ARCH-LIFE-08` | `initialize_magnetic_quantities` is a local boolean defaulting to `yes` and guards magnetic scheduling. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=initialize_magnetic_quantities` |

## Details

### HydroBase_Initial dispatch

The first conditional chain compares `initial_hydro` with three keyword
values:

1. `HydroTest1D` declares `GRHayLID_1D_tests_hydro_data` in
   `HydroBase_Initial`. It reads `Grid::coordinates` and writes HydroBase
   `rho`, `press`, `eps`, and `vel` everywhere.
2. Inside that arm, `initialize_magnetic_quantities` guards
   `GRHayLID_1D_tests_magnetic_data` in the same group, textually ordered
   `after GRHayLID_1D_tests_hydro_data`. It reads `Grid::coordinates` and
   writes HydroBase `Avec` and `Bvec` everywhere.
3. `IsotropicGas` declares `GRHayLID_IsotropicGas` in `HydroBase_Initial`.
   It reads `Grid::coordinates` and writes HydroBase `rho`, `press`, `eps`,
   `vel`, `Y_e`, and `temperature` everywhere.
4. `ConstantDensitySphere` declares `GRHayLID_ConstantDensitySphere` in the
   same group with the same declared read/write set as IsotropicGas.

These CCL blocks say `schedule`; this page makes no claim that a condition was
selected or a routine ran.

### Standalone beta-equilibrium declaration

Under `impose_beta_equilibrium`, `GRHayLID_BetaEquilibrium` is declared
`at CCTK_INITIAL after HydroBase_Initial before
HydroBase_Prim2ConInitial as impose_beta_equilibrium`. It reads
`HydroBase::rho` and writes `press`, `eps`, `Y_e`, and `temperature`
everywhere. Its description also names entropy, although entropy is absent
from both declared `WRITES` and the visible body; GID-0008 records that
mismatch.

### Entropy declaration and alias-name ordering

An outer comparison requires `initial_entropy="GRHayLID"`. Its nested chain
declares exactly two EOS arms:

- `EOS_type="Hybrid"` declares `GRHayLID_compute_entropy_hybrid`, reading
  `rho` and `press` and writing `entropy` everywhere.
- `EOS_type="Tabulated"` declares
  `GRHayLID_compute_entropy_tabulated`, reading `rho`, `Y_e`, and
  `temperature` and writing `rho`, `press`, `eps`, `entropy`, `Y_e`, and
  `temperature` everywhere.

Both are declared at `CCTK_INITIAL after (HydroBase_Initial
impose_beta_equilibrium) before HydroBase_Prim2ConInitial`. The token
`impose_beta_equilibrium` is the alias introduced by the beta-equilibrium
schedule block, not its C routine name. Local text establishes this naming
relationship but not framework interpretation. No third EOS arm or final
error arm is declared.

## Caveats

- CCL ordering is a declaration, not evidence of execution or external
  scheduler behavior.
- Existence and semantics of `HydroBase_Initial`, `CCTK_INITIAL`, and
  `HydroBase_Prim2ConInitial` remain external.
- Beta-equilibrium description includes entropy while its declared writes do
  not; see [GID-0008](../contradictions.md#gid-0008).
- Entropy dispatch has no arm outside Hybrid and Tabulated and no error arm;
  see [GID-0010](../contradictions.md#gid-0010).
- Both three-dimensional setup descriptions say "1D test"; see
  [GID-0011](../contradictions.md#gid-0011).

## Sources

- [Schedule declarations](../../../GRHayLID/schedule.ccl)
- [Parameter declarations](../../../GRHayLID/param.ccl)

## Related Pages

- [Purpose and Build Surface](purpose-and-build-surface.md)
- [Initial Data](../initial-data/index.md)
- [HydroBase Keyword Extensions](../integration/hydrobase-keyword-extensions.md)
- [GID-0008](../contradictions.md#gid-0008)
- [GID-0010](../contradictions.md#gid-0010)
- [GID-0011](../contradictions.md#gid-0011)
