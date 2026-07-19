# Beta Equilibrium

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns local declarations and visible dataflow for imposing neutrino-
free beta equilibrium on pre-existing initial data. It does not establish EOS
table contents, GRHayLib algorithms, error termination, successful scheduling,
or numerical validity.

## Summary

README describes beta equilibrium as usable with any initial-data thorn when
the surrounding GRHayLib and schedule setup is suitable. Schedule CCL places
`GRHayLID_BetaEquilibrium` at `CCTK_INITIAL` after `HydroBase_Initial` and
before `HydroBase_Prim2ConInitial` when `impose_beta_equilibrium` is true. The
body checks for a tabulated EOS, visibly compares `beq_temperature` with table
bounds before applying `CHECK_PARAMETER`, initializes a GRHayLib beta-
equilibrium helper, and then writes pressure, specific internal energy,
electron fraction, and temperature. It does not write entropy.

## Mode Applicability

| Applicability | Visible behavior |
| --- | --- |
| BetaEquilibrium | Boolean `impose_beta_equilibrium` guards the schedule; the body requires the local `ghl_eos` handle to report tabulated EOS. |
| Common | README declares that the feature need not use GRHayLID as the initial-data source; local compatibility and external initialization remain prerequisites. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-BEQ-01` | `impose_beta_equilibrium` is a boolean switch defaulting to `no`. | declared | Local parameter declaration | `ccl:GRHayLID/param.ccl#parameter=impose_beta_equilibrium` |
| `ID-BEQ-02` | Schedule CCL declares the routine at `CCTK_INITIAL`, ordered after `HydroBase_Initial` and before `HydroBase_Prim2ConInitial`, with writes to pressure, specific internal energy, electron fraction, and temperature. | declared | Beta-equilibrium schedule block | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |
| `ID-BEQ-03` | The function visibly compares `beq_temperature` with handle fields before invoking `CHECK_PARAMETER`, then branches at `1.01*ghl_eos->rho_atm` and writes four HydroBase quantities. | visible-implementation | Complete local function body | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` |
| `ID-BEQ-04` | The setup call's computation and resulting handle state are delegated to GRHayLib and unverified locally. | out-of-scope | Visible external call boundary | `c:GRHayLID/src/BetaEquilibrium.c#call=ghl_tabulated_compute_Ye_of_rho_beq_constant_T?function=GRHayLID_BetaEquilibrium` |
| `ID-BEQ-05` | README declares beta equilibrium usable with initial data supplied by another thorn, subject to stated GRHayLib and schedule conditions. | declared | Purpose statement | `doc:GRHayLID/README#section=1. Purpose` |
| `ID-BEQ-06` | `beq_temperature` is a nonnegative real whose forbidden default is `-1`. | declared | Temperature declaration | `ccl:GRHayLID/param.ccl#parameter=beq_temperature` |

## Details

### Declared lifecycle

`impose_beta_equilibrium` defaults to `no`. When selected, schedule CCL reads
`HydroBase::rho` and declares writes to `press`, `eps`, `Y_e`, and
`temperature`. Its description additionally names entropy, but neither its
`WRITES` clause nor the C body writes entropy. This mismatch is tracked as
[GID-0008](../contradictions.md#gid-0008).

### Visible pre-loop ordering

The body first compares `ghl_eos->eos_type` with a tabulated enum. It then
compares `beq_temperature` with `ghl_eos->table_T_min` and
`ghl_eos->table_T_max`; only after that comparison does it invoke
`CHECK_PARAMETER(beq_temperature)`. Thus the declared `-1` sentinel visibly
reaches the bounds comparison first. Whether either error macro terminates and
what table bounds contain are external, so no reachability or outcome follows.
See [GID-0007](../contradictions.md#gid-0007).

### Visible pointwise branches

Before iteration, the function calls
`ghl_tabulated_compute_Ye_of_rho_beq_constant_T` and passes its return code to
`ghl_abort_if_error`. At points where `rho <= 1.01*ghl_eos->rho_atm`, it copies
the handle's atmosphere pressure, energy, electron fraction, and temperature.
Elsewhere it calls helpers to obtain electron fraction and pressure/energy,
checks returned codes locally, and writes those outputs with
`beq_temperature`. These are visible reads, calls, comparisons, and writes—not
claims about the called routines' results.

## Caveats

- The shared `ghl_eos` handle and every `ghl_*` result are external boundaries.
- CCL order establishes a declared lifecycle, not actual schedule execution.
- Local source does not prove error-macro termination or that the sentinel is
  accepted or rejected by external table bounds.
- No checked-in test or oracle demonstrates this path.

## Sources

- [Purpose statement](../../../GRHayLID/README)
- [Parameter declarations](../../../GRHayLID/param.ccl)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)
- [Beta-equilibrium implementation](../../../GRHayLID/src/BetaEquilibrium.c)

## Related Pages

- [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [Entropy Computation](entropy-computation.md)
- [Parameters and Configurations](../integration/parameters-and-configurations.md)
- [GID-0007: temperature check ordering](../contradictions.md#gid-0007)
- [GID-0008: entropy write mismatch](../contradictions.md#gid-0008)
