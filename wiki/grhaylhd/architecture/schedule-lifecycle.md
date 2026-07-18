# Declared Schedule Lifecycle

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page reconstructs schedule declarations, guards, order constraints, and
declared read/write/sync sets from local CCL. It describes declared lifecycle,
not observed execution, framework scheduling semantics, or successful MoL
integration.

## Summary

GRHayLHD declares setup at `MoL_Register` and `BASEGRID`, initial conversion in
`HydroBase_Prim2ConInitial`, recovery and primitive boundaries in
`HydroBase_Con2Prim`, optional stress-energy contribution in `AddToTmunu`,
source then flux RHS groups in `MoL_CalcRHS`, and optional HydroBase output at
`CCTK_ANALYSIS`. Four EOS/entropy branches populate seven operation families.

## Variant Applicability

| Applicability | Scheduled concrete prefix | Conditional state in declared reads/writes/syncs |
| --- | --- | --- |
| Common | Setup, converters, Tmunu, and abstract operation groups | Core conservative, velocity, metric, and base HydroBase fields |
| Hybrid/Simple | `GRHayLHD_hybrid_` | Core fields only |
| Hybrid/Simple+Entropy | `GRHayLHD_hybrid_entropy_` | Adds entropy primitive, conservative, RHS, and flux fields |
| Tabulated | `GRHayLHD_tabulated_` | Adds electron fraction, temperature, `Ye_star`, RHS, and flux |
| Tabulated+Entropy | `GRHayLHD_tabulated_entropy_` | Adds both entropy and electron-fraction state |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ARCH-LIFE-01` | Schedule declares MoL registration. | declared | `MoL_Register` schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_RegisterVars` |
| `ARCH-LIFE-02` | Schedule declares symmetry initialization at `BASEGRID`. | declared | Setup schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_InitSymBound` |
| `ARCH-LIFE-03` | Schedule declares initial conversion group. | declared | Prim2Con2Prim group | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_Prim2Con2Prim` |
| `ARCH-LIFE-04` | Schedule declares evolution recovery group. | declared | Con2Prim group | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_Con2Prim` |
| `ARCH-LIFE-05` | Schedule declares optional stress-energy computation. | declared | Tmunu schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_compute_Tmunu` |
| `ARCH-LIFE-06` | Schedule declares RHS group. | declared | MoL RHS schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_RHS` |
| `ARCH-LIFE-07` | Leakage-active branch schedules a HydroBase converter without a positive cadence guard. | declared | Leakage schedule branch | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=GRHayLHD_RHS` |

## Details

### Phase graph

1. `MoL_Register` schedules `GRHayLHD_RegisterVars` with `META`; `BASEGRID`
   schedules `GRHayLHD_InitSymBound`.
2. `HydroBase_Prim2ConInitial` declares `GRHayLHD_Prim2Con2Prim`:
   `convert_HydroBase_to_GRHayLHD`; optional primitive perturbation after that
   converter and before Prim2Con; Prim2Con after the converter; Con2Prim after
   Prim2Con; optional reverse conversion after Con2Prim when
   `Convert_to_HydroBase_every` is nonzero.
3. `HydroBase_Con2Prim` declares `GRHayLHD_Con2Prim`: empty sync schedule for
   core, `Ye_star`, and `ent_star`; optional conservative perturbation after
   sync and before recovery; recovery after sync; primitive outer boundaries
   after recovery.
4. `AddToTmunu` schedules `GRHayLHD_compute_Tmunu` only under
   `update_Tmunu`.
5. `MoL_CalcRHS` declares `GRHayLHD_RHS`: source group first, flux group
   explicitly `after` sources, then an NRPyLeakageET-active conversion after
   fluxes.
6. `CCTK_ANALYSIS` optionally schedules reverse conversion under nonzero
   `Convert_to_HydroBase_every`, before named diagnostic groups and after
   `ML_BSSN_evolCalcGroup`, with `GLOBAL-EARLY,LOOP-LOCAL` options.

### Declared data movement

Initial HydroBase-to-native conversion reads lapse, shift, and
`HydroBase::vel`; writes and syncs native velocities. Prim2Con variants read
metric/gauge, HydroBase mode primitives, and native velocities; write `u0`,
core conservatives, normalized primitives, and optional evolved state; their
declared sync sets follow the same optional fields. Recovery variants read
coordinates, metric/gauge, core and optional conservatives; write corrected
conservatives, primitives, `u0`, and `failure_checker`.

Outer-boundary variants read metric/gauge and their mode primitives; write
limited primitives, `u0`, core conservatives, and optional state. Their sync
sets cover primitive fields, not the newly recomputed conservative groups.
Source variants read metric/gauge/curvature and primitives and write mode RHS
groups. Flux variants read metric/gauge, primitives, and interior RHS; write
flux temporaries and accumulated RHS fields.

Tmunu reads metric/gauge, base thermodynamics, native velocity, and `u0`, then
declares writes to all three TmunuBase stress-energy groups. Reverse conversion
at initial and analysis phases declares both HydroBase velocity and Lorentz
factor writes; leakage conversion declares only HydroBase velocity writes.

### Four-way dispatch

Nested conditions choose Hybrid or Simple versus Tabulated, then entropy false
versus true. Each branch schedules exactly one concrete function into each of
seven abstract operation groups: Prim2Con, Con2Prim, outer boundary, source
RHS, flux RHS, primitive perturbation, and conservative perturbation. Exact
field differences are owned by the canonical EOS and entropy variant matrix.

## Caveats

- This page uses CCL verbs: declares, schedules, reads, writes, and syncs. It
  makes no observed-execution claim.
- Generic Cactus schedule-bin ordering and MoL time-integration semantics are
  external and were not needed to inventory local declarations.
- Leakage scheduling is gated only by thorn activity, while converter code
  applies modulo `Convert_to_HydroBase_every`; see
  [GRH-0006](../contradictions.md#grh-0006).
- Leakage declaration omits body-visible metric reads and Lorentz-factor
  write; see [GRH-0011](../contradictions.md#grh-0011).
- `update_Tmunu` steering versus setup-conditional registration is documented
  by Integration owner page; this page states only local schedule condition.

## Sources

- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Parameters](../../../GRHayLHD/param.ccl)
- [HydroBase reverse converter](../../../GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c)

## Related Pages

- [Purpose and Build Surface](purpose-build-surface.md)
- [Variables and Storage](variables-and-storage.md)
