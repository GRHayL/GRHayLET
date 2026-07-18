# Schedule Lifecycle

> Status: confirmed · Last reconciled: 07-18-2026
> Up: [Architecture](index.md)

## Summary

`schedule.ccl` declares registration, initial conversion, recurring recovery,
optional stress-energy coupling, MoL RHS, analysis conversion, four conditional
variant handoffs, and a backward-compatible branch. `MoL_registration.c`
declares evolved/RHS, constrained, and save/restore registrations. These are
schedule and call-shape claims only; no runtime execution was observed.

## Detail

### Storage and Registration

Unconditional storage covers HydroBase magnetic/vector-potential groups,
three-timelevel hydrodynamic conservatives and EM potential state, primitives,
centered/staggered B, RHS, reconstructed/characteristic/flux temporaries,
`Ye_star_flux`, `ent_star_flux`, and `failure_checker`. Tabulated selection adds
three-timelevel `Ye_star` and `Ye_star_rhs` while repeating `Ye_star_flux` in
its conditional storage clause; entropy evolution likewise adds `ent_star` and
`ent_star_rhs` while repeating `ent_star_flux`.

Setup declares:

1. `IllinoisGRMHD_RegisterVars` in `MoL_Register`;
2. Driver BC registration in `Driver_BoundarySelect`;
3. `IllinoisGRMHD_InitSymBound` at `BASEGRID`.

`IllinoisGRMHD_RegisterVars` registers `Ax`, `Ay`, `Az`, and `phitilde` against
their scalar RHS variables, and registers conservative groups against their RHS
groups. Entropy and electron-fraction groups are conditional. HydroBase rho,
pressure, and eps are constrained; entropy and Y_e/temperature are conditional;
Tmunu groups are constrained when enabled, but the three Tmunu registration
return values are not accumulated into `ierr`. Failures returned by those calls
therefore do not reach the later `CCTK_ERROR` checks through `ierr`. ADMBase
lapse, shift, metric, and curvature are save-and-restore groups. Every other
registration call accumulates its return value into `ierr` before a check.

Claim evidence:
- Claim: Local code requests these MoL registration classes conditionally; the three enabled Tmunu calls do not accumulate their return values into `ierr`, so the local error checks do not test those returned statuses; this does not establish MoL internals or successful registration at runtime.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/src/MoL_registration.c`, `IllinoisGRMHD_RegisterVars`
- Corroboration: registered `IllinoisGRMHD/interface.ccl`, MoL alias declarations; `IllinoisGRMHD/schedule.ccl`, `MoL_Register` entry
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=three Tmunu return statuses inspected-not-run; options=evolve_entropy,tabulated EOS,update_Tmunu branches; date=07-18-2026`

### Declared Phase Graph

1. **Initial data:** `IllinoisGRMHD_Prim2Con2Prim` runs in
   `HydroBase_Prim2ConInitial`: HydroBase ingress; optional primitive
   perturbation and A-to-B reconstruction are each declared after ingress,
   but `schedule.ccl` declares no relative order between them. Both are
   declared before Prim2Con, followed by Con2Prim and optional HydroBase
   egress.
2. **Post-initial symmetry:** `IllinoisGRMHD_set_gz_symmetries` is declared at
   `CCTK_POSTPOSTINITIAL` after `Con2Prim`.
3. **Recurring recovery:** `IllinoisGRMHD_Con2Prim` runs in
   `HydroBase_Con2Prim`: declared sync, A outer boundary, B reconstruction,
   optional conservative perturbation, Con2Prim, then matter outer boundary.
4. **Optional coupling:** `IllinoisGRMHD_compute_Tmunu` is in `AddToTmunu` when
   `update_Tmunu` is true.
5. **RHS:** `IllinoisGRMHD_RHS` is in `MoL_CalcRHS`; sources precede fluxes;
   an active `NRPyLeakageET` condition adds a HydroBase velocity conversion;
   EM gauge RHS is after flux evaluation.
6. **Analysis:** positive `Convert_to_HydroBase_every` declares egress at
   `CCTK_ANALYSIS` with named before/after constraints.
7. **Variants and compatibility:** EOS/entropy conditions schedule one family
   into each Prim2Con, Con2Prim, hydro-boundary, source, flux, and perturbation
   handoff. Exact selection and state extras belong to
   [State and EOS Modes](../evolution/state-and-eos-modes.md). An active
   `ID_converter_ILGRMHD` condition declares compatibility storage,
   initialization/data copies, forced conversions, and Hybrid-family handoffs.

The compatibility branch is subject to
[`CONTR-0001`](../contradictions.md#contr-0001); current-tree scheduling does
not establish external release policy.

### Sync Handoffs

The initial B reconstruction declares synchronization of centered and staggered
B. Recurring `IllinoisGRMHD_sync` unconditionally lists conservatives,
`Ye_star`, `ent_star` (groups whose storage is conditional), and EM potential
state in its `SYNC` clause before A boundaries. Its local C function body is
empty; synchronization is expressed by the CCL `SYNC` clause. Variant Prim2Con
and hydro-boundary entries add mode-specific `SYNC` lists. These are declared
handoffs, not observed communication.

Claim evidence:
- Claim: `IllinoisGRMHD_sync` has an empty C body while its schedule entry carries the conservative/EM `SYNC` declaration; static inspection does not prove synchronization execution.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/schedule.ccl`, `IllinoisGRMHD_sync` schedule entry and `SYNC` clause
- Corroboration: registered `IllinoisGRMHD/src/sync.c`, `IllinoisGRMHD_sync`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-applicable; options=declared sync list; date=07-17-2026`

### Ownership Boundaries

This page owns phase/order and handoff declarations. It does not restate exact
2x2 EOS/entropy dispatch, conversion kernels, recovery ladders, matter-boundary
return rules, magnetic algorithms, or external scheduler semantics. Schedule
verbs here mean “declares” or “schedules,” never “was observed to execute.”

## Sources

- [`schedule.ccl`](../../IllinoisGRMHD/schedule.ccl) — storage conditions,
  groups, functions, order constraints, conditional families, and compatibility
  branch.
- [`interface.ccl`](../../IllinoisGRMHD/interface.ccl) — declared state and MoL/
  Driver alias surface used by scheduled routines.
- [`MoL_registration.c`](../../IllinoisGRMHD/src/MoL_registration.c) —
  `IllinoisGRMHD_RegisterVars`.
- [`sync.c`](../../IllinoisGRMHD/src/sync.c) — empty `IllinoisGRMHD_sync` body.

## See Also

- Parent: [Architecture](index.md)
- Depends on: [State and EOS Modes](../evolution/state-and-eos-modes.md)
- See also: [Cactus Surface and Build](cactus-surface-and-build.md)
- Implements: [HydroBase, GRHayLib, and Tmunu](../integration/hydrobase-grhaylib-and-tmunu.md)
- See also: [Matter Boundaries and Perturbations](../evolution/matter-boundaries-and-perturbations.md)
