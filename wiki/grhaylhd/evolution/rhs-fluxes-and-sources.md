# RHS Fluxes and Sources

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page traces declared RHS ordering and visible source/flux construction in
all four variant families plus common metric helpers. GRHayLib numerical
semantics, formal accuracy beyond local labels, stability, and conservation
remain external.

## Summary

Schedule CCL declares source RHS before flux RHS. Source routines initialize
advected-source fields, build metric, extrinsic curvature, zero-field
primitives and metric derivatives, then call a GRHayLib source helper. Flux
routines reconstruct six-point primitive stencils at faces, call directional
characteristic-speed and HLLE helpers, store fluxes, and add directional flux
differences to RHS.

## Variant Applicability

| Applicability | Zero-source fields | Reconstructed extras | Extra flux/RHS accumulation |
| --- | --- | --- | --- |
| Common | `rho_star_rhs` | Density, pressure, three velocities | Five core conservative fluxes and RHS fields |
| Hybrid/Simple | Core only | Hybrid density steepening uses local effective Gamma helper | None |
| Hybrid/Simple+Entropy | Core plus `ent_star_rhs` | Adds entropy | `ent_star_flux` to `ent_star_rhs` |
| Tabulated | Core plus `Ye_star_rhs` | Adds `Y_e`; tabulated bound and EOS helper calls | `Ye_star_flux` to `Ye_star_rhs` |
| Tabulated+Entropy | Core plus both optional RHS fields | Adds entropy and `Y_e`; tabulated helper calls | Both optional flux/RHS pairs |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-RHS-01` | Schedule declares flux evaluation after source evaluation. | declared | Flux group schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_evaluate_fluxes_rhs` |
| `EV-RHS-02` | Header labels `COMPUTE_DERIV` as fourth-order derivative. | visible-implementation | Macro and adjacent local label | `macro:GRHayLHD/src/GRHayLHD.h#name=COMPUTE_DERIV` |
| `EV-RHS-03` | Face interpolation uses four points with explicit local coefficients. | visible-implementation | Face macro | `macro:GRHayLHD/src/GRHayLHD.h#name=COMPUTE_FCVAL` |
| `EV-RHS-04` | Hybrid RHS visibly builds sources and core flux divergence. | visible-implementation | Variant source function | `c:GRHayLHD/src/Hybrid/evaluate_sources_rhs.c#symbol=GRHayLHD_hybrid_evaluate_sources_rhs` |
| `EV-RHS-05` | HybridEntropy flux visibly reconstructs and advects entropy. | visible-implementation | Variant flux function | `c:GRHayLHD/src/HybridEntropy/evaluate_fluxes_rhs.c#symbol=GRHayLHD_hybrid_entropy_evaluate_fluxes_rhs` |
| `EV-RHS-06` | Tabulated flux visibly reconstructs and advects electron fraction. | visible-implementation | Variant flux function | `c:GRHayLHD/src/Tabulated/evaluate_fluxes_rhs.c#symbol=GRHayLHD_tabulated_evaluate_fluxes_rhs` |
| `EV-RHS-07` | TabulatedEntropy flux visibly reconstructs and advects both extras. | visible-implementation | Variant flux function | `c:GRHayLHD/src/TabulatedEntropy/evaluate_fluxes_rhs.c#symbol=GRHayLHD_tabulated_entropy_evaluate_fluxes_rhs` |
| `EV-RHS-08` | Metric-derivative helper visibly samples four offsets and applies local derivative macro. | visible-implementation | Common derivative function | `c:GRHayLHD/src/compute_metric_derivs.c#symbol=GRHayLHD_compute_metric_derivs` |
| `EV-RHS-09` | Face helper visibly applies local interpolation macro to four input points. | visible-implementation | Common face function | `c:GRHayLHD/src/interpolate_metric_to_face.c#symbol=GRHayLHD_interpolate_metric_to_face` |
| `EV-RHS-10` | Setup routine visibly calls an error API when any ghost-zone count is below three, with a message that local PPM requires three ghost zones. | visible-implementation | Setup precondition check | `c:GRHayLHD/src/InitSymBound.c#symbol=GRHayLHD_InitSymBound` |

## Details

### Declared ordering

`GRHayLHD_RHS` is scheduled in `MoL_CalcRHS`. Within it,
`GRHayLHD_evaluate_sources_rhs` has no local `after` clause and
`GRHayLHD_evaluate_fluxes_rhs` is explicitly scheduled after source group.
Thus source-before-flux is declared ordering only.

### Ghost-zone precondition

`GRHayLHD_InitSymBound` visibly calls `CCTK_VERROR` when any
`cctk_nghostzones` component is below three; its message attributes the
requirement to this thorn's PPM implementation. This is a visible local
precondition check, not proof that setup ran or rejected a configuration.

### Source construction

Each source function loops over interior bounds defined by ghost-zone counts.
It sets density RHS to zero and also sets active entropy/electron-fraction RHS
to zero. It builds ADM metric and extrinsic curvature from local gridfunctions,
loads mode primitives with `BU[0..2] = 0`, limits velocity and computes `u0`,
then obtains metric derivatives in all three directions through
`GRHayLHD_compute_metric_derivs`. Finally it calls
`ghl_calculate_source_terms` and assigns returned tau and momentum sources.

`GRHayLHD_compute_metric_derivs` samples offsets minus two, minus one, plus one,
plus two and applies `COMPUTE_DERIV` times inverse spacing to lapse, shift, and
six spatial-metric components. Header comment explicitly labels this macro
fourth-order. Claim is limited to local label and formula; no measured order is
asserted.

### Face metric and primitive reconstruction

`GRHayLHD_interpolate_metric_to_face` applies `COMPUTE_FCVAL` to four values at
offsets minus two, minus one, zero, plus one with coefficients `-0.0625`,
`0.5625`, `0.5625`, `-0.0625`. This page calls it four-point interpolation and
makes no unsupported order claim.

For each of three flux directions, variant code collects six values at offsets
minus three through plus two for density, pressure, velocities, and active
extras. It calls `ghl_compute_ftilde`, density PPM reconstruction with
steepening, ordinary PPM reconstruction for pressure/velocity/extras, and
zeros left/right `BU`. Hybrid families obtain effective Gamma from a local
helper; tabulated families pass `1.0` to density steepening and visibly invoke
tabulated bound and pressure-to-energy/temperature helpers.

Directional function pointers select characteristic-speed and mode-specific
HLLE helpers for axes 0, 1, 2. Returned core and optional conservative fluxes
are stored. Interior RHS then accumulates, for each direction,
`inverse_spacing * (flux[index] - flux[index+direction])` for every active
field.

### Independent inspection notes

All four source files construct matching mode primitive fields and zero matching
advected-source RHS fields. All four flux files share core stencil shape but
use distinct mode-specific HLLE symbols and optional reconstruction/flux
fields. Similar structure is not treated as identity.

## Caveats

- PPM, steepening, characteristic-speed, HLLE, bounds, and EOS helper semantics
  are delegated to GRHayLib.
- Source-before-flux is a CCL declaration, not observed scheduler execution.
- Local comments still use some MHD wording; zero `BU` and hydrodynamic state
  do not create magnetic evolution.
- No build, run, convergence, stability, or conservation claim is made.

## Sources

- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Setup precondition](../../../GRHayLHD/src/InitSymBound.c)
- [Common header macros](../../../GRHayLHD/src/GRHayLHD.h)
- [Metric derivative helper](../../../GRHayLHD/src/compute_metric_derivs.c)
- [Face metric interpolation](../../../GRHayLHD/src/interpolate_metric_to_face.c)
- Four source and four flux implementations under [variant source tree](../../../GRHayLHD/src)

## Related Pages

- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
