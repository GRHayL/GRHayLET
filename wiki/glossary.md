# IllinoisGRMHD Glossary

> Canonical recurring routing terms and one owner each. · Status: confirmed · Last reconciled: 07-17-2026

## Terms

| Term | Routing meaning | Owner |
| --- | --- | --- |
| IllinoisGRMHD | Local Cactus thorn and this KB's sole domain tree. | [IllinoisGRMHD Overview](architecture/overview.md) |
| GRMHD | General-relativistic magnetohydrodynamic evolution role documented for this thorn. | [IllinoisGRMHD Overview](architecture/overview.md) |
| Cactus thorn | Locally declared implementation unit with CCL surfaces. | [Cactus Surface and Build](architecture/cactus-surface-and-build.md) |
| CCL | Configuration, interface, parameter, and schedule declaration files. | [Cactus Surface and Build](architecture/cactus-surface-and-build.md) |
| GRHayL/GRHayLib | External library/interface boundary visible through inheritance, headers, and `ghl_*` calls. | [HydroBase, GRHayLib, and Tmunu Boundary](integration/hydrobase-grhaylib-and-tmunu.md) |
| HydroBase | Inherited primitive/vector-potential exchange boundary. | [HydroBase, GRHayLib, and Tmunu Boundary](integration/hydrobase-grhaylib-and-tmunu.md) |
| ADMBase | Inherited lapse, shift, metric, and curvature boundary. | [HydroBase, GRHayLib, and Tmunu Boundary](integration/hydrobase-grhaylib-and-tmunu.md) |
| TmunuBase | Inherited additive stress-energy destination. | [HydroBase, GRHayLib, and Tmunu Boundary](integration/hydrobase-grhaylib-and-tmunu.md) |
| MoL | Declared Method-of-Lines registration and RHS schedule boundary. | [Schedule Lifecycle](architecture/schedule-lifecycle.md) |
| primitive | Hydrodynamic state used as conversion/recovery input or output. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| conservative | Densitized evolved hydrodynamic state. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| Prim2Con | Primitive-to-conservative conversion path. | [Primitive-Conservative Conversion](evolution/primitive-conservative-conversion.md) |
| Con2Prim | Conservative-to-primitive recovery path. | [Con2Prim Recovery and Diagnostics](evolution/con2prim-recovery-and-diagnostics.md) |
| hybrid/Simple EOS | Selectors routed to Hybrid-family source handoffs. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| tabulated EOS | Selector routed to Tabulated-family source handoffs. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| entropy evolution | `evolve_entropy` branch adding entropy state and source-family handoffs. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| electron fraction | Tabulated-family `Y_e`/`Ye_star` state. | [State and EOS Modes](evolution/state-and-eos-modes.md) |
| PPM | Reconstruction path and ghost-zone constraint used locally. | [Reconstruction, Fluxes, and Sources](evolution/reconstruction-fluxes-and-sources.md) |
| HLLE/HLL | Local external-call boundaries for hydro and induction fluxes. | [Reconstruction, Fluxes, and Sources](evolution/reconstruction-fluxes-and-sources.md) |
| vector potential | Evolved semi-staggered `Ax/Ay/Az` magnetic state. | [Staggered State and Magnetic Reconstruction](magnetics/staggered-state-and-magnetic-reconstruction.md) |
| constrained transport | Repository-attributed magnetic design intent, not a new validation guarantee. | [Staggered State and Magnetic Reconstruction](magnetics/staggered-state-and-magnetic-reconstruction.md) |
| phitilde | Fully staggered densitized scalar potential used by gauge evolution. | [Staggered State and Magnetic Reconstruction](magnetics/staggered-state-and-magnetic-reconstruction.md) |
| staggered/semi-staggered grid | Declared A, scalar-potential, and B placement. | [Staggered State and Magnetic Reconstruction](magnetics/staggered-state-and-magnetic-reconstruction.md) |
| densitized B | Face-centered curl result before determinant division. | [Staggered State and Magnetic Reconstruction](magnetics/staggered-state-and-magnetic-reconstruction.md) |
| Lorenz gauge | Gauge RHS handoff using configured damping input. | [Induction and Lorenz-Gauge RHS](magnetics/induction-and-lorenz-gauge-rhs.md) |
| atmosphere reset | Recovery terminal/nonpositive-density fallback behavior. | [Con2Prim Recovery and Diagnostics](evolution/con2prim-recovery-and-diagnostics.md) |
| failure_checker | Per-point recovery repair encoding with active hundreds-marker conflict. | [Con2Prim Recovery and Diagnostics](evolution/con2prim-recovery-and-diagnostics.md) |
| test oracle | Checked-in generated comparison evidence, not current pass proof. | [Test Harness and Oracles](validation/test-harness-and-oracles.md) |
| Balsara | Five visible one-dimensional case configurations. | [Balsara and TOV Cases](validation/balsara-and-tov-cases.md) |
| magnetized TOV | Visible stellar case files and naming gap. | [Balsara and TOV Cases](validation/balsara-and-tov-cases.md) |
