# GRHayLHD Glossary

> Page status: reviewed · Last reviewed: 07-17-2026

Each recurring term has one canonical owner. Meanings route queries; they do
not add claims beyond owner-page evidence.

| Term | Routing meaning | Owner |
| --- | --- | --- |
| GRHayLHD | Local hydrodynamics thorn and sole domain tree for this branch. | [Purpose and Build Surface](architecture/purpose-build-surface.md) |
| GRHD | General-relativistic hydrodynamics purpose stated for GRHayLHD. | [Purpose and Build Surface](architecture/purpose-build-surface.md) |
| non-magnetic | Documented omission of magnetic evolution and visible zero-`BU` setup. | [Purpose and Build Surface](architecture/purpose-build-surface.md) |
| CCL | Local interface, parameter, schedule, and configuration declarations. | [Declared Schedule Lifecycle](architecture/schedule-lifecycle.md) |
| storage | Declared core and mode-conditional gridfunction allocation surface. | [Variables and Storage](architecture/variables-and-storage.md) |
| Hybrid/Simple | Shared local source family selected by Hybrid or Simple `EOS_type`. | [EOS and Entropy Variants](evolution/eos-entropy-variants.md) |
| Tabulated | Local source family adding electron-fraction and temperature state. | [EOS and Entropy Variants](evolution/eos-entropy-variants.md) |
| entropy evolution | `evolve_entropy` branch adding entropy state and operation family. | [EOS and Entropy Variants](evolution/eos-entropy-variants.md) |
| primitive | Hydrodynamic state used by conversion, recovery, flux, boundary, and Tmunu paths. | [Primitive-Conservative Conversion](evolution/primitive-conservative-conversion.md) |
| conservative | Densitized evolved hydrodynamic state and optional entropy/electron-fraction state. | [Variables and Storage](architecture/variables-and-storage.md) |
| Prim2Con | Primitive-to-conservative construction path. | [Primitive-Conservative Conversion](evolution/primitive-conservative-conversion.md) |
| Con2Prim | Conservative-to-primitive recovery path. | [Conservative Recovery](evolution/conservative-recovery.md) |
| atmosphere | Non-positive-density or terminal-recovery fallback primitive state. | [Conservative Recovery](evolution/conservative-recovery.md) |
| Font1D | Explicit Hybrid-family post-retry recovery call; external semantics remain delegated. | [Conservative Recovery](evolution/conservative-recovery.md) |
| failure_checker | Per-point recovery diagnostic with open legend/writeback mismatch. | [Conservative Recovery](evolution/conservative-recovery.md) |
| PPM | External reconstruction call family used by local flux routines. | [RHS Fluxes and Sources](evolution/rhs-fluxes-and-sources.md) |
| HLLE | External directional flux call family used by local flux routines. | [RHS Fluxes and Sources](evolution/rhs-fluxes-and-sources.md) |
| Matter_BC | Local selector for copy, outflow, or frozen matter boundaries. | [Matter Boundaries and Symmetry](evolution/matter-boundaries-and-symmetry.md) |
| perturbation | Optional multiplicative primitive or conservative source path. | [Perturbations and Diagnostics](evolution/perturbations-and-diagnostics.md) |
| HydroBase velocity | External Valencia-form velocity exchanged with native GRHayLHD velocity. | [HydroBase Velocity Conversion](integration/hydrobase-velocity-conversion.md) |
| GRHayLib | External library boundary visible through handles, headers, and `ghl_*` calls. | [GRHayLib Contract](integration/grhaylib-contract.md) |
| MoL | External Method-of-Lines registration and RHS schedule boundary. | [ADM, MoL, and Tmunu Contracts](integration/adm-mol-tmunu-contracts.md) |
| TmunuBase | External destination for optional additive stress-energy writes. | [ADM, MoL, and Tmunu Contracts](integration/adm-mol-tmunu-contracts.md) |
| authored test input | Top-level test parfile kept distinct from companion configuration and oracle. | [Test Inventory and Oracles](validation/test-inventory-and-oracles.md) |
| companion configuration | Checked-in generated-configuration assertion without proven oracle lineage. | [Parameters and Configurations](integration/parameters-and-configurations.md) |
| oracle | Checked-in numeric observation, not current pass or validity proof. | [Test Inventory and Oracles](validation/test-inventory-and-oracles.md) |
| coverage gap | Missing local evidence dimension, not proof of defect. | [Coverage Gaps](validation/coverage-gaps.md) |
