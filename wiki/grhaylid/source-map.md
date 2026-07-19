# GRHayLID Source Map

> Page status: reviewed · Last reviewed: 07-19-2026

This is the sole canonical source-to-page edge table. Registry ingest state is
owned separately and is not promoted merely by adding an edge; governance
issue locators do not stand in for domain claim edges.

Every registered source has at least one domain consumer, and every reviewed
domain leaf has at least one source edge. Governance pages route and constrain
those edges; they do not create independent domain-evidence relationships.

| Source ID | Page ID | Claim kind | Typed locator | Claim status | Next action |
| --- | --- | --- | --- | --- | --- |
| `grhaylid-readme` | `grhaylid.architecture.purpose-and-build-surface` | stated-purpose | `doc:GRHayLID/README#section=1. Purpose` | declared | Recheck when purpose text changes. |
| `grhaylid-thornguide` | `grhaylid.architecture.purpose-and-build-surface` | stated-purpose | `doc:GRHayLID/doc/documentation.tex#section=Introduction` | declared | Keep setup enumeration source-limited. |
| `grhaylid-ccl` | `grhaylid.architecture.purpose-and-build-surface` | cactus-interface | `ccl:GRHayLID/interface.ccl#implementation=GRHayLID` | declared | Reconcile interface boundaries on declaration changes. |
| `grhaylid-ccl` | `grhaylid.architecture.purpose-and-build-surface` | cactus-interface | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` | declared | Keep requirement distinct from build or discovery success. |
| `grhaylid-build` | `grhaylid.architecture.purpose-and-build-surface` | build-surface | `build:GRHayLID/src/make.code.defn#field=SRCS` | declared | Reconcile six source units on manifest changes. |
| `grhaylid-header` | `grhaylid.architecture.purpose-and-build-surface` | visible-formula | `macro:GRHayLID/src/GRHayLID.h#name=CHECK_PARAMETER` | visible-implementation | Keep macro semantics limited to visible expansion. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_hydro_data` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_1D_tests_magnetic_data` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_IsotropicGas` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_ConstantDensitySphere` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` | declared | Reconcile guard, ordering, and field declarations on schedule changes. |
| `grhaylid-ccl` | `grhaylid.architecture.schedule-lifecycle` | parameter | `ccl:GRHayLID/param.ccl#parameter=initialize_magnetic_quantities` | declared | Reconcile nested magnetic guard on parameter changes. |
| `grhaylid-c` | `grhaylid.initial-data.one-d-tests-hydro` | visible-call-order | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | visible-implementation | Recheck live dispatch arms and loop ordering on source changes. |
| `grhaylid-c` | `grhaylid.initial-data.one-d-tests-hydro` | visible-call-order | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | unresolved | Track the commented sound-wave dispatch arm in GID-0005. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-hydro` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | declared | Reconcile selector meaning and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-hydro` | parameter | `ccl:GRHayLID/param.ccl#parameter=shock_direction` | declared | Reconcile selector meaning and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-hydro` | parameter | `ccl:GRHayLID/param.ccl#parameter=discontinuity_position` | declared | Reconcile selector meaning and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-hydro` | parameter | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` | declared | Reconcile selector meaning and visible consumer. |
| `grhaylid-c` | `grhaylid.initial-data.one-d-tests-hydro` | external-behavior | `c:GRHayLID/src/1D_tests_hydro_data.c#call=ghl_hybrid_find_polytropic_index?function=GRHayLID_1D_tests_hydro_data` | out-of-scope | Keep external call semantics delegated. |
| `grhaylid-c` | `grhaylid.initial-data.one-d-tests-magnetic` | visible-dataflow | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` | visible-implementation | Recheck coordinates, rotations, and writes on source changes. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-magnetic` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` | declared | Reconcile selector guard and write surface. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-magnetic` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_Bvec` | declared | Reconcile selector guard and write surface. |
| `grhaylid-ccl` | `grhaylid.initial-data.one-d-tests-magnetic` | parameter | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` | unresolved | Track absent visible consumer in GID-0004. |
| `grhaylid-c` | `grhaylid.initial-data.isotropic-gas` | visible-dataflow | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` | visible-implementation | Recheck uniform state and writes on source changes. |
| `grhaylid-ccl` | `grhaylid.initial-data.isotropic-gas` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_rho` | declared | Reconcile parameter and visible sentinel consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.isotropic-gas` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_Y_e` | declared | Reconcile parameter and visible sentinel consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.isotropic-gas` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_temperature` | declared | Reconcile parameter and visible sentinel consumer. |
| `grhaylid-thornguide` | `grhaylid.initial-data.isotropic-gas` | stated-purpose | `doc:GRHayLID/doc/documentation.tex#section=IsotropicGas` | declared | Keep stated purpose distinct from numerical validation. |
| `grhaylid-c` | `grhaylid.initial-data.constant-density-sphere` | visible-dataflow | `c:GRHayLID/src/ConstantDensitySphere.c#symbol=GRHayLID_ConstantDensitySphere` | visible-implementation | Recheck radius branch, state, and writes on source changes. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_sphere_radius` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vx_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vy_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vz_interior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_exterior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_exterior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.constant-density-sphere` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_exterior` | declared | Reconcile parameter and visible consumer. |
| `grhaylid-thornguide` | `grhaylid.initial-data.constant-density-sphere` | stated-purpose | `doc:GRHayLID/doc/documentation.tex#section=ConstantDensitySphere` | declared | Keep stated purpose distinct from numerical validation. |
| `grhaylid-c` | `grhaylid.initial-data.beta-equilibrium` | visible-call-order | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` | visible-implementation | Recheck bounds, sentinel order, and writes. |
| `grhaylid-ccl` | `grhaylid.initial-data.beta-equilibrium` | parameter | `ccl:GRHayLID/param.ccl#parameter=impose_beta_equilibrium` | declared | Reconcile parameter and schedule/body consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.beta-equilibrium` | parameter | `ccl:GRHayLID/param.ccl#parameter=beq_temperature` | declared | Reconcile parameter and schedule/body consumer. |
| `grhaylid-ccl` | `grhaylid.initial-data.beta-equilibrium` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` | declared | Reconcile schedule description and declared writes. |
| `grhaylid-c` | `grhaylid.initial-data.beta-equilibrium` | external-behavior | `c:GRHayLID/src/BetaEquilibrium.c#call=ghl_tabulated_compute_Ye_of_rho_beq_constant_T?function=GRHayLID_BetaEquilibrium` | out-of-scope | Keep external table-call semantics delegated. |
| `grhaylid-readme` | `grhaylid.initial-data.beta-equilibrium` | stated-purpose | `doc:GRHayLID/README#section=1. Purpose` | declared | Keep standalone-use wording source-limited. |
| `grhaylid-c` | `grhaylid.initial-data.entropy-computation` | visible-dataflow | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_hybrid` | visible-implementation | Recheck visible inputs, in-place writes, and outputs. |
| `grhaylid-c` | `grhaylid.initial-data.entropy-computation` | visible-dataflow | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` | visible-implementation | Recheck visible inputs, in-place writes, and outputs. |
| `grhaylid-ccl` | `grhaylid.initial-data.entropy-computation` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` | declared | Reconcile keyword extension and schedule guard. |
| `grhaylid-ccl` | `grhaylid.initial-data.entropy-computation` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` | declared | Reconcile EOS dispatch and declared writes. |
| `grhaylid-ccl` | `grhaylid.initial-data.entropy-computation` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` | declared | Reconcile EOS dispatch and declared writes. |
| `grhaylid-readme` | `grhaylid.initial-data.entropy-computation` | stated-purpose | `doc:GRHayLID/README#section=1. Purpose` | unresolved | Track any-EOS wording in GID-0009. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_hydro` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_Y_e` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_temperature` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_Bvec` | declared | Reconcile extension value and local consumer. |
| `grhaylid-ccl` | `grhaylid.integration.hydrobase-keyword-extensions` | cactus-interface | `ccl:GRHayLID/interface.ccl#implementation=GRHayLID` | declared | Reconcile HydroBase sharing and write surface. |
| `grhaylid-thornguide` | `grhaylid.integration.hydrobase-keyword-extensions` | stated-purpose | `doc:GRHayLID/doc/documentation.tex#section=Parameters` | declared | Keep documented selections aligned with declarations. |
| `grhaylid-c` | `grhaylid.integration.hydrobase-keyword-extensions` | visible-dataflow | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` | visible-implementation | Recheck one in-body keyword precondition. |
| `grhaylid-ccl` | `grhaylid.integration.grhaylib-contract` | parameter | `ccl:GRHayLID/param.ccl#parameter=EOS_type` | declared | Reconcile shared keyword consumers. |
| `grhaylid-header` | `grhaylid.integration.grhaylib-contract` | visible-dataflow | `macro:GRHayLID/src/GRHayLID.h#include=GRHayLib.h` | visible-implementation | Recheck include boundary. |
| `grhaylid-ccl` | `grhaylid.integration.grhaylib-contract` | cactus-interface | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` | declared | Keep requirement distinct from discovery and link success. |
| `grhaylid-c` | `grhaylid.integration.grhaylib-contract` | visible-call-order | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` | visible-implementation | Recheck handle reads and local call/error order. |
| `grhaylid-c` | `grhaylid.integration.grhaylib-contract` | external-behavior | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` | out-of-scope | Keep library and EOS-table semantics delegated. |
| `grhaylid-c` | `grhaylid.integration.grhaylib-contract` | visible-call-order | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` | visible-implementation | Recheck handle reads and local call/error order. |
| `grhaylid-c` | `grhaylid.integration.grhaylib-contract` | external-behavior | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` | out-of-scope | Keep library and EOS-table semantics delegated. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=initialize_magnetic_quantities` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=shock_direction` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=discontinuity_position` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_rho` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_Y_e` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_temperature` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_sphere_radius` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vx_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vy_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_vz_interior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_rho_exterior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_Y_e_exterior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=ConstantDensitySphere_T_exterior` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=impose_beta_equilibrium` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=beq_temperature` | declared | Reconcile declaration, default, meaning, and visible consumer. |
| `grhaylid-ccl` | `grhaylid.integration.parameters-and-configurations` | parameter | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` | unresolved | Track absent visible consumer in GID-0004. |
| `grhaylid-header` | `grhaylid.integration.parameters-and-configurations` | visible-formula | `macro:GRHayLID/src/GRHayLID.h#name=CHECK_PARAMETER` | visible-implementation | Recheck sentinel macro and call sites. |
| `grhaylid-ccl` | `grhaylid.validation.coverage-gaps` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | coverage-gap | Add checked regression evidence across declared selections. |
| `grhaylid-ccl` | `grhaylid.validation.coverage-gaps` | parameter | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` | coverage-gap | Add checked staggering evidence. |
| `grhaylid-ccl` | `grhaylid.validation.coverage-gaps` | schedule-intent | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` | coverage-gap | Add checked beta-equilibrium configuration and regression evidence. |
| `grhaylid-c` | `grhaylid.validation.coverage-gaps` | visible-dataflow | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | coverage-gap | Add checked dispatch evidence, including sound-wave selection. |
| `grhaylid-build` | `grhaylid.validation.coverage-gaps` | build-surface | `build:GRHayLID/src/make.code.defn#field=SRCS` | coverage-gap | Re-audit when checked-in test artifacts appear. |
| `grhaylid-ccl` | `grhaylid.validation.coverage-gaps` | parameter | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` | coverage-gap | Add checked entropy-dispatch configuration evidence. |
