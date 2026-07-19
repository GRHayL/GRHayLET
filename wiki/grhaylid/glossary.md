# GRHayLID Glossary

> Page status: reviewed · Last reviewed: 07-19-2026

Each recurring term has one canonical owner. Meanings route queries; they do
not add claims beyond owner-page evidence.

| Term | Routing meaning | Owner |
| --- | --- | --- |
| GRHayLID | Local initial-data thorn and sole domain tree for this branch. | [Purpose and Build Surface](architecture/purpose-and-build-surface.md) |
| initial data (ID) | State setup purpose documented for the local thorn. | [Purpose and Build Surface](architecture/purpose-and-build-surface.md) |
| no owned storage | Absence of locally declared gridfunction groups, storage, and MoL registration. | [Purpose and Build Surface](architecture/purpose-and-build-surface.md) |
| HydroBase_Initial | Declared group receiving initial-hydro schedules. | [Declared Schedule Lifecycle](architecture/schedule-lifecycle.md) |
| HydroBase_Prim2ConInitial | Declared ordering target after beta-equilibrium setup. | [Declared Schedule Lifecycle](architecture/schedule-lifecycle.md) |
| mode | Feature or selection axis used by applicability rows and schedule guards. | [Declared Schedule Lifecycle](architecture/schedule-lifecycle.md) |
| HydroTest1D | HydroBase `initial_hydro` keyword extension selecting local one-dimensional tests. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| Balsara tests | Five named one-dimensional live dispatch arms. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| shock tube | One-dimensional live dispatch arm with a selectable discontinuity. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| sound wave | Declared selection whose pre-loop dispatch arm is commented out. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| equilibrium test | Live one-dimensional arm with its local state setup. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| shock direction | Parameter selecting coordinate rotation for one-dimensional states. | [One-Dimensional Hydro Tests](initial-data/one-d-tests-hydro.md) |
| vector potential | Locally written magnetic potential fields. | [One-Dimensional Magnetic Tests](initial-data/one-d-tests-magnetic.md) |
| staggering | Declared half-cell intent and visible coordinate-use mismatch. | [One-Dimensional Magnetic Tests](initial-data/one-d-tests-magnetic.md) |
| Avec | HydroBase vector-potential destination written by magnetic setup. | [One-Dimensional Magnetic Tests](initial-data/one-d-tests-magnetic.md) |
| Bvec | HydroBase magnetic-field destination written by magnetic setup. | [One-Dimensional Magnetic Tests](initial-data/one-d-tests-magnetic.md) |
| IsotropicGas | Uniform three-dimensional gas setup. | [Isotropic Gas](initial-data/isotropic-gas.md) |
| ConstantDensitySphere | Interior/exterior three-dimensional sphere setup. | [Constant-Density Sphere](initial-data/constant-density-sphere.md) |
| beta equilibrium | Constant-temperature electron-fraction setup for tabulated EOS selection. | [Beta Equilibrium](initial-data/beta-equilibrium.md) |
| atmosphere | Branch using visible atmosphere-handle fields below the local density threshold. | [Beta Equilibrium](initial-data/beta-equilibrium.md) |
| entropy computation | Hybrid and tabulated entropy setup functions and their declared dispatch. | [Entropy Computation](initial-data/entropy-computation.md) |
| keyword extension | Shared HydroBase selector value added by a local `EXTENDS KEYWORD` clause. | [HydroBase Keyword Extensions](integration/hydrobase-keyword-extensions.md) |
| initial_hydro | HydroBase selection extended with three local setup families. | [HydroBase Keyword Extensions](integration/hydrobase-keyword-extensions.md) |
| GRHayLib | External library boundary visible through declarations, handles, and `ghl_*` calls. | [GRHayLib Contract](integration/grhaylib-contract.md) |
| ghl_eos handle | External handle whose fields are visibly read by local functions. | [GRHayLib Contract](integration/grhaylib-contract.md) |
| EOS_type | Shared GRHayLib keyword used by local guards and entropy schedule dispatch. | [GRHayLib Contract](integration/grhaylib-contract.md) |
| CHECK_PARAMETER | Header macro that checks the local `-1` sentinel convention. | [Parameters and Configurations](integration/parameters-and-configurations.md) |
| -1 sentinel | Forbidden real-parameter value used to mark required unset inputs. | [Parameters and Configurations](integration/parameters-and-configurations.md) |
| coverage gap | Missing checked-in evidence dimension, not proof of a defect. | [Coverage Gaps](validation/coverage-gaps.md) |
