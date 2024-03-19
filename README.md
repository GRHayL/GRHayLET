# GRHayLET

---

This repository houses a collection of Einstein Toolkit code modules
("thorns") which depend on GRHayL and are maintained by the GRHayL
developers. Many thorns have two variants--one for the Carpet driver
and another for the CarpetX driver. Those for CarpetX are denoted by
an 'X' suffix.

The initial data thorns (GRHayLID, GRHayLIDX) provide simple initial
data setups, such as the Balsara tests and a constant density sphere.

The purely hydrodynamic evolution thorns (GRHayLHD, GRHayLHDX) can
evolve the GRHD evolution equations with B explicitly set to zero
and are a simplification of the IllinoisGRMHD code.

The GRMHD thorn IllinoisGRMHD solves the GRMHD evolution equations
and is a refactoring of the original IllinoisGRMHD thorn to use
GRHayL. The original thorn is itself a rewrite of the Illinois
Numerical Relativity (ILNR) group's GRMHD code. Details about how
to update old parfiles to use the new thorn are included in the
thorn's documentation.

The NRPyLeakageET thorn (Carpet only) implements a neutrino leakage
scheme in the Toolkit.

More details about individual thorns are given in their
doc/documentation.tex files. To see the compiled versions of this
documentation, please visit the thorn
[documentation page](https://einsteintoolkit.org/documentation/ThornGuide.php)
of the Einstein Toolkit website.
