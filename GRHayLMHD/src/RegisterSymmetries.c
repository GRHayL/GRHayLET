#include "GRHayLMHD.h"
#include "Symmetry.h"

#define SET_GROUP_SYMMETRY_TO_SCALAR(group) \
    {                                       \
        const int sym[3] = { 1, 1, 1 };     \
        SetCartSymGN(cctkGH, sym, group);   \
    }

#define SET_GF_SYMMETRY_TO_VEC(group, _x, _y, _z) \
    {                                             \
        const int sym[3] = { _x, _y, _z };        \
        SetCartSymVN(cctkGH, sym, group);         \
    }

#define SET_GF_SYMMETRY_TO_VECX(group) SET_GF_SYMMETRY_TO_VEC(group, -1, +1, +1)
#define SET_GF_SYMMETRY_TO_VECY(group) SET_GF_SYMMETRY_TO_VEC(group, +1, -1, +1)
#define SET_GF_SYMMETRY_TO_VECZ(group) SET_GF_SYMMETRY_TO_VEC(group, +1, +1, -1)

#define SET_GROUP_BC_TO_COPY(group)                                                 \
    if(Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, group, "copy") < 0) { \
        CCTK_VERROR("Failed to register BC for group %s", group);                   \
    }

void GRHayLMHD_RegisterSymmetries(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Register symmetries for hydrodynamics quantities
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::rho");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::press");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::eps");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::Y_e");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::temperature");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::entropy");
    SET_GROUP_SYMMETRY_TO_SCALAR("HydroBase::w_lorentz");
    SET_GROUP_SYMMETRY_TO_SCALAR("GRHayLMHD::hydro_conservatives");
    SET_GF_SYMMETRY_TO_VECX("GRHayLMHD::vx");
    SET_GF_SYMMETRY_TO_VECY("GRHayLMHD::vy");
    SET_GF_SYMMETRY_TO_VECZ("GRHayLMHD::vz");

    // Not sure we can support symmetries for staggered GFs
    // SET_GROUP_SYMMETRY_TO_SCALAR("GRHayLMHD::Phi_tilde");
    // SET_GROUP_SYMMETRY_TO_SCALAR("GRHayLMHD::A_x_tilde");
    // SET_GROUP_SYMMETRY_TO_SCALAR("GRHayLMHD::A_y_tilde");
    // SET_GROUP_SYMMETRY_TO_SCALAR("GRHayLMHD::A_z_tilde");

    // TODO: test if this works
    // Register boundary conditions for hydrodynamics quantities
    // if(!CCTK_IsFunctionAliased("Driver_SelectVarForBC")) {
    //     CCTK_WARN(CCTK_WARN_ALERT, ">>>>> Driver_SelectVarForBC not aliased, SO NO BCS <<<<<");
    //     return;
    // }
    // SET_GROUP_BC_TO_COPY("HydroBase::rho");
    // SET_GROUP_BC_TO_COPY("HydroBase::press");
    // SET_GROUP_BC_TO_COPY("HydroBase::eps");
    // SET_GROUP_BC_TO_COPY("HydroBase::Y_e");
    // SET_GROUP_BC_TO_COPY("HydroBase::temperature");
    // SET_GROUP_BC_TO_COPY("HydroBase::entropy");
    // SET_GROUP_BC_TO_COPY("HydroBase::w_lorentz");
    // SET_GROUP_BC_TO_COPY("GRHayLMHD::hydro_conservatives");
    // SET_GROUP_BC_TO_COPY("GRHayLMHD::hydro_velocities");
}
