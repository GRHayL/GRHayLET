#include "GRHayLMHD.h"
#include "Symmetry.h"

#define SET_SYMMETRY_NONE(group)          \
    {                                     \
        const int sym[3] = { 1, 1, 1 };   \
        SetCartSymGN(cctkGH, sym, group); \
    }

void GRHayLMHD_RegisterSymmetries(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: add other symmetries
    SET_SYMMETRY_NONE("GRHayLMHD::hydro_conservatives");
    SET_SYMMETRY_NONE("GRHayLMHD::Phi_tilde");
    SET_SYMMETRY_NONE("GRHayLMHD::A_x_tilde");
    SET_SYMMETRY_NONE("GRHayLMHD::A_y_tilde");
    SET_SYMMETRY_NONE("GRHayLMHD::A_z_tilde");
}
