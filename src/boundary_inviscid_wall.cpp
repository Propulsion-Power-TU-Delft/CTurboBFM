#include "boundary_inviscid_wall.hpp"
#include "math_utils.hpp"


StateVector BoundaryInviscidWall::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
                
    StateVector primitive = getPrimitiveVariablesFromConservative(internalConservative);
    Vector3D velocity({primitive[1], primitive[2], primitive[3]});
    FloatType pressure = _fluid.computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    StateVector flux {{0,0,0,0,0}};
    auto Sdir = surface / surface.magnitude();
    FloatType un = velocity.dot(Sdir);
    FloatType a = _fluid.computeSoundSpeed_p_rho(pressure, primitive[0]);
    FloatType pWall = pressure + primitive[0] * un * a;
    if (pWall < 0.01 * pressure) {
        pWall = 0.01 * pressure;
    }
    flux[0] = 0.0;
    flux[1] = pWall * Sdir.x();
    flux[2] = pWall * Sdir.y();
    flux[3] = pWall * Sdir.z();
    flux[4] = 0.0;
    return flux;
}