#include "boundary_inlet_supersonic.hpp"
#include "math_utils.hpp"

StateVector BoundaryInletSupersonic::computeBoundaryFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const Vector3D& midPoint, 
    const std::array<size_t, 3>& indices, 
    const FlowSolution& flowSolution, 
    const size_t& iterCounter) {

    FloatType pressureBoundary = _boundaryValues[0];
    FloatType temperatureBoundary = _boundaryValues[1];

    if (isNodeOnNoSlipWall(indices[0], indices[1], indices[2])) {
        Vector3D uWall = getWallVelocity(indices[0], indices[1], indices[2]);
        return computeStagnationWallFlux(surface, pressureBoundary, temperatureBoundary, uWall);
    }

    Vector3D velocityBoundary({_boundaryValues[2], _boundaryValues[3], _boundaryValues[4]});
    FloatType densityBoundary = _fluid.computeDensity_p_T(pressureBoundary, temperatureBoundary);
    FloatType energyBoundary = _fluid.computeStaticEnergy_p_rho(pressureBoundary, densityBoundary);
    FloatType totEnergyBoundary = energyBoundary + 0.5 * velocityBoundary.dot(velocityBoundary);
    StateVector primitiveBoundary({
        densityBoundary, 
        velocityBoundary.x(), 
        velocityBoundary.y(), 
        velocityBoundary.z(), 
        totEnergyBoundary
    });

    auto flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}
