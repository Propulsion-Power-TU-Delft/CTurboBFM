#include "boundary_outlet_mass_flow.hpp"
#include "math_utils.hpp"

StateVector BoundaryOutletMassFlow::computeBoundaryFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const Vector3D& midPoint, 
    const std::array<size_t, 3>& indices, 
    const FlowSolution& flowSolution, 
    const size_t& iterCounter) {

    size_t idx = indices[1];
    if (idx >= _pressureProfile.size()) {
        idx = _pressureProfile.size() - 1;
    }
    FloatType outletPressure = _pressureProfile[idx];
    auto flux = computeOutletFlux(internalConservative, surface, outletPressure, iterCounter);
    return flux;
}
