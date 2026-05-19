#include "boundary_inlet.hpp"
#include "math_utils.hpp"

StateVector BoundaryInlet::computeBoundaryFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const Vector3D& midPoint, 
    const std::array<size_t, 3> &indices, 
    const FlowSolution& flowSolution, 
    const size_t& iterCounter) {    
    
    // gather boundary values from config
    FloatType totPressure = _boundaryValues[0];
    FloatType totTemperature = _boundaryValues[1];
    Vector3D flowDirection = {_boundaryValues[2], _boundaryValues[3], _boundaryValues[4]};

    if (_referenceFrame == ReferenceFrame::CYLINDRICAL){
        FloatType theta = _mesh.getTheta(indices[0], indices[1], indices[2]);
        flowDirection = computeCartesianComponentsFromCylindrical(flowDirection, theta);
    }

    if (flowDirection.x() == 1.0 && flowDirection.y() == 1.0 && flowDirection.z() == 1.0) {
        // in this case the direction is normal to the surface (just code convention)
        flowDirection = -surface;
    }

    flowDirection /= flowDirection.magnitude();

    StateVector flux = computeSubsonicInletFlux(
        internalConservative, 
        surface, 
        totPressure, 
        totTemperature, 
        flowDirection);
    
    return flux;
}
