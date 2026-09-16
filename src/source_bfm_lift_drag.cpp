#include "source_bfm_lift_drag.hpp"

StateVector SourceBFMLiftDrag::computeBodyForceSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    FlowSolution &conservativeVars,
    FloatType &dt,
    FloatType &timePhysical) {

    computeFlowState(i, j, k, primitive, conservativeVars, timePhysical);
    
    FloatType leIdx = _config.getLeadingEdgeIndex();
    FloatType teIdx = _config.getTrailingEdgeIndex();
    FloatType chord = (
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, teIdx, j, k) - 
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, leIdx, j, k)
    );
    _solidity = chord / _tangentialPitch;

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    
    return inviscidComponent + viscousComponent;
}


StateVector SourceBFMLiftDrag::computeInviscidComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce) {
    
    FloatType w = _relVelCartesian.magnitude();
    FloatType beta0 = _mesh.getInputFields(InputField::LIFT_DRAG_BFM_BETA0, i, j, k);
    FloatType forceMag = (
        0.5 * w*w * 2.0 * M_PI / _staggeredPitch * (-_flowAngle + beta0)
    );
    
    Vector3D forceCylindrical = _inviscidForceDirCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
    
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

StateVector SourceBFMLiftDrag::computeViscousComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &viscousForce) {

    FloatType w = _relVelCartesian.magnitude();
    FloatType beta0 = _mesh.getInputFields(InputField::LIFT_DRAG_BFM_BETA0, i, j, k);
    FloatType Kp0 = _mesh.getInputFields(InputField::LIFT_DRAG_BFM_KP_ETAMAX, i, j, k);
    FloatType Kp = Kp0 + 2.0 * M_PI * _solidity * std::pow(_flowAngle-beta0, 2);
    FloatType forceMag = Kp * w*w / _staggeredPitch;

    Vector3D forceCylindrical = _viscousForceDirCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
    
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}
