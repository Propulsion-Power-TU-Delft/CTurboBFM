#include "source_bfm_thollet.hpp"

StateVector SourceBFMThollet::computeBodyForceSource(
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

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    
    return inviscidComponent + viscousComponent;
}


StateVector SourceBFMThollet::computeInviscidComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce) {

    FloatType Kmach = computeCompressibilityCorrection(_relVelCylindric, primitive);
    FloatType forceMag = _Kn * Kmach * _relVelCylindric.dot(_relVelCylindric) * M_PI * _deviationAngle /
        _tangentialPitch / std::abs(_normalCamberTangential) / _blockage;
    
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

FloatType SourceBFMThollet::computeCompressibilityCorrection(
    const Vector3D& relativeVelocityCylindric, 
    const StateVector& primitive) {

    FloatType soundSpeed = _fluid.computeSoundSpeed_rho_u_et(
        primitive[0], 
        {primitive[1], primitive[2], primitive[3]}, 
        primitive[4]);
    FloatType relativeMach = relativeVelocityCylindric.magnitude() / soundSpeed;

    if (std::abs(relativeMach-1.0) < 1E-03)
        relativeMach = 0.999;
    
    FloatType kPrime = 0.0;
    if (relativeMach<1.0)
        kPrime = 1.0 / std::sqrt(1.0 - relativeMach*relativeMach);
    else
        kPrime = 2.0 / M_PI / std::sqrt(relativeMach*relativeMach - 1.0);

    return std::min(3.0, kPrime);
}

StateVector SourceBFMThollet::computeViscousComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &viscousForce) {
    
    if (_isStalledBfmActive){
        FloatType deviationAngleStall = _mesh.getInputFields(InputField::DEVIATION_ANGLE_STALL, i, j, k);
        if (_deviationAngle > deviationAngleStall){
            return StateVector({0,0,0,0,0});
        }
    }
    
    FloatType nu = _config.getFluidKinematicViscosity();
    FloatType stwl = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i, j, k);

    if (stwl < 1E-06){
        if (i + 1 < _mesh.getNumberPointsI()){
            stwl = (_mesh.getInputFields(InputField::STREAMWISE_LENGTH, i+1, j, k) + 
                    _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i, j, k)) / 2.0;
        }
        if (stwl < 1E-06){
            stwl = 1E-06;
        }
    }

    FloatType ReX = _relVelCylindric.magnitude() * stwl / nu;
    FloatType Cf = 0.0592 * std::pow(ReX, -0.2);

    FloatType forceMag = 0.0;
    if (!_isOffDesignActive){
        forceMag = _relVelCylindric.dot(_relVelCylindric) / (
            _tangentialPitch * _blockage * std::abs(_normalCamberTangential)) * (_Kf * Cf);
        }
    else{
        FloatType deviationAnglePivot = _mesh.getInputFields(InputField::DEVIATION_ANGLE_PIVOT, i, j, k);
        FloatType Kmach = computeCompressibilityCorrection(_relVelCylindric, primitive);
        forceMag = _relVelCylindric.dot(_relVelCylindric) / (
            _tangentialPitch * _blockage * std::abs(_normalCamberTangential)) * (
            _Kf * Cf + M_PI * Kmach * std::pow(std::abs(_deviationAngle - deviationAnglePivot), 2.0 * _Kd));
    }

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
