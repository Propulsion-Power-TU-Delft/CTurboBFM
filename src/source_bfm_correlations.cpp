#include "source_bfm_correlations.hpp"

SourceBFMCorrelations::SourceBFMCorrelations(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh) : SourceBFMBase(config, fluid, mesh) {
        _leadingEdgeIdx = _config.getLeadingEdgeIndex();
        _trailingEdgeIdx = _config.getTrailingEdgeIndex();
    }


StateVector SourceBFMCorrelations::computeBodyForceSource(
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
    computeCorrelationParameters(i, j, k, conservativeVars);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);

    return inviscidComponent + viscousComponent;
}


void SourceBFMCorrelations::computeCorrelationParameters(
    size_t i, 
    size_t j, 
    size_t k, 
    FlowSolution &conservativeSolution) {

    _solidity = (_mesh.getInputFields(InputField::STREAMWISE_LENGTH, _trailingEdgeIdx, j, k) - 
                 _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _leadingEdgeIdx, j, k)) / _tangentialPitch;
    
    StateVector conservativeInlet = conservativeSolution.at(_leadingEdgeIdx, j, k);
    StateVector conservativeOutlet = conservativeSolution.at(_trailingEdgeIdx, j, k);

    StateVector primitiveInlet = getPrimitiveVariablesFromConservative(conservativeInlet);
    StateVector primitiveOutlet = getPrimitiveVariablesFromConservative(conservativeOutlet);

    Vector3D velocityInlet = {primitiveInlet[1], primitiveInlet[2], primitiveInlet[3]};
    Vector3D velocityOutlet = {primitiveOutlet[1], primitiveOutlet[2], primitiveOutlet[3]};

    Vector3D velocityInletCyl = computeCylindricalComponentsFromCartesian(velocityInlet, _theta);
    Vector3D velocityOutletCyl = computeCylindricalComponentsFromCartesian(velocityOutlet, _theta);

    FloatType radiusInlet = _mesh.getRadius(_leadingEdgeIdx, j, k);
    FloatType radiusOutlet = _mesh.getRadius(_trailingEdgeIdx, j, k);

    Vector3D dragVelocityInlet = {0.0, 0.0, _omega * radiusInlet};
    Vector3D dragVelocityOutlet = {0.0, 0.0, _omega * radiusOutlet};
    
    _relVelInlet = velocityInletCyl - dragVelocityInlet;
    _relVelOutlet = velocityOutletCyl - dragVelocityOutlet;

    _diffusionFactor = 1.0 - _relVelOutlet.magnitude() / _relVelInlet.magnitude() + std::abs(
                                _relVelInlet.z() - _relVelOutlet.z()) / 
                                (2.0 * _solidity * _relVelInlet.magnitude());
    
    FloatType velMeridionalInlet = std::sqrt(_relVelInlet.x() * _relVelInlet.x() + 
                                             _relVelInlet.y() * _relVelInlet.y());

    _inletFlowAngleDeg = std::atan2(_relVelInlet.z(), velMeridionalInlet) * 180.0 / M_PI;
    
    _bladeMeridionalLength = _mesh.getInputFields(
        InputField::STREAMWISE_LENGTH, _trailingEdgeIdx, j, k) - 
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _leadingEdgeIdx, j, k);

    _camberAngleDegAbs = std::abs(
        _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, _trailingEdgeIdx, j, k) -
         _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, _leadingEdgeIdx, j, k))  * 180 / M_PI;
}




StateVector SourceBFMCorrelations::computeInviscidComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce) {
    
    // formulation found in Aungier book, corrected with Cetin compressibility effects
    FloatType kappa1Deg = _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, _leadingEdgeIdx, j, k) * 180 / M_PI; 
    FloatType kappa2Deg = _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, _trailingEdgeIdx, j, k) * 180 / M_PI; 
    
    FloatType incidenceInletDeg = _inletFlowAngleDeg - kappa1Deg; 

    // Keep in mind that correlations formula and values are all supposed to have positive values 
    FloatType _flowAngleInletDegAbs = std::abs(_inletFlowAngleDeg);
    
    FloatType n = 0.025*_solidity - 0.06 - (
        std::pow(_flowAngleInletDegAbs/90.0, 1.0+1.2*_solidity)) / (1.5 + 0.43*_solidity);
    
    FloatType p = 0.914 + _solidity*_solidity*_solidity / 160.0;
    
    _incidenceZeroTenThkStar = std::pow(_flowAngleInletDegAbs,p) / (
        5.0 + 46.0 * std::exp(-2.3*_solidity)) - 0.1*std::pow(_solidity, 3.0) *
         std::exp((_flowAngleInletDegAbs-70.0) / 4.0);

    // find thickness max along the stream
    FloatType thkMax = 0.0;
    for (size_t streamIdx = _leadingEdgeIdx; streamIdx <= _trailingEdgeIdx; ++streamIdx) {
        FloatType r = _mesh.getRadius(streamIdx, j, k);
        FloatType pitch = 2.0 * M_PI * r / _bladesNumber;
        FloatType blk = _mesh.getInputFields(InputField::BLOCKAGE, streamIdx, j, k);
        FloatType thkTmp = pitch * (1.0 - blk);
        if (thkTmp > thkMax) {
            thkMax = thkTmp;
        }
    }
    
    FloatType tcRatio = thkMax / _bladeMeridionalLength;
    FloatType Kti = std::pow(10.0*tcRatio, 0.28/(0.1 + std::pow(tcRatio,0.3)));
    FloatType Ksh = 0.8; // blade shape factor
    FloatType istar = Ksh*Kti*_incidenceZeroTenThkStar + n*_camberAngleDegAbs;
    
    Vector3D absVel({primitive[1], primitive[2], primitive[3]});
    FloatType soundSpeed = _fluid.computeSoundSpeed_rho_u_et(primitive[0], absVel, primitive[4]);
    FloatType inletMach = _relVelCartesian.magnitude() / soundSpeed;

    FloatType istarCorrected = istar + 1.3026*inletMach + 5.7380; // compressibility correction by Cetin

    FloatType deviation10OutStar = 0.01*_solidity*_flowAngleInletDegAbs + (
        0.74*std::pow(_solidity,1.9) + 3.0*_solidity) *
        std::pow(_flowAngleInletDegAbs/90.0, 1.67+1.09*_solidity);
    
    FloatType Ktd = 6.25*tcRatio + 37.5*tcRatio*tcRatio;
    FloatType acRatio = 0.3; // max thickness location / chord length --> guess it
    FloatType deviationStar = (0.92*acRatio*acRatio + 0.002*std::abs(kappa2Deg)) / (
        1-0.002*_camberAngleDegAbs/std::sqrt(_solidity)) * 
        _camberAngleDegAbs/std::sqrt(_solidity) + (Ksh*Ktd-1.0)*deviation10OutStar;
    
    FloatType deviationStarCorrected = -1.099379 + 3.0186*deviationStar - 0.1988*deviationStar*deviationStar;

    // off-design calculation taken by aungier
    FloatType dd_distar = (1.0 + (_solidity+0.25*std::pow(_solidity, 4.0)) * 
        std::pow((_flowAngleInletDegAbs/53.0), 2.5)) / std::exp(3.1*_solidity);
    
    FloatType umInlet = std::sqrt(_relVelInlet.x()*_relVelInlet.x() + _relVelInlet.y()*_relVelInlet.y());
    FloatType umOutlet = std::sqrt(_relVelOutlet.x()*_relVelOutlet.x() + _relVelOutlet.y()*_relVelOutlet.y());

    FloatType incidenceInletDegMag = 0.0;
    // incidence considered positive when positive flow turning needs to be applied
    if (_omega>0.0) {
        incidenceInletDegMag = incidenceInletDeg*(-1.0); 
    }
    else {
        incidenceInletDegMag = incidenceInletDeg;
    }
    FloatType deviationOutletDeg = deviationStarCorrected + dd_distar*(
        incidenceInletDegMag - istarCorrected) + 10.0*(1.0 - umOutlet/umInlet);

    // curvature radius of camber line (clipped to avoid too high values)
    FloatType camberCurvature = std::abs(1.0/_mesh.getInputFields(InputField::BLADE_CAMBER_CURVATURE, i, j, k)); 
    FloatType curvatureRadius = 1.0 / camberCurvature;
    if (curvatureRadius>2) {
        curvatureRadius = 2.0; 
    }

    FloatType normalizedLocation = (
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i,j,k) - 
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _leadingEdgeIdx,j,k)) / (
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _trailingEdgeIdx,j,k) - 
        _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _leadingEdgeIdx,j,k));
    
    FloatType deviationExpectedDeg = incidenceInletDegMag + normalizedLocation * (
        deviationOutletDeg - incidenceInletDegMag);
    FloatType deviationExpectedRad = deviationExpectedDeg * M_PI / 180.0;

    FloatType uThetaCurrent = _relVelCylindric.z();
    FloatType kappa = _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, i, j, k); 
    FloatType Wn, uThetaExpected;
    
    if (_omega>0.0) {
        uThetaExpected = _velMeridional * std::tan(-deviationExpectedRad + kappa);
        Wn = uThetaExpected - uThetaCurrent; 
    }
    else {
        uThetaExpected = _velMeridional * std::tan(deviationExpectedRad + kappa);
        Wn = -uThetaExpected + uThetaCurrent;
    }
    
    FloatType KN = _config.getKnCorrelationBfmCoefficient(); 
    FloatType fnMag = KN * _velMeridional * Wn / (_tangentialPitch * std::abs(std::cos(_metalAngle)))  ;

    Vector3D forceCartesian = _inviscidForceDirCartesian * fnMag;
    Vector3D forceCylindrical = _inviscidForceDirCylindrical * fnMag;

    // update the force in memory
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}


StateVector SourceBFMCorrelations::computeViscousComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &viscousForce) {
    
    // zero losses for the moment
    Vector3D forceCartesian = {0,0,0};
    Vector3D forceCylindrical = {0,0,0};
    
    // update the force in memory
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}
