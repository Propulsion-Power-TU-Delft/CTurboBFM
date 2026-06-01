#include "turbulence_model_sa.hpp"

TurbulenceModelSA::TurbulenceModelSA(
    const Config &config, 
    const FluidBase &fluid, 
    const Mesh &mesh, 
    const std::vector<Boundary> &boundaries,
    const Matrix3D<FloatType> &wallDistance,
    const FlowSolution &initialSolution) 
    : TurbulenceModelBase(config, fluid, mesh, boundaries, initialSolution), 
    _wallDistance(wallDistance) {

        setupModelEquations();
        setupInitialValues();
    }

void TurbulenceModelSA::setupModelEquations() {
    _nuHat.resize(_ni, _nj, _nk);
    _nuHatGrad.resize(_ni, _nj, _nk);
    _nuLaminar.resize(_ni, _nj, _nk);
    _fv1.resize(_ni, _nj, _nk);
    _initNuHat.resize(_ni, _nj, _nk);
}

void TurbulenceModelSA::setupInitialValues() {
    if (_config.restartSolution()) {
        initializeFromRestartFile();
    }
    else{
        initializeFromZero();
    }

}

void TurbulenceModelSA::initializeFromZero() {
    FloatType initTemperature = _config.getInitTemperature();
    FloatType initPressure = _config.getInitPressure();
    FloatType initDensity = _fluid.computeDensity_p_T(initPressure, initTemperature);
    FloatType initNu = _fluid.computeMolecularDynamicViscosity(initTemperature) / initDensity;

    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                _initNuHat(i, j, k) = 0.1 * initNu;
            }
        }
    }
    _nuHat = _initNuHat;
}

void TurbulenceModelSA::initializeFromRestartFile() {

    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                StateVector conservative = _initialSolution.at(i, j, k);
                StateVector primitive = getPrimitiveVariablesFromConservative(conservative);
                FloatType temperature = _fluid.computeTemperature_rho_u_et(
                    primitive[0], 
                    {primitive[1], primitive[2], primitive[3]}, 
                    primitive[4]);
                FloatType nu = _fluid.computeMolecularDynamicViscosity(temperature) / primitive[0];
                _initNuHat(i, j, k) = 0.1 * nu;
            }
        }
    }
    _nuHat = _initNuHat;
}



void TurbulenceModelSA::updateBoundaryValues(const FlowSolution &sol) {
    
    for (auto &boundary : _boundaries) {

        // zero nu hat on the no-slip walls
        if (boundary.type == BoundaryType::NO_SLIP_WALL) {
            enforceNoSlipWallCondition(boundary, sol);
        }
        else if (
            boundary.type == BoundaryType::INLET || 
            boundary.type == BoundaryType::INLET_2D || 
            boundary.type == BoundaryType::INLET_SUPERSONIC){
            enforceInletCondition(boundary, sol);
        }
        else if (boundary.type == BoundaryType::OUTLET || 
            boundary.type == BoundaryType::RADIAL_EQUILIBRIUM || 
            boundary.type == BoundaryType::THROTTLE ||
            boundary.type == BoundaryType::OUTLET_SUPERSONIC){ 
            enforceOutletCondition(boundary, sol);
        }
        else {
            continue;
        }
    }
}

void TurbulenceModelSA::enforceNoSlipWallCondition(const Boundary& boundary, const FlowSolution &sol) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = 0.0;
            }
        }
    }
}

void TurbulenceModelSA::enforceOutletCondition(const Boundary& boundary, const FlowSolution &sol) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);
    BoundaryOrientation orientation = boundary.orientation;
    
    int donorStepI, donorStepJ, donorStepK;
    switch (orientation)
    {
    case BoundaryOrientation::I_START:
        donorStepI = 1;
        donorStepJ = 0;
        donorStepK = 0;
        break;
    case BoundaryOrientation::I_END:
        donorStepI = -1;
        donorStepJ = 0;
        donorStepK = 0;
        break;
    case BoundaryOrientation::J_START:
        donorStepI = 0;
        donorStepJ = 1;
        donorStepK = 0;
        break;
    case BoundaryOrientation::J_END:
        donorStepI = 0;
        donorStepJ = -1;
        donorStepK = 0;
        break;
    case BoundaryOrientation::K_START:
        donorStepI = 0;
        donorStepJ = 0;
        donorStepK = 1;
        break;
    case BoundaryOrientation::K_END:
        donorStepI = 0;
        donorStepJ = 0;
        donorStepK = -1;
        break;
    default:
        break;
    }

    // extract nuHat from the donor cell (internal)
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = _nuHat(i + donorStepI, j + donorStepJ, k + donorStepK);
            }
        }
    }
}

void TurbulenceModelSA::enforceInletCondition(const Boundary& boundary, const FlowSolution &sol) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);

    // inlet values same of the initial solution (or restart solution)
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = _initNuHat(i, j, k);
            }
        }
    }
}

void TurbulenceModelSA::solve(
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad, 
    const Matrix3D<FloatType> &dt){
    
    updateBoundaryValues(sol);
    
    updateSolutionGradient();
    
    updateMeanFlowTerms(sol);

    Matrix3D<FloatType> residual(_ni, _nj, _nk);
    
    computeFluxTerms(residual, sol);
    
    computeSourceTerms(residual, sol, solGrad);
    
    updateSolution(residual, dt);
}

void TurbulenceModelSA::updateSolutionGradient() {
    computeGradientGreenGauss(
        _mesh.getSurfacesI(), 
        _mesh.getSurfacesJ(), 
        _mesh.getSurfacesK(), 
        _mesh.getMidPointsI(), 
        _mesh.getMidPointsJ(), 
        _mesh.getMidPointsK(), 
        _mesh.getVertices(), 
        _mesh.getVolumes(), 
        _nuHat, 
        _nuHatGrad);
}

void TurbulenceModelSA::updateMeanFlowTerms(const FlowSolution &sol) {
    Matrix3D<FloatType> density = sol.getDensity();
    Matrix3D<FloatType> velocityX = sol.getVelocityX();
    Matrix3D<FloatType> velocityY = sol.getVelocityY();
    Matrix3D<FloatType> velocityZ = sol.getVelocityZ();
    Matrix3D<FloatType> totalEnergy = sol.getTotalEnergy();
    
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                FloatType temperature = _fluid.computeTemperature_rho_u_et(
                    density(i, j, k),
                    {velocityX(i, j, k),
                    velocityY(i, j, k),
                    velocityZ(i, j, k)},
                    totalEnergy(i, j, k));
                
                _nuLaminar(i, j, k) = _fluid.computeMolecularDynamicViscosity(temperature);
            }
        }
    }
}

void TurbulenceModelSA::computeFluxTerms(Matrix3D<FloatType> &residual, const FlowSolution &sol) {
    
    computeFluxContribution(FluxDirection::I, sol, residual);
    
    if (_topology!=Topology::ONE_DIMENSIONAL){
        computeFluxContribution(FluxDirection::J, sol, residual);
    }
    
    if (_topology==Topology::THREE_DIMENSIONAL){
        computeFluxContribution(FluxDirection::K, sol, residual);
    }
}

void TurbulenceModelSA::computeFluxContribution(
    FluxDirection direction, 
    const FlowSolution &sol, 
    Matrix3D<FloatType> &residual) {
    
    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);

    Vector3D surface {};
    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();
    size_t dirFace = 0;
    size_t stopFace = 0;
    FloatType advFlux, viscFlux;
    
    for (size_t iFace = 1; iFace < ni-1; ++iFace) {
        for (size_t jFace = 1; jFace < nj-1; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                
                switch (direction)
                {
                case (FluxDirection::I):
                    dirFace = iFace;
                    stopFace = ni-1;
                    break;
                case (FluxDirection::J):
                    dirFace = jFace;
                    stopFace = nj-1;
                    break;
                case (FluxDirection::K):
                    dirFace = kFace;
                    stopFace = nk-1;
                    break;
                default:
                    throw std::runtime_error("Invalid FluxDirection.");
                }
                
                if (dirFace == 0) { // starting boundary interfaces
                    continue;
                }
                else if (dirFace == stopFace) { // ending boundary interfaces
                    continue;
                } 
                else { // flux across internal faces
                    surface = surfaces(iFace, jFace, kFace);
                    advFlux = computeAdvectionFlux(
                        {iFace - 1*stepMask[0], jFace - 1*stepMask[1], kFace - 1*stepMask[2]},
                        {iFace, jFace, kFace},
                        surface,
                        sol);
                    viscFlux = computeViscousFlux(
                        {iFace - 1*stepMask[0], jFace - 1*stepMask[1], kFace - 1*stepMask[2]},
                        {iFace, jFace, kFace},
                        surface);
                    
                    residual(
                        iFace - 1*stepMask[0], 
                        jFace - 1*stepMask[1], 
                        kFace - 1*stepMask[2]
                        ) += (advFlux - viscFlux) * surface.magnitude();
                    
                    residual(iFace, jFace, kFace) -= (advFlux - viscFlux) * surface.magnitude();
                }
            }
        }
    }
}


FloatType TurbulenceModelSA::computeAdvectionFlux(
        const std::array<size_t, 3> &idxLeft, 
        const std::array<size_t, 3> &idxRight,
        const Vector3D &surface,
        const FlowSolution &solution) {
    
    StateVector Uleft = solution.at(idxLeft[0], idxLeft[1], idxLeft[2]);
    StateVector Uright = solution.at(idxRight[0], idxRight[1], idxRight[2]);

    StateVector Wleft = getPrimitiveVariablesFromConservative(Uleft);
    StateVector Wright = getPrimitiveVariablesFromConservative(Uright);

    Vector3D velLeft = Vector3D(Wleft[1], Wleft[2], Wleft[3]);
    Vector3D velRight = Vector3D(Wright[1], Wright[2], Wright[3]);

    FloatType phiLeft = _nuHat(idxLeft[0], idxLeft[1], idxLeft[2]);
    FloatType phiRight = _nuHat(idxRight[0], idxRight[1], idxRight[2]);

    Vector3D n = surface.normalized();
    FloatType un = 0.5 * (velLeft + velRight).dot(n);
    FloatType phiUpwind = (un >= 0.0) ? phiLeft : phiRight;
    FloatType flux = phiUpwind * un;

    return flux;
}

FloatType TurbulenceModelSA::computeViscousFlux(
        const std::array<size_t, 3> &idxLeft, 
        const std::array<size_t, 3> &idxRight,
        const Vector3D &surface) {

    // follow notation blazek pag 228
    Vector3D n = surface.normalized();
    
    Vector3D tauLeft = _nuHatGrad(idxLeft[0], idxLeft[1], idxLeft[2]) / _sigma * (
        _nuLaminar(idxLeft[0], idxLeft[1], idxLeft[2]) + _nuHat(idxLeft[0], idxLeft[1], idxLeft[2]));
    
    Vector3D tauRight = _nuHatGrad(idxRight[0], idxRight[1], idxRight[2]) / _sigma * (
        _nuLaminar(idxRight[0], idxRight[1], idxRight[2]) + _nuHat(idxRight[0], idxRight[1], idxRight[2]));
    
    FloatType flux = 0.5 * (tauLeft + tauRight).dot(n);
    return flux;
}

void TurbulenceModelSA::computeSourceTerms(
        Matrix3D<FloatType> &residual, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad){
    
    for (size_t i = 0; i < _ni; i++){
        for (size_t j = 0; j < _nj; j++){
            for (size_t k = 0; k < _nk; k++){
                FloatType Qt = computeSource({i, j, k}, sol, solGrad);        
                FloatType dV = _mesh.getVolume(i, j, k);
                residual(i, j, k) -= Qt * dV;    
            }
        }
    }
}

FloatType TurbulenceModelSA::computeSource(
    const std::array<size_t, 3> &idx, 
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad){
    
    if (_wallDistance(idx[0], idx[1], idx[2]) == 0.0){
        return 0.0;
    }
    
    size_t i = idx[0];
    size_t j = idx[1];
    size_t k = idx[2];
    
    FloatType chi = _nuHat(i,j,k) / _nuLaminar(i,j,k);
    
    _fv1(i,j,k) = chi*chi*chi / (chi*chi*chi + _cv1*_cv1*_cv1);
    
    FloatType fv2 = std::pow(1.0+chi/_cv2, -3.0);
    
    FloatType fv3 = (1.0 + chi*_fv1(i,j,k)) * (1.0 - fv2) / std::max(chi, 0.001);
    
    FloatType S = computeRotationRateMagnitude(
        solGrad.at(SolutionName::VELOCITY_X)(i,j,k),
        solGrad.at(SolutionName::VELOCITY_Y)(i,j,k),
        solGrad.at(SolutionName::VELOCITY_Z)(i,j,k)
    );
    
    FloatType Shat = fv3*S + _nuHat(i,j,k)*fv2 / (_kappa*_kappa * _wallDistance(i,j,k)*_wallDistance(i,j,k));
    
    FloatType ft2 = _ct3 * std::exp(-_ct4*chi*chi);

    FloatType r = _nuHat(i,j,k) / (Shat * _kappa*_kappa * _wallDistance(i,j,k)*_wallDistance(i,j,k));

    FloatType g = r + _cw2*(std::pow(r,6) - r);

    FloatType fw = g * std::pow(
        (1.0 + std::pow(_cw3, 6.0)) / (std::pow(g, 6.0) + std::pow(_cw3, 6.0)),
        1.0 / 6.0);

    FloatType Qt = _cb1*(1.0-ft2)*Shat*_nuHat(i,j,k) + _cb2/_sigma*(
        _nuHatGrad(i,j,k).x() * _nuHatGrad(i,j,k).x() +
        _nuHatGrad(i,j,k).y() * _nuHatGrad(i,j,k).y() +
        _nuHatGrad(i,j,k).z() * _nuHatGrad(i,j,k).z()
        ) - (_cw1*fw - _cb1*ft2/_kappa/_kappa) * std::pow(_nuHat(i,j,k)/_wallDistance(i,j,k), 2.0); 
    
    return Qt;
}



void TurbulenceModelSA::updateSolution(const Matrix3D<FloatType> &residual, const Matrix3D<FloatType> &dt){
    Matrix3D<FloatType> volume = _mesh.getVolumes();
    _nuHat -= residual * dt / volume;
}


FloatType TurbulenceModelSA::getEddyViscosity(const FloatType &density, size_t i, size_t j, size_t k) const {
    FloatType mu = _fv1(i,j,k) * density * _nuHat(i,j,k);
    return mu;
}