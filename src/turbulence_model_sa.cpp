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
}

void TurbulenceModelSA::setupInitialValues() {
    FloatType initTemperature = _config.getInitTemperature();
    FloatType initPressure = _config.getInitPressure();
    FloatType initDensity = _fluid.computeDensity_p_T(initPressure, initTemperature);
    _initNu = _fluid.computeMolecularDynamicViscosity(initTemperature) / initDensity;

    if (_config.restartSolution()) {
        initializeFromRestartFile();
    }
    else{
        initializeFromZero();
    }

}

void TurbulenceModelSA::initializeFromZero() {

    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                _nuHat(i, j, k) = _initialNuHatScaling * _initNu;
            }
        }
    }
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
                _nuHat(i, j, k) = _initialNuHatScaling * nu;
            }
        }
    }
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
            boundary.type == BoundaryType::INLET_SUPERSONIC ||
            boundary.type == BoundaryType::FARFIELD){
            enforceInletCondition(boundary, sol);
        }
        else if (boundary.type == BoundaryType::OUTLET || 
            boundary.type == BoundaryType::RADIAL_EQUILIBRIUM || 
            boundary.type == BoundaryType::THROTTLE ||
            boundary.type == BoundaryType::OUTLET_SUPERSONIC ||
            boundary.type == BoundaryType::INVISCID_WALL ||
            boundary.type == BoundaryType::TRANSPARENT){ 
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
    
    int donorStepI=0, donorStepJ=0, donorStepK=0;
    
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
                _nuHat(i, j, k) = _initNu * _farfieldNuHatScaling;
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
                
                _nuLaminar(i, j, k) = _fluid.computeMolecularDynamicViscosity(temperature) / density(i, j, k);
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
    
    for (size_t iFace = 0; iFace < ni; ++iFace) {
        for (size_t jFace = 0; jFace < nj; ++jFace) {
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
                
                if (dirFace == 0 || dirFace == stopFace) { // skip flux on boundaries
                    continue;
                } 
                else { 
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

    Vector3D n = surface.normalized();
    
    // Left state
    FloatType nuL_left = _nuLaminar(idxLeft[0], idxLeft[1], idxLeft[2]);
    FloatType nuHat_left = _nuHat(idxLeft[0], idxLeft[1], idxLeft[2]);
    FloatType nuTotal_left;
    if (nuHat_left >= 0.0) {
        nuTotal_left = nuL_left + nuHat_left;
    } else {
        FloatType chi = (nuL_left > 0.0) ? (nuHat_left / nuL_left) : 0.0;
        FloatType chi3 = chi * chi * chi;
        FloatType fn = (_cn1 + chi3) / (_cn1 - chi3);
        nuTotal_left = nuL_left + nuHat_left * fn;
    }
    Vector3D tauLeft = (_nuHatGrad(idxLeft[0], idxLeft[1], idxLeft[2]) / _sigma) * nuTotal_left;

    // Right state
    FloatType nuL_right = _nuLaminar(idxRight[0], idxRight[1], idxRight[2]);
    FloatType nuHat_right = _nuHat(idxRight[0], idxRight[1], idxRight[2]);
    FloatType nuTotal_right;
    if (nuHat_right >= 0.0) {
        nuTotal_right = nuL_right + nuHat_right;
    } else {
        FloatType chi = (nuL_right > 0.0) ? (nuHat_right / nuL_right) : 0.0;
        FloatType chi3 = chi * chi * chi;
        FloatType fn = (_cn1 + chi3) / (_cn1 - chi3);
        nuTotal_right = nuL_right + nuHat_right * fn;
    }
    Vector3D tauRight = (_nuHatGrad(idxRight[0], idxRight[1], idxRight[2]) / _sigma) * nuTotal_right;

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
    
    size_t i = idx[0];
    size_t j = idx[1];
    size_t k = idx[2];

    FloatType d = _wallDistance(i, j, k);
    if (d <= 0.0){
        return 0.0;
    }
    
    FloatType nuHat = _nuHat(i, j, k);
    FloatType nuL = _nuLaminar(i, j, k);
    if (nuL <= 0.0) {
        nuL = 1e-14;
    }

    FloatType S_vort = computeRotationRateMagnitude(
        solGrad.at(SolutionName::VELOCITY_X)(i, j, k),
        solGrad.at(SolutionName::VELOCITY_Y)(i, j, k),
        solGrad.at(SolutionName::VELOCITY_Z)(i, j, k)
    );

    FloatType gradNuHatSq = _nuHatGrad(i, j, k).x() * _nuHatGrad(i, j, k).x() +
                            _nuHatGrad(i, j, k).y() * _nuHatGrad(i, j, k).y() +
                            _nuHatGrad(i, j, k).z() * _nuHatGrad(i, j, k).z();

    FloatType crossDiff = (_cb2 / _sigma) * gradNuHatSq;

    if (nuHat >= 0.0) {
        FloatType chi = nuHat / nuL;
        FloatType chi3 = chi * chi * chi;
        FloatType cv1_3 = _cv1 * _cv1 * _cv1;
        FloatType fv1 = chi3 / (chi3 + cv1_3);
        _fv1(i, j, k) = fv1;

        FloatType fv2 = 1.0 - chi / (1.0 + chi * fv1);

        FloatType dSq = d * d;
        FloatType kappaSq_dSq = _kappa * _kappa * dSq;

        FloatType Sbar = (nuHat / kappaSq_dSq) * fv2;
        FloatType Shat;
        if (Sbar >= -_c2 * S_vort) {
            Shat = S_vort + Sbar;
        } else {
            Shat = S_vort + (S_vort * (_c2 * _c2 * S_vort + _c3 * Sbar)) / ((_c3 - 2.0 * _c2) * S_vort - Sbar);
        }

        FloatType r;
        if (Shat * kappaSq_dSq <= 0.0) {
            r = 10.0;
        } else {
            r = std::min(nuHat / (Shat * kappaSq_dSq), (FloatType)10.0);
        }

        FloatType r2 = r * r;
        FloatType r6 = r2 * r2 * r2;
        FloatType g = r + _cw2 * (r6 - r);
        FloatType g6 = std::pow(g, 6.0);
        FloatType cw3_6 = std::pow(_cw3, 6.0);
        FloatType fw = g * std::pow((1.0 + cw3_6) / (g6 + cw3_6), 1.0 / 6.0);

        FloatType ft2 = _ct3 * std::exp(-_ct4 * chi * chi);

        FloatType P = _cb1 * (1.0 - ft2) * Shat * nuHat;
        FloatType D = (_cw1 * fw - (_cb1 / (_kappa * _kappa)) * ft2) * (nuHat * nuHat / dSq);

        FloatType Qt = P - D + crossDiff;
        return Qt;
    } else {
        _fv1(i, j, k) = 0.0;

        FloatType dSq = d * d;
        FloatType P_neg = _cb1 * (1.0 - _ct3) * S_vort * nuHat;
        FloatType D_neg = _cw1 * (nuHat * nuHat / dSq);

        FloatType Qt = P_neg + D_neg + crossDiff;
        return Qt;
    }
}

void TurbulenceModelSA::updateSolution(const Matrix3D<FloatType> &residual, const Matrix3D<FloatType> &dt){
    Matrix3D<FloatType> volume = _mesh.getVolumes();
    _nuHat -= residual * dt / volume;
}

FloatType TurbulenceModelSA::getEddyViscosity(const FloatType &density, size_t i, size_t j, size_t k) const {
    FloatType nuHat = _nuHat(i, j, k);
    if (nuHat <= 0.0) {
        return 0.0;
    }
    FloatType nuL = _nuLaminar(i, j, k);
    if (nuL <= 0.0) {
        return 0.0;
    }
    FloatType chi = nuHat / nuL;
    FloatType chi3 = chi * chi * chi;
    FloatType cv1_3 = _cv1 * _cv1 * _cv1;
    FloatType fv1 = chi3 / (chi3 + cv1_3);
    FloatType mu = fv1 * density * nuHat;
    return mu;
}