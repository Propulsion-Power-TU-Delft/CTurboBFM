#include "turbulence_model_sa.hpp"

TurbulenceModelSA::TurbulenceModelSA(
    const Config &config, 
    const FluidBase &fluid, 
    const Mesh &mesh, 
    const std::vector<Boundary> &boundaries,
    const Matrix3D<FloatType> &wallDistance) 
    : TurbulenceModelBase(config, fluid, mesh, boundaries), 
    _wallDistance(wallDistance) {

        setupModelEquations();
        setupInitialValues();
    }

void TurbulenceModelSA::setupModelEquations() {
    _nuHat.resize(_ni, _nj, _nk);
}

void TurbulenceModelSA::setupInitialValues() {
    if (_config.restartSolution()) {
        throw std::runtime_error("Restart is not supported yet for SA turbulence model.");
    }

    FloatType initTemperature = _config.getInitTemperature();
    FloatType initPressure = _config.getInitPressure();
    FloatType initDensity = _fluid.computeDensity_p_T(initPressure, initTemperature);
    FloatType initNu = _fluid.computeMolecularDynamicViscosity(initTemperature) / initDensity;

    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                _nuHat(i, j, k) = 0.1 * initNu;
            }
        }
    }
}

void TurbulenceModelSA::setupBoundaryValues(
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) {
    
    for (auto &boundary : _boundaries) {

        // zero nu hat on the no-slip walls
        if (boundary.type == BoundaryType::NO_SLIP_WALL) {
            
            enforceNoSlipWallCondition(boundary, sol, solutionGrad);
        }
        else if (
            boundary.type == BoundaryType::INLET || 
            boundary.type == BoundaryType::INLET_2D || 
            boundary.type == BoundaryType::INLET_SUPERSONIC){
            
            enforceInletCondition(boundary, sol, solutionGrad);
        }
        else if (boundary.type == BoundaryType::OUTLET || 
            boundary.type == BoundaryType::RADIAL_EQUILIBRIUM || 
            boundary.type == BoundaryType::THROTTLE ||
            boundary.type == BoundaryType::OUTLET_SUPERSONIC){ 

            enforceOutletCondition(boundary, sol, solutionGrad);
        }
    }
}

void TurbulenceModelSA::enforceNoSlipWallCondition(
    const Boundary& boundary, 
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = 0.0;
            }
        }
    }
}

void TurbulenceModelSA::enforceInletCondition(
    const Boundary& boundary, 
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = 0.0;
            }
        }
    }
}

void TurbulenceModelSA::enforceOutletCondition(
    const Boundary& boundary, 
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) {

    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _ni, _nj, _nk);
    for (size_t i = range.iStart; i < range.iLast; ++i) {
        for (size_t j = range.jStart; j < range.jLast; ++j) {
            for (size_t k = range.kStart; k < range.kLast; ++k) {
                _nuHat(i, j, k) = 0.0;
            }
        }
    }
}

void TurbulenceModelSA::solve(
    const FlowSolution &sol, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad){
        
    std::cout << "TurbulenceModelSA::solve" << std::endl;
    }