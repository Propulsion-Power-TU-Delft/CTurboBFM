#include "solver_base.hpp"
#include <iostream>

SolverBase::SolverBase(Config& config, Mesh& mesh)
    : _config(config), _mesh(mesh)
{
    _nDimensions = _mesh.getNumberDimensions();
    _nPointsI = _mesh.getNumberPointsI();
    _nPointsJ = _mesh.getNumberPointsJ();
    _nPointsK = _mesh.getNumberPointsK();

    _timeStep.resize(_nPointsI, _nPointsJ, _nPointsK);
    _time.push_back(0.0);

    _topology = _config.getTopology();

    _residualsDropConvergence = _config.getResidualsDropConvergence();  

    _fluidModel = _config.getFluidModel();
    if (_fluidModel == FluidModel::IDEAL){
        _fluid = std::make_unique<FluidIdeal>(_config.getFluidGamma(), _config.getFluidGasConstant());
    }
    else if (_fluidModel == FluidModel::REAL){
        throw std::runtime_error("Real fluid model not implemented yet.");
    }
    else{
        throw std::runtime_error("Unsupported fluid model selected.");
    }

    AdvectionScheme advectionScheme = _config.getAdvectionScheme();
    switch (advectionScheme)
    {
    case AdvectionScheme::JST:
        _advection = std::make_unique<AdvectionJst>(_config, *_fluid);
        break;
    case AdvectionScheme::ROE:
        _advection = std::make_unique<AdvectionRoe>(_config, *_fluid);
        break;
    default:
        throw std::runtime_error("Unsupported convection scheme selected.");
    }

    readBoundaryConditions();
    computeWallDistance();
    
    
}

void SolverBase::readBoundaryConditions(){

    std::array<BoundaryIndex, 6> bounds = {
        BoundaryIndex::I_START,
        BoundaryIndex::I_END,
        BoundaryIndex::J_START,
        BoundaryIndex::J_END,
        BoundaryIndex::K_START,
        BoundaryIndex::K_END};
    
    // read the boundaries type
    for (auto& bound : bounds) {
        _boundaryTypes[bound] = _config.getBoundaryType(bound);
    }

    // read the boundary conditions values
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::INLET || _boundaryTypes[bound] == BoundaryType::INLET_SUPERSONIC) {
            _boundaryValues[bound] = _config.getInletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_2D){
            _boundaryValues[bound] = std::vector<FloatType> {};
            _inlet2DfilePath = _config.getInlet2DfilePath();
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET 
                 || _boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC 
                 || _boundaryTypes[bound] == BoundaryType::THROTTLE){

            _boundaryValues[bound] = _config.getOutletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryValues[bound] = _config.getOutletBCValues();
            _hubStaticPressure = _boundaryValues[bound][0];
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _mesh.checkPeriodicity();
            _boundaryValues[bound] = _config.getPeriodicityInfo();
        }
        else if (_boundaryTypes[bound] == BoundaryType::NO_SLIP_WALL){
            _boundaryVelocities[bound] = _config.getNoSlipWallVelocity(bound);
        }
        else {
            // others type of boundaries don't need any other information -> zero length vector
            _boundaryValues[bound] = std::vector<FloatType> {};
        }
    }

    // instantiate boundary objects
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::INVISCID_WALL || _boundaryTypes[bound] == BoundaryType::NO_SLIP_WALL){
            _boundaryConditions[bound] = std::make_unique<BoundaryInviscidWall>(
                _config, 
                _mesh, 
                *_fluid, 
                bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET){
            _boundaryConditions[bound] = std::make_unique<BoundaryInlet>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_2D){
            _boundaryConditions[bound] = std::make_unique<BoundaryInlet2D>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _inlet2DfilePath);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryInletSupersonic>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutlet>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletRadialEquilibrium>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _radialProfilePressure);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletSupersonic>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::THROTTLE){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletThrottle>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                _radialProfilePressure);
        }
        else if (_boundaryTypes[bound] == BoundaryType::WEDGE){
            _boundaryConditions[bound] = std::make_unique<BoundaryFake>(
                _config, 
                _mesh, 
                *_fluid, 
                bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryFake>(
                _config, 
                _mesh, 
                *_fluid, 
                bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::TRANSPARENT){
            _boundaryConditions[bound] = std::make_unique<BoundaryTransparent>(
                _config, 
                _mesh, 
                *_fluid, 
                bound, 
                *_advection);
        }
    }
}


const std::array<int, 3> SolverBase::getStepMask(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return {1, 0, 0};
        case FluxDirection::J: return {0, 1, 0};
        case FluxDirection::K: return {0, 0, 1};
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}


void SolverBase::getBoundarySliceIndices(
    BoundaryIndex boundaryIdx, 
    size_t &iStart, 
    size_t &iLast, 
    size_t &jStart, 
    size_t &jLast, 
    size_t &kStart, 
    size_t &kLast) const{

    iStart=0, 
    iLast=_nPointsI, 
    jStart=0, 
    jLast=_nPointsJ, 
    kStart=0, 
    kLast=_nPointsK;
    
    switch (boundaryIdx)
    {
    case BoundaryIndex::I_START:
        iLast = 1;
        break;
    case BoundaryIndex::I_END:
        iStart = _nPointsI-1;
        break;
    case BoundaryIndex::J_START:
        jLast = 1;
        break;
    case BoundaryIndex::J_END:
        jStart = _nPointsJ-1;
        break;
    case BoundaryIndex::K_START:
        kLast = 1;
        break;
    case BoundaryIndex::K_END:
        kStart = _nPointsK-1;
        break;
    }
}


void SolverBase::computeWallDistance() {
    _wallDistance.resize(_nPointsI, _nPointsJ, _nPointsK);
    for (size_t i = 0; i < _nPointsI; ++i) {
        for (size_t j = 0; j < _nPointsJ; ++j) {
            for (size_t k = 0; k < _nPointsK; ++k) {
                FloatType minDistance = 1.0E9;
                FloatType distanceTmp = 1.0E9;
                for (auto& bound : _boundaryTypes) {
                    if (bound.second == BoundaryType::INVISCID_WALL || bound.second == BoundaryType::NO_SLIP_WALL){ 
                        distanceTmp = computeMinimumDistanceToBoundary(i, j, k, bound.first);
                    }
                    if (distanceTmp < minDistance){
                        minDistance = distanceTmp;
                    }
                }
                _wallDistance(i,j,k) = minDistance;
            }
        }
    }
}


FloatType SolverBase::computeMinimumDistanceToBoundary(size_t i, size_t j, size_t k, BoundaryIndex boundaryIdx) const {
    size_t iStart, iLast, jStart, jLast, kStart, kLast;
    getBoundarySliceIndices(boundaryIdx, iStart, iLast, jStart, jLast, kStart, kLast);
    
    FloatType dx, dy, dz;
    FloatType minDistance = 1.0E9;
    for (size_t ib = iStart; ib < iLast; ++ib) {
        for (size_t jb = jStart; jb < jLast; ++jb) {
            for (size_t kb = kStart; kb < kLast; ++kb) {
                dx = _mesh.getVertex(i,j,k).x() - _mesh.getVertex(ib,jb,kb).x();
                dy = _mesh.getVertex(i,j,k).y() - _mesh.getVertex(ib,jb,kb).y();
                dz = _mesh.getVertex(i,j,k).z() - _mesh.getVertex(ib,jb,kb).z();
                FloatType distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (distance < minDistance){
                    minDistance = distance;
                }
            }
        }
    }
    return minDistance;
}