#include "solver_base.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

SolverBase::SolverBase(Config& config, Mesh& mesh)
    : _config(config), _mesh(mesh)
{
    setupSolverInfo();
    buildFluidModel();
    buildAdvectionModel();
    readBoundaryConditions();
    computeWallDistance();
}

void SolverBase::setupSolverInfo() {
    _nDimensions = _mesh.getNumberDimensions();
    _nPointsI = _mesh.getNumberPointsI();
    _nPointsJ = _mesh.getNumberPointsJ();
    _nPointsK = _mesh.getNumberPointsK();

    _timeStep.resize(_nPointsI, _nPointsJ, _nPointsK);
    _time.push_back(0.0);

    _topology = _config.getTopology();

    _residualsDropConvergence = _config.getResidualsDropConvergence();  
}

void SolverBase::buildFluidModel() {
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
}

void SolverBase::buildAdvectionModel() {
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
}

void SolverBase::readBoundaryConditions(){
    readBoundaryFile();
    buildBoundaryDataStructures();
    buildBoundaryFluxes();
    buildBoundaryConditionsMap();
}

void SolverBase::readBoundaryFile() {
    std::string filename = _config.getBoundaryConditionsFilePath();
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    bool inDataSection = false;

    while (std::getline(file, line)) {

        // Skip empty lines
        if (line.empty())
            continue;

        // Skip metadata lines (NDIMENSIONS=2, NI=177, etc.)
        if (line.find('=') != std::string::npos)
            continue;

        // Detect CSV header — marks start of data section
        if (line.find("PATCH_NAME") != std::string::npos) {
            inDataSection = true;
            continue;
        }

        if (!inDataSection)
            continue;

        // Parse patch data row
        std::stringstream ss(line);
        std::string token;
        Boundary patch;

        // PATCH_NAME
        std::getline(ss, token, ',');
        token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
        patch.name = token;

        // I_MIN, I_MAX, J_MIN, J_MAX, K_MIN, K_MAX
        auto readInt = [&]() {
            std::getline(ss, token, ',');
            token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
            return std::stoi(token);
        };

        patch.i_min = readInt();
        patch.i_max = readInt();
        patch.j_min = readInt();
        patch.j_max = readInt();
        patch.k_min = readInt();
        patch.k_max = readInt();

        patch.type   = _config.getBoundaryType(patch.name);
        patch.values = _config.getBoundaryValues(patch.name);

        _boundaries.push_back(patch);
    }

    std::cout << "Boundary conditions read from file: " << filename << std::endl;
}

void SolverBase::buildBoundaryDataStructures() {
    for (auto& bound : _boundaries) {
        if (bound.type == BoundaryType::RADIAL_EQUILIBRIUM || bound.type == BoundaryType::THROTTLE){
            RadialEquilibriumProfile profile;
            profile.boundary = bound;

            if(bound.i_min != bound.i_max){
                throw std::runtime_error("Radial equilibrium only supported on j-k patches");
            };

            if(bound.i_min == 0){
                throw std::runtime_error("Radial equilibrium only supported to the last i-position");
            };

            size_t nPoints = bound.j_max - bound.j_min;
            if (nPoints <= 0) {
                throw std::runtime_error("Invalid boundary indices for radial equilibrium patch: " + bound.name);
            }
            
            for (size_t j=bound.j_min; j<bound.j_max; j++){
                FloatType radius = _mesh.getRadius(bound.i_max-1, j, 0);
                profile.radius.push_back(radius);
            }
            profile.pressure.resize(nPoints);
            _radialEquilibriumProfiles.push_back(profile);
        }
        else if (bound.type == BoundaryType::PERIODIC){
            FloatType boundaryId = bound.values[0];

            if (boundaryId == 0){
                _periodicityTranslation = bound.values[1];
                _periodicityAngleDeg = bound.values[2];
                if (_periodicityAngleDeg == 360.0){
                    _periodicityAngleDeg = 0.0;
                }
                _periodicityAngleRad = _periodicityAngleDeg * M_PI / 180.0;
                _mesh.checkPeriodicity(_periodicityTranslation, _periodicityAngleRad);
            }
            
            // some restrictions for now, to be relaxed in the future
            if (bound.k_min != bound.k_max){
                throw std::runtime_error("Periodic boundary only supported on i-j patches");
            };
            if (bound.i_min != 0 || bound.j_min != 0 || bound.i_max != _nPointsI || bound.j_max != _nPointsJ){
                throw std::runtime_error("Periodic boundary only supported on full i-j patches");
            };
            if (bound.k_min == 0 && boundaryId != 0){
                throw std::runtime_error("Periodic boundary with i_min=0 must have boundary id 0");
            }

        }
    }
}

void SolverBase::buildBoundaryFluxes() {
    for (auto& bound : _boundaries) {
        if (bound.type == BoundaryType::INVISCID_WALL || bound.type == BoundaryType::NO_SLIP_WALL){
            bound.fluxMethod = std::make_unique<BoundaryInviscidWall>(
                _config, 
                _mesh, 
                *_fluid);
        }
        else if (bound.type == BoundaryType::INLET){
            bound.fluxMethod = std::make_unique<BoundaryInlet>(
                _config, 
                _mesh, 
                *_fluid, 
                bound.values);
        }
        else if (bound.type == BoundaryType::INLET_2D){
            bound.fluxMethod = std::make_unique<BoundaryInlet2D>(
                _config, 
                _mesh, 
                *_fluid, 
                _inlet2DfilePath);
        }
        else if (bound.type == BoundaryType::INLET_SUPERSONIC){
            bound.fluxMethod = std::make_unique<BoundaryInletSupersonic>(
                _config, 
                _mesh, 
                *_fluid, 
                bound.values);
        }
        else if (bound.type == BoundaryType::OUTLET){
            bound.fluxMethod = std::make_unique<BoundaryOutlet>(
                _config, 
                _mesh, 
                *_fluid, 
                bound.values);
        }
        else if (bound.type == BoundaryType::RADIAL_EQUILIBRIUM){
            bound.fluxMethod = std::make_unique<BoundaryOutletRadialEquilibrium>(
                _config, 
                _mesh, 
                *_fluid, 
                _radialEquilibriumProfiles.back().pressure);
        }
        else if (bound.type == BoundaryType::OUTLET_SUPERSONIC){
            bound.fluxMethod = std::make_unique<BoundaryOutletSupersonic>(
                _config, 
                _mesh, 
                *_fluid, 
                bound.values);
        }
        else if (bound.type == BoundaryType::THROTTLE){
            bound.fluxMethod = std::make_unique<BoundaryOutletThrottle>(
                _config, 
                _mesh, 
                *_fluid, 
                _radialEquilibriumProfiles.back().pressure);
        }
        else if (bound.type == BoundaryType::WEDGE){
            bound.fluxMethod = std::make_unique<BoundaryFake>(
                _config, 
                _mesh, 
                *_fluid);
        }
        else if (bound.type == BoundaryType::PERIODIC){
            bound.fluxMethod = std::make_unique<BoundaryFake>(
                _config, 
                _mesh, 
                *_fluid);
        }
        else if (bound.type == BoundaryType::TRANSPARENT){
            bound.fluxMethod = std::make_unique<BoundaryTransparent>(
                _config, 
                _mesh, 
                *_fluid, 
                *_advection);
        }
    }
}


void SolverBase::buildBoundaryConditionsMap() {
    // build the structure for the objects referencing boundary conditions
    size_t niFaces, njFaces, nkFaces;

    njFaces = _mesh.getSurfacesI().sizeJ();
    nkFaces = _mesh.getSurfacesI().sizeK();
    _boundaryConditionsMapI.resize(2, njFaces, nkFaces);

    nkFaces = _mesh.getSurfacesJ().sizeK();
    niFaces = _mesh.getSurfacesJ().sizeI();
    _boundaryConditionsMapJ.resize(niFaces, 2, nkFaces);

    niFaces = _mesh.getSurfacesK().sizeI();
    njFaces = _mesh.getSurfacesK().sizeJ();
    _boundaryConditionsMapK.resize(niFaces, njFaces, 2);

    // now associate every boundary face to the flux method pointer stored in _boundaries
    for (const auto& bound : _boundaries) {
        size_t i_min = bound.i_min;
        size_t i_max = bound.i_max;
        size_t j_min = bound.j_min;
        size_t j_max = bound.j_max;
        size_t k_min = bound.k_min;
        size_t k_max = bound.k_max;
        
        // Translate node numbering of bound object to dual node numbering of boundary conditions map 
        size_t itmp, jtmp, ktmp;
        if (i_min == i_max) { 
            if (i_min == 0) {
                itmp = 0;
            }
            else {
                itmp = 1;
            }
            for (size_t j = j_min; j < j_max; ++j) {
                for (size_t k = k_min; k < k_max; ++k) {
                    _boundaryConditionsMapI(itmp, j, k) = bound.fluxMethod;
                }
            }
        } 
        else if (j_min == j_max) { // J boundaries
            if (j_min == 0) {
                jtmp = 0;
            }
            else {
                jtmp = 1;
            }
            for (size_t i = i_min; i < i_max; ++i) {
                for (size_t k = k_min; k < k_max; ++k) {
                    _boundaryConditionsMapJ(i, jtmp, k) = bound.fluxMethod;
                }
            }
        }
        else if (k_min == k_max) { // K boundaries
            if (k_min == 0) {
                ktmp = 0;
            }
            else {
                ktmp = 1;
            }
            for (size_t i = i_min; i < i_max; ++i) {
                for (size_t j = j_min; j < j_max; ++j) {
                    _boundaryConditionsMapK(i, j, ktmp) = bound.fluxMethod;
                }
            }
        }
        else {
            throw std::runtime_error("Invalid boundary definition: " + bound.name);
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

const Matrix3D<std::shared_ptr<BoundaryBase>>& SolverBase::getBoundaryConditionsMap(FluxDirection direction) const {
    switch (direction)
    {
    case FluxDirection::I:
        return _boundaryConditionsMapI;
        break;
    case FluxDirection::J:
        return _boundaryConditionsMapJ;
        break;
    case FluxDirection::K:
        return _boundaryConditionsMapK;
        break;
    default:
        throw std::runtime_error("Invalid flux direction.");
    }
}

void SolverBase::computeWallDistance() {
    _wallDistance.resize(_nPointsI, _nPointsJ, _nPointsK);
    for (size_t i = 0; i < _nPointsI; ++i) {
        for (size_t j = 0; j < _nPointsJ; ++j) {
            for (size_t k = 0; k < _nPointsK; ++k) {
                FloatType minDistance = 1.0E9;
                FloatType distanceTmp = 1.0E9;
                for (auto& bound : _boundaries) {
                    if (bound.type == BoundaryType::NO_SLIP_WALL){ 
                        distanceTmp = computeMinimumDistanceToBoundary(i, j, k, bound);
                        if (distanceTmp < minDistance){
                            minDistance = distanceTmp;
                        }
                    }
                }

                _wallDistance(i,j,k) = minDistance;
            }
        }
    }
}

FloatType SolverBase::computeMinimumDistanceToBoundary(size_t i, size_t j, size_t k, Boundary boundary) const {
    
    BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(boundary, _nPointsI, _nPointsJ, _nPointsK);
    
    FloatType dx, dy, dz;
    FloatType minDistance = 1.0E9;
    for (size_t ib = range.iStart; ib < range.iLast; ++ib) {
        for (size_t jb = range.jStart; jb < range.jLast; ++jb) {
            for (size_t kb = range.kStart; kb < range.kLast; ++kb) {
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