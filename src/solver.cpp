#include "solver.hpp"
#include "kdtree.hpp"
#include "math_utils.hpp"
#include "types.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include <cstdlib>
#include <sstream>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

Solver::Solver(Config& config, Mesh& mesh)
    : _config(config), _mesh(mesh)
{
    setupSolverInfo();
    buildFluidModel();
    buildAdvectionModel();
    readBoundaryConditions();
    computeWallDistance();
    buildBfmInfo();
    initializeSolutionArrays();
    buildTurbulenceModel();
    buildOutputStructure();
}

void Solver::setupSolverInfo() {
    _nDimensions = _mesh.getNumberDimensions();
    _nPointsI = _mesh.getNumberPointsI();
    _nPointsJ = _mesh.getNumberPointsJ();
    _nPointsK = _mesh.getNumberPointsK();

    _timeStep.resize(_nPointsI, _nPointsJ, _nPointsK);
    _time.push_back(0.0);
    _currentTime = 0.0;
    _historyBufferSize = _config.getHistoryBufferSize();

    _topology = _config.getTopology();

    _residualsDropConvergence = _config.getResidualsDropConvergence();  
}

void Solver::buildFluidModel() {
    _fluidModel = _config.getFluidModel();
    if (_fluidModel == FluidModel::IDEAL){
        _fluid = std::make_unique<FluidIdeal>(_config.getFluidGamma(), _config.getFluidGasConstant());
    }
    else if (_fluidModel == FluidModel::REAL){
        std::string tableFile = _config.getFluidTableFile();
        bool forceRegen = _config.isFluidTableForceRegenerate();
        bool tableExists = std::filesystem::exists(tableFile);

        if (!tableExists || forceRegen) {
            std::cout << "Fluid table '" << tableFile << "' not found or regeneration requested. Auto-generating..." << std::endl;
            generateFluidTable(tableFile);
        }

        std::cout << "Loading fluid look-up table: " << tableFile << std::endl;
        _fluid = std::make_unique<FluidReal>(tableFile);
    }
    else{
        throw std::runtime_error("Unsupported fluid model selected.");
    }

    if (_config.isViscosityActive()){
        _fluid->setTransportProperties(_config);
    }
}

void Solver::generateFluidTable(const std::string& tableFile) {
    // 1. Locate python generator script
    std::filesystem::path scriptPath;
    std::string configScript = _config.getFluidTableGeneratorScript();
    if (!configScript.empty() && std::filesystem::exists(configScript)) {
        scriptPath = configScript;
    } else {
        std::vector<std::filesystem::path> candidates = {
            "python/generate_fluid_table.py",
            "../python/generate_fluid_table.py",
            "../../python/generate_fluid_table.py",
            "../../../python/generate_fluid_table.py"
        };
        if (const char* root = std::getenv("TURBOBFM_ROOT")) {
            candidates.push_back(std::filesystem::path(root) / "python" / "generate_fluid_table.py");
        }
#if defined(__APPLE__)
        char buf[2048];
        uint32_t bufSize = sizeof(buf);
        if (_NSGetExecutablePath(buf, &bufSize) == 0) {
            std::filesystem::path exePath = std::filesystem::weakly_canonical(buf);
            candidates.push_back(exePath.parent_path().parent_path() / "python" / "generate_fluid_table.py");
        }
#elif defined(__linux__)
        if (std::filesystem::exists("/proc/self/exe")) {
            std::filesystem::path exePath = std::filesystem::weakly_canonical("/proc/self/exe");
            candidates.push_back(exePath.parent_path().parent_path() / "python" / "generate_fluid_table.py");
        }
#endif
        for (const auto& cand : candidates) {
            if (std::filesystem::exists(cand)) {
                scriptPath = std::filesystem::absolute(cand);
                break;
            }
        }
    }

    if (scriptPath.empty() || !std::filesystem::exists(scriptPath)) {
        throw std::runtime_error("Could not find table generator script 'python/generate_fluid_table.py'. "
                                 "Please set TURBOBFM_ROOT or FLUID_TABLE_GENERATOR in configuration.");
    }

    std::string pythonExe = _config.getFluidTablePythonExecutable();
    std::string fluidName = _config.getFluidName();

    FloatType pMin = 0.0, pMax = 0.0, tMin = 0.0, tMax = 0.0;
    if (_config.hasFluidTablePMin() && _config.hasFluidTablePMax() &&
        _config.hasFluidTableTMin() && _config.hasFluidTableTMax()) {
        pMin = _config.getFluidTablePMin();
        pMax = _config.getFluidTablePMax();
        tMin = _config.getFluidTableTMin();
        tMax = _config.getFluidTableTMax();
    } else if (_config.has("INIT_PRESSURE") && _config.has("INIT_TEMPERATURE")) {
        FloatType pInit = _config.getInitPressure();
        FloatType tInit = _config.getInitTemperature();
        pMin = _config.hasFluidTablePMin() ? _config.getFluidTablePMin() : std::max(1.0e4, pInit * 0.5);
        pMax = _config.hasFluidTablePMax() ? _config.getFluidTablePMax() : pInit * 1.5;
        tMin = _config.hasFluidTableTMin() ? _config.getFluidTableTMin() : std::max(100.0, tInit * 0.85);
        tMax = _config.hasFluidTableTMax() ? _config.getFluidTableTMax() : tInit * 1.35;
    } else {
        pMin = _config.hasFluidTablePMin() ? _config.getFluidTablePMin() : 5.0e6;
        pMax = _config.hasFluidTablePMax() ? _config.getFluidTablePMax() : 1.2e7;
        tMin = _config.hasFluidTableTMin() ? _config.getFluidTableTMin() : 305.0;
        tMax = _config.hasFluidTableTMax() ? _config.getFluidTableTMax() : 450.0;
    }

    int nRho = _config.getFluidTableNRho();
    int nE = _config.getFluidTableNE();
    int nP = _config.getFluidTableNP();
    int nT = _config.getFluidTableNT();
    int nS = _config.getFluidTableNS();

    // Ensure parent directory of tableFile exists
    std::filesystem::path tablePath(tableFile);
    if (tablePath.has_parent_path()) {
        std::filesystem::create_directories(tablePath.parent_path());
    }

    std::ostringstream cmd;
    cmd << pythonExe << " \"" << scriptPath.string() << "\""
        << " --fluid \"" << fluidName << "\""
        << " --output \"" << tablePath.string() << "\""
        << " --p_min " << pMin
        << " --p_max " << pMax
        << " --T_min " << tMin
        << " --T_max " << tMax
        << " --n_rho " << nRho
        << " --n_e " << nE
        << " --n_p " << nP
        << " --n_T " << nT
        << " --n_s " << nS;

    std::cout << "[FluidReal] Auto-generating thermodynamic table for '" << fluidName << "'..." << std::endl;
    std::cout << "  Script:  " << scriptPath.string() << std::endl;
    std::cout << "  Output:  " << tablePath.string() << std::endl;
    std::cout << "  Ranges:  P in [" << pMin << ", " << pMax << "] Pa, T in [" << tMin << ", " << tMax << "] K" << std::endl;
    std::cout << "  Points:  (n_rho=" << nRho << ", n_e=" << nE << ", n_p=" << nP << ", n_t=" << nT << ", n_s=" << nS << ")" << std::endl;

    int ret = std::system(cmd.str().c_str());
    if (ret != 0 || !std::filesystem::exists(tablePath)) {
        throw std::runtime_error("[FluidReal] Failed to generate table file '" + tablePath.string() + "'. Command exit code: " + std::to_string(ret));
    }
    std::cout << "[FluidReal] Thermodynamic table generated successfully." << std::endl;
}

void Solver::buildAdvectionModel() {
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

void Solver::readBoundaryConditions(){
    readBoundaryFile();
    buildBoundaryDataStructures();
    buildBoundaryFluxes();
    buildBoundaryConditionsMap();
    checkBoundaryConditionsMap();
}

void Solver::readBoundaryFile() {
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
        patch.computeOrientation();

        _boundaries.push_back(patch);
    }

    std::cout << "Boundary conditions read from file: " << filename << std::endl;
}

void Solver::buildBoundaryDataStructures() {
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

            if (bound.j_max <= bound.j_min) {
                throw std::runtime_error("Invalid boundary indices for radial equilibrium patch: " + bound.name);
            }
            size_t nPoints = bound.j_max - bound.j_min;
            
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

void Solver::buildBoundaryFluxes() {
    size_t radialProfileIdx = 0;
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
                _radialEquilibriumProfiles[radialProfileIdx++].pressure);
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
                _radialEquilibriumProfiles[radialProfileIdx++].pressure);
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
        else if (bound.type == BoundaryType::FARFIELD){
            bound.fluxMethod = std::make_unique<BoundaryFarfield>(
                _config,
                _mesh,
                *_fluid,
                bound.values);
        }
    }
}


void Solver::buildBoundaryConditionsMap() {
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

void Solver::checkBoundaryConditionsMap() const {

    std::vector<Matrix3D<std::shared_ptr<BoundaryBase>>> boundaryConditionsMaps;
    boundaryConditionsMaps.push_back(_boundaryConditionsMapI);
    
    if (_config.getTopology() != Topology::ONE_DIMENSIONAL) {
        boundaryConditionsMaps.push_back(_boundaryConditionsMapJ);
    }

    if (_config.getTopology() == Topology::THREE_DIMENSIONAL) {
        boundaryConditionsMaps.push_back(_boundaryConditionsMapK);
    }
    
    // check that at least a boundary flux method is associated to every boundary node
    for (const auto& bcMap : boundaryConditionsMaps) {
        for (size_t i = 0; i < bcMap.sizeI(); ++i) {
            for (size_t j = 0; j < bcMap.sizeJ(); ++j) {
                for (size_t k = 0; k < bcMap.sizeK(); ++k) {
                    if (bcMap(i, j, k) == nullptr) {
                        throw std::runtime_error("Missing boundary condition for i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k));
                    }
                }
            }
        }
    }
}

const std::array<int, 3> Solver::getStepMask(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return {1, 0, 0};
        case FluxDirection::J: return {0, 1, 0};
        case FluxDirection::K: return {0, 0, 1};
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}


void Solver::getBoundarySliceIndices(
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

const Matrix3D<std::shared_ptr<BoundaryBase>>& Solver::getBoundaryConditionsMap(FluxDirection direction) const {
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

void Solver::computeWallDistance() {
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

FloatType Solver::computeMinimumDistanceToBoundary(size_t i, size_t j, size_t k, Boundary boundary) const {
    
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

void Solver::buildBfmInfo(){
    BodyForceModel bfmModel = _config.getBFMModel();
    if (bfmModel == BodyForceModel::HALL) {
        _bfmSource = std::make_unique<SourceBFMHall>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BodyForceModel::HALL_THOLLET) {
        _bfmSource = std::make_unique<SourceBFMThollet>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BodyForceModel::CHIMA) {
        _bfmSource = std::make_unique<SourceBFMChima>(_config, *_fluid, _mesh, _turboPerformance);
    }
    else if (bfmModel == BodyForceModel::ONLY_BLOCKAGE) {
        _bfmSource = std::make_unique<SourceBFMBase>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BodyForceModel::CORRELATIONS){
        _bfmSource = std::make_unique<SourceBFMCorrelations>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BodyForceModel::LIFT_DRAG){
        _bfmSource = std::make_unique<SourceBFMLiftDrag>(_config, *_fluid, _mesh);
    }
    else if (bfmModel == BodyForceModel::NONE) {
        _bfmSource = std::make_unique<SourceBFMBase>(_config, *_fluid, _mesh);
    }
    else {
        throw std::runtime_error("Unsupported BFM model selected.");
    }
    
    _isBfmActive = _config.isBFMActive();
    
    _isGongFormulationActive = _config.isGongFormulationActive();

    _isGreitzerModelingActive = _config.enableGreitzerModeling();
    if (_isGreitzerModelingActive) {
        FloatType throttleCoeff;
        bool hasThrottleCoeff = false;
        for (auto& bound : _boundaries) {
            if (bound.type == BoundaryType::THROTTLE){ 
                hasThrottleCoeff = true;
                throttleCoeff = bound.values[0];
            }
        }

        if (!hasThrottleCoeff) {
            throw std::runtime_error("Greitzer modeling enabled but no throttle boundary condition found!");
        }

        _greitzerModel = std::make_unique<GreitzerModel>(_config, *_fluid, throttleCoeff);
    }
}

void Solver::buildTurbulenceModel(){
    _isTurbulenceActive = _config.isTurbulenceActive();
    TurbulenceModel turbulenceModel = _config.getTurbulenceModel();
    if (turbulenceModel == TurbulenceModel::SPALART_ALLMARAS) {
        _turbulenceModel = std::make_unique<TurbulenceModelSA>(
            _config,
            *_fluid, 
            _mesh, 
            _boundaries, 
            _wallDistance, 
            _conservativeSolution);
    }
    else if (turbulenceModel == TurbulenceModel::NONE) {
        _turbulenceModel = std::make_unique<TurbulenceModelNone>(
            _config, 
            *_fluid, 
            _mesh, 
            _boundaries, 
            _wallDistance, 
            _conservativeSolution);
    }
    else {
        throw std::runtime_error("Unsupported turbulence model selected.");
    }
}

void Solver::initializeSolutionArrays(){
    _conservativeSolution.resize(_nPointsI, _nPointsJ, _nPointsK);
    _inviscidForce.resize(_nPointsI, _nPointsJ, _nPointsK);
    _viscousForce.resize(_nPointsI, _nPointsJ, _nPointsK);
    _deviationAngle.resize(_nPointsI, _nPointsJ, _nPointsK);

    bool restartSolution = _config.restartSolution();
    if (restartSolution) {
        initializeSolutionFromRestart();
    }
    else {
        initializeSolutionFromScratch();
    }

    _solutionGrad[SolutionName::DENSITY] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
    _solutionGrad[SolutionName::VELOCITY_X] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
    _solutionGrad[SolutionName::VELOCITY_Y] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
    _solutionGrad[SolutionName::VELOCITY_Z] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
    _solutionGrad[SolutionName::TOTAL_ENERGY] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
    _solutionGrad[SolutionName::TEMPERATURE] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);

    computeGradientOfField(_conservativeSolution.getDensity(), _solutionGrad[SolutionName::DENSITY]);
    computeGradientOfField(_conservativeSolution.getVelocityX(), _solutionGrad[SolutionName::VELOCITY_X]);
    computeGradientOfField(_conservativeSolution.getVelocityY(), _solutionGrad[SolutionName::VELOCITY_Y]);
    computeGradientOfField(_conservativeSolution.getVelocityZ(), _solutionGrad[SolutionName::VELOCITY_Z]);
    computeGradientOfField(_conservativeSolution.getTotalEnergy(), _solutionGrad[SolutionName::TOTAL_ENERGY]);

    if (_isGreitzerModelingActive){
        updateMassFlows(_conservativeSolution);
        updateTurboPerformance(_conservativeSolution);
        FloatType massflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
        StateVector conservativeHubOutlet = _conservativeSolution.at(_nPointsI-1, 0, 0);
        StateVector primitiveHubOutlet = getPrimitiveVariablesFromConservative(conservativeHubOutlet);
        FloatType outletPressureHub = _fluid->computePressure_primitive(primitiveHubOutlet);
        FloatType outletTemperatureHub = _fluid->computeTemperature_rho_u_et(
            primitiveHubOutlet[0], 
            {primitiveHubOutlet[1], primitiveHubOutlet[2], primitiveHubOutlet[3]}, 
            primitiveHubOutlet[4]);
        _greitzerModel->initializeState(outletPressureHub, massflow, massflow, outletTemperatureHub);
    }

    for (auto& radialProfile: _radialEquilibriumProfiles) {
        Boundary bound = radialProfile.boundary;
        size_t k = 0;
        size_t i = bound.i_min-1;
        for (size_t j=bound.j_min; j<bound.j_max; j++){
            StateVector conservative = _conservativeSolution.at(i, j, 0);
            StateVector primitive = getPrimitiveVariablesFromConservative(conservative);
            FloatType pressure = _fluid->computePressure_primitive(primitive);
            radialProfile.pressure[k] = pressure;
            k++;
        }
    }

}

void Solver::buildOutputStructure(){
    _output = std::make_unique<Output>(
        _config, 
        _mesh, 
        _conservativeSolution, 
        _solutionGrad,
        *_fluid, 
        _boundaries,
        *_turbulenceModel,
        _inviscidForce, 
        _viscousForce, 
        _deviationAngle,
        _wallDistance);
}

void Solver::initializeSolutionFromScratch(){
    FloatType initMach = _config.getInitMachNumber();
    FloatType initTemperature = _config.getInitTemperature();
    FloatType initPressure = _config.getInitPressure();
    Vector3D initDirection = _config.getInitDirection();

    Matrix3D<Vector3D> flowDirection(_nPointsI, _nPointsJ, _nPointsK);
    if (initDirection == Vector3D{0.0, 0.0, 0.0}) { // alias for adaptive scenario
        _mesh.computeAdaptiveFlowDirection(flowDirection);
    }
    else {
        _mesh.computeUniformFlowDirection(initDirection, flowDirection);    
    }

    FloatType density {0.0}, totEnergy {0.0};
    Vector3D velocity {0.0, 0.0, 0.0};
    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                _fluid->computeInitFields(
                    initMach, 
                    initTemperature, 
                    initPressure, 
                    flowDirection(i,j,k), 
                    density, 
                    velocity, 
                    totEnergy);
                _conservativeSolution._rho(i,j,k) = density;
                _conservativeSolution._rhoU(i,j,k) = density * velocity.x();
                _conservativeSolution._rhoV(i,j,k) = density * velocity.y();
                _conservativeSolution._rhoW(i,j,k) = density * velocity.z();
                _conservativeSolution._rhoE(i,j,k) = density * totEnergy;
            }
        }
    }
}

void Solver::initializeSolutionFromRestart(){
    std::string restartFileName = _config.getRestartFilepath();
    
    size_t NI=0, NJ=0, NK=0;    
    Matrix3D<Vector3D> inputCoordinates;
    Matrix3D<FloatType> inputDensity;
    Matrix3D<FloatType> inputVelX;
    Matrix3D<FloatType> inputVelY;
    Matrix3D<FloatType> inputVelZ;
    Matrix3D<FloatType> inputTemperature;
    Matrix3D<Vector3D> inputForceViscous;
    Matrix3D<Vector3D> inputForceInviscid;
    
    readRestartFile(
        restartFileName, 
        NI, 
        NJ, 
        NK, 
        inputCoordinates,
        inputDensity, 
        inputVelX, 
        inputVelY, 
        inputVelZ, 
        inputTemperature, 
        inputForceViscous, 
        inputForceInviscid);

    std::string restartType = _config.getRestartType();

    if (restartType == "nearest_neighbor" || restartType == "nearest_neighbor_axisymmetric") {
        nearestNeighborRestart(
            inputCoordinates,
            inputDensity, 
            inputVelX, 
            inputVelY, 
            inputVelZ, 
            inputTemperature, 
            inputForceViscous, 
            inputForceInviscid);
    }
    else if (NI != _nPointsI || NJ != _nPointsJ || NK != _nPointsK) {
        if (NI == _nPointsI && NJ == _nPointsJ && restartType == "axisymmetric") {
            axisymmetricRestart(inputDensity, inputVelX, inputVelY, inputVelZ, inputTemperature);
        }
        else if (restartType == "standard") {
            throw std::runtime_error(
                "Restart file dimensions (" + std::to_string(NI) + "," + std::to_string(NJ) + "," + std::to_string(NK) +
                ") do not match solver dimensions (" + std::to_string(_nPointsI) + "," + std::to_string(_nPointsJ) + "," + std::to_string(_nPointsK) + 
                ") and RESTART_TYPE=standard was requested.");
        }
        else {
            std::cout << "Restart file dimensions (" << NI << "," << NJ << "," << NK 
                      << ") do not match solver dimensions (" << _nPointsI << "," << _nPointsJ << "," << _nPointsK 
                      << "). Using nearest-neighbor interpolation.\n";
            nearestNeighborRestart(
                inputCoordinates,
                inputDensity, 
                inputVelX, 
                inputVelY, 
                inputVelZ, 
                inputTemperature, 
                inputForceViscous, 
                inputForceInviscid);
        }
    }
    else {
        standardRestart(
            inputDensity, 
            inputVelX, 
            inputVelY, 
            inputVelZ, 
            inputTemperature, 
            inputForceViscous, 
            inputForceInviscid);
    }
}


void Solver::standardRestart(
    Matrix3D<FloatType> &inputDensity, 
    Matrix3D<FloatType> &inputVelX, 
    Matrix3D<FloatType> &inputVelY, 
    Matrix3D<FloatType> &inputVelZ, 
    Matrix3D<FloatType> &inputTemperature, 
    Matrix3D<Vector3D> &inputForceViscous, 
    Matrix3D<Vector3D> &inputForceInviscid) {

    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                _conservativeSolution._rho(i,j,k) = inputDensity(i,j,k);
                _conservativeSolution._rhoU(i,j,k) = inputDensity(i,j,k) * inputVelX(i,j,k);
                _conservativeSolution._rhoV(i,j,k) = inputDensity(i,j,k) * inputVelY(i,j,k);
                _conservativeSolution._rhoW(i,j,k) = inputDensity(i,j,k) * inputVelZ(i,j,k);

                FloatType pressure = _fluid->computePressure_rho_T(inputDensity(i,j,k), inputTemperature(i,j,k));
                FloatType staticEnergy = _fluid->computeStaticEnergy_p_rho(pressure, inputDensity(i,j,k));
                FloatType totalEnergy = staticEnergy + 0.5*(
                    inputVelX(i,j,k)*inputVelX(i,j,k) 
                    + inputVelY(i,j,k)*inputVelY(i,j,k) 
                    + inputVelZ(i,j,k)*inputVelZ(i,j,k));
                _conservativeSolution._rhoE(i,j,k) = inputDensity(i,j,k) * totalEnergy;
            }
        }
    }
    std::cout << "Standard initialization done.\n";

    _viscousForce = inputForceViscous;    
    _inviscidForce = inputForceInviscid;
}


void Solver::axisymmetricRestart(
    Matrix3D<FloatType> &inputDensity, 
    Matrix3D<FloatType> &inputVelX, 
    Matrix3D<FloatType> &inputVelY, 
    Matrix3D<FloatType> &inputVelZ, 
    Matrix3D<FloatType> &inputTemperature) {

    FloatType thetaInitial, thetaPoint, thetaRotation;
    Vector3D velocityInitial, velocityPoint;
    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                thetaInitial = atan2FromZeroTo2pi(_mesh.getVertex(i,j,0).z(), _mesh.getVertex(i,j,0).y());
                thetaPoint = atan2FromZeroTo2pi(_mesh.getVertex(i,j,k).z(), _mesh.getVertex(i,j,k).y());
                thetaRotation = thetaPoint - thetaInitial;
                
                velocityInitial = Vector3D(inputVelX(i,j,0), inputVelY(i,j,0), inputVelZ(i,j,0));
                velocityPoint = rotateVectorAlongXAxis(velocityInitial, thetaRotation);

                _conservativeSolution._rho(i,j,k) = inputDensity(i,j,0);
                _conservativeSolution._rhoU(i,j,k) = inputDensity(i,j,0) * velocityPoint.x();
                _conservativeSolution._rhoV(i,j,k) = inputDensity(i,j,0) * velocityPoint.y();
                _conservativeSolution._rhoW(i,j,k) = inputDensity(i,j,0) * velocityPoint.z();

                FloatType pressure = _fluid->computePressure_rho_T(inputDensity(i,j,0), inputTemperature(i,j,0));
                FloatType staticEnergy = _fluid->computeStaticEnergy_p_rho(pressure, inputDensity(i,j,0));
                FloatType totalEnergy = staticEnergy + 0.5*(velocityPoint.dot(velocityPoint));
                _conservativeSolution._rhoE(i,j,k) = inputDensity(i,j,0) * totalEnergy;
            }
        }
    }
    std::cout << "Axisymmetric initialization done.\n";
}


void Solver::nearestNeighborRestart(
    const Matrix3D<Vector3D> &inputCoordinates,
    const Matrix3D<FloatType> &inputDensity, 
    const Matrix3D<FloatType> &inputVelX, 
    const Matrix3D<FloatType> &inputVelY, 
    const Matrix3D<FloatType> &inputVelZ, 
    const Matrix3D<FloatType> &inputTemperature, 
    const Matrix3D<Vector3D> &inputForceViscous, 
    const Matrix3D<Vector3D> &inputForceInviscid) {

    size_t NI = inputDensity.sizeI();
    size_t NJ = inputDensity.sizeJ();
    size_t NK = inputDensity.sizeK();
    size_t totalPoints = NI * NJ * NK;

    if (totalPoints == 0) {
        throw std::runtime_error("Cannot perform nearest-neighbor restart: restart data is empty.");
    }

    std::string restartType = _config.getRestartType();
    bool isAxisymmetric = (restartType == "axisymmetric" ||
                           restartType == "nearest_neighbor_axisymmetric" ||
                           _topology == Topology::AXISYMMETRIC ||
                           (NK == 1 && _nPointsK > 1 && _mesh.isPeriodicityActive()));

    if (isAxisymmetric) {
        std::vector<std::array<FloatType, 2>> coarsePointsXR(totalPoints);
        for (size_t idx = 0; idx < totalPoints; ++idx) {
            Vector3D pt = inputCoordinates[idx];
            FloatType r = std::sqrt(pt.y() * pt.y() + pt.z() * pt.z());
            coarsePointsXR[idx] = {pt.x(), r};
        }

        KDTree<2> tree(coarsePointsXR);

        #pragma omp parallel for
        for (int64_t i = 0; i < static_cast<int64_t>(_nPointsI); ++i) {
            for (size_t j = 0; j < _nPointsJ; ++j) {
                for (size_t k = 0; k < _nPointsK; ++k) {
                    Vector3D pt = _mesh.getVertex(i, j, k);
                    FloatType rTarget = std::sqrt(pt.y() * pt.y() + pt.z() * pt.z());
                    FloatType thetaTarget = atan2FromZeroTo2pi(pt.z(), pt.y());

                    size_t bestIdx = tree.findNearest({pt.x(), rTarget});

                    Vector3D ptCoarse = inputCoordinates[bestIdx];
                    FloatType thetaCoarse = atan2FromZeroTo2pi(ptCoarse.z(), ptCoarse.y());
                    FloatType thetaRot = thetaTarget - thetaCoarse;

                    FloatType rho = inputDensity[bestIdx];
                    FloatType T   = inputTemperature[bestIdx];

                    Vector3D velCoarse(inputVelX[bestIdx], inputVelY[bestIdx], inputVelZ[bestIdx]);
                    Vector3D velPoint = rotateVectorAlongXAxis(velCoarse, thetaRot);

                    _conservativeSolution._rho(i, j, k)  = rho;
                    _conservativeSolution._rhoU(i, j, k) = rho * velPoint.x();
                    _conservativeSolution._rhoV(i, j, k) = rho * velPoint.y();
                    _conservativeSolution._rhoW(i, j, k) = rho * velPoint.z();

                    FloatType pressure = _fluid->computePressure_rho_T(rho, T);
                    FloatType staticEnergy = _fluid->computeStaticEnergy_p_rho(pressure, rho);
                    FloatType totalEnergy = staticEnergy + 0.5 * velPoint.dot(velPoint);
                    _conservativeSolution._rhoE(i, j, k) = rho * totalEnergy;

                    if (_isBfmActive) {
                        _inviscidForce(i, j, k) = rotateVectorAlongXAxis(inputForceInviscid[bestIdx], thetaRot);
                        _viscousForce(i, j, k)  = rotateVectorAlongXAxis(inputForceViscous[bestIdx], thetaRot);
                    }
                }
            }
        }
        std::cout << "Axisymmetric nearest-neighbor initialization done.\n";
    }
    else {
        std::vector<std::array<FloatType, 3>> coarsePointsXYZ(totalPoints);
        for (size_t idx = 0; idx < totalPoints; ++idx) {
            Vector3D pt = inputCoordinates[idx];
            coarsePointsXYZ[idx] = {pt.x(), pt.y(), pt.z()};
        }

        KDTree<3> tree(coarsePointsXYZ);

        #pragma omp parallel for
        for (int64_t i = 0; i < static_cast<int64_t>(_nPointsI); ++i) {
            for (size_t j = 0; j < _nPointsJ; ++j) {
                for (size_t k = 0; k < _nPointsK; ++k) {
                    Vector3D pt = _mesh.getVertex(i, j, k);
                    size_t bestIdx = tree.findNearest({pt.x(), pt.y(), pt.z()});

                    FloatType rho = inputDensity[bestIdx];
                    FloatType u   = inputVelX[bestIdx];
                    FloatType v   = inputVelY[bestIdx];
                    FloatType w   = inputVelZ[bestIdx];
                    FloatType T   = inputTemperature[bestIdx];

                    _conservativeSolution._rho(i, j, k)  = rho;
                    _conservativeSolution._rhoU(i, j, k) = rho * u;
                    _conservativeSolution._rhoV(i, j, k) = rho * v;
                    _conservativeSolution._rhoW(i, j, k) = rho * w;

                    FloatType pressure = _fluid->computePressure_rho_T(rho, T);
                    FloatType staticEnergy = _fluid->computeStaticEnergy_p_rho(pressure, rho);
                    FloatType totalEnergy = staticEnergy + 0.5 * (u * u + v * v + w * w);
                    _conservativeSolution._rhoE(i, j, k) = rho * totalEnergy;

                    if (_isBfmActive) {
                        _inviscidForce(i, j, k) = inputForceInviscid[bestIdx];
                        _viscousForce(i, j, k)  = inputForceViscous[bestIdx];
                    }
                }
            }
        }
        std::cout << "Cartesian nearest-neighbor initialization done.\n";
    }
}


void Solver::readRestartFile(
    const std::string &restartFileName, 
    size_t &NI, 
    size_t &NJ, 
    size_t &NK,
    Matrix3D<Vector3D> &inputCoordinates,
    Matrix3D<FloatType> &inputDensity, 
    Matrix3D<FloatType> &inputVelX, 
    Matrix3D<FloatType> &inputVelY, 
    Matrix3D<FloatType> &inputVelZ, 
    Matrix3D<FloatType> &inputTemperature, 
    Matrix3D<Vector3D> &inputForceViscous, 
    Matrix3D<Vector3D> &inputForceInviscid) {

    std::ifstream file(restartFileName);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open restart file: " + restartFileName);
    }

    std::string line;

    // Read NI, NJ, NK
    std::getline(file, line); NI = std::stoi(line.substr(line.find('=') + 1));
    std::getline(file, line); NJ = std::stoi(line.substr(line.find('=') + 1));
    std::getline(file, line); NK = std::stoi(line.substr(line.find('=') + 1));

    inputCoordinates.resize(NI, NJ, NK);
    inputDensity.resize(NI, NJ, NK);
    inputVelX.resize(NI, NJ, NK);
    inputVelY.resize(NI, NJ, NK);
    inputVelZ.resize(NI, NJ, NK);
    inputTemperature.resize(NI, NJ, NK);
    inputForceViscous.resize(NI, NJ, NK);
    inputForceInviscid.resize(NI, NJ, NK);

    bool isBfmSimulation = _config.isBFMActive();

    std::getline(file, line);
    std::istringstream headerStream(line);
    std::string column;
    std::unordered_map<std::string, int> columnIndex;
    int idx = 0;
    while (std::getline(headerStream, column, ',')) {
        size_t first = column.find_first_not_of(" \t\r\n");
        size_t last = column.find_last_not_of(" \t\r\n");
        if (first != std::string::npos) {
            column = column.substr(first, last - first + 1);
        }
        columnIndex[column] = idx++;
    }

    // Optional coordinate columns
    bool hasX = (columnIndex.find("x") != columnIndex.end() || columnIndex.find("X") != columnIndex.end());
    bool hasY = (columnIndex.find("y") != columnIndex.end() || columnIndex.find("Y") != columnIndex.end());
    bool hasZ = (columnIndex.find("z") != columnIndex.end() || columnIndex.find("Z") != columnIndex.end());
    bool hasCoords = (hasX && hasY && hasZ);
    constexpr size_t UNSET = std::numeric_limits<size_t>::max();
    size_t ix = columnIndex.count("x") ? columnIndex["x"] : (columnIndex.count("X") ? columnIndex["X"] : UNSET);
    size_t iy = columnIndex.count("y") ? columnIndex["y"] : (columnIndex.count("Y") ? columnIndex["Y"] : UNSET);
    size_t iz = columnIndex.count("z") ? columnIndex["z"] : (columnIndex.count("Z") ? columnIndex["Z"] : UNSET);

    // Get indexes of the fields of interest (primary fields should always be available in the restart)
    size_t iDensity   = columnIndex.at("Density");
    size_t iVelX      = columnIndex.at("Velocity X");
    size_t iVelY      = columnIndex.at("Velocity Y");
    size_t iVelZ      = columnIndex.at("Velocity Z");
    size_t iTotEnergy = columnIndex.at("Total Energy");

    // these could also not be there, not a problem
    size_t iForceInviscidX{UNSET};
    size_t iForceInviscidY{UNSET};
    size_t iForceInviscidZ{UNSET};
    size_t iForceViscousX{UNSET};
    size_t iForceViscousY{UNSET};
    size_t iForceViscousZ{UNSET};
    if (isBfmSimulation){
        auto getCol = [&columnIndex](const std::string& name) -> size_t {
            auto it = columnIndex.find(name);
            return (it != columnIndex.end()) ? it->second : UNSET;
        };
        iForceInviscidX = getCol("Inviscid Body Force X");
        iForceInviscidY = getCol("Inviscid Body Force Y");
        iForceInviscidZ = getCol("Inviscid Body Force Z");
        iForceViscousX  = getCol("Viscous Body Force X");
        iForceViscousY  = getCol("Viscous Body Force Y");
        iForceViscousZ  = getCol("Viscous Body Force Z");
    }
    
    // Read data
    size_t i = 0, j = 0, k = 0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string value;
        std::vector<FloatType> row;
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }

        // Store coordinates if available
        if (hasCoords && ix < row.size() && iy < row.size() && iz < row.size()) {
            inputCoordinates(i, j, k) = Vector3D(row[ix], row[iy], row[iz]);
        }

        // Store using current indices
        inputDensity(i,j,k)     = row[iDensity];
        inputVelX(i,j,k)        = row[iVelX];
        inputVelY(i,j,k)        = row[iVelY];
        inputVelZ(i,j,k)        = row[iVelZ];
        FloatType totalEnergy = row[iTotEnergy];
        inputTemperature(i,j,k) = _fluid->computeTemperature_rho_u_et(
            inputDensity(i,j,k), 
            {inputVelX(i,j,k), inputVelY(i,j,k), inputVelZ(i,j,k)}, 
            totalEnergy);

        if (isBfmSimulation) {
            inputForceInviscid(i,j,k).x() = (iForceInviscidX != UNSET) ? row[iForceInviscidX] : 0.0;
            inputForceInviscid(i,j,k).y() = (iForceInviscidY != UNSET) ? row[iForceInviscidY] : 0.0;
            inputForceInviscid(i,j,k).z() = (iForceInviscidZ != UNSET) ? row[iForceInviscidZ] : 0.0;

            inputForceViscous(i,j,k).x()  = (iForceViscousX  != UNSET) ? row[iForceViscousX]  : 0.0;
            inputForceViscous(i,j,k).y()  = (iForceViscousY  != UNSET) ? row[iForceViscousY]  : 0.0;
            inputForceViscous(i,j,k).z()  = (iForceViscousZ  != UNSET) ? row[iForceViscousZ]  : 0.0;
        }

        // Update indices: k fastest, then j, then i
        if (++k == NK) {
            k = 0;
            if (++j == NJ) {
                j = 0;
                ++i;
            }
        }
    }

    file.close();
    std::cout << "Data read successfully.\n";
}


void Solver::solve(){
    size_t nIterMax = _config.getMaxIterations();
    Matrix3D<FloatType> timestep(_nPointsI, _nPointsJ, _nPointsK);                          
    std::vector<FloatType> timeIntegrationCoeffs = _config.getTimeIntegrationCoeffs();      
    FlowSolution residuals(_nPointsI, _nPointsJ, _nPointsK);                                
    size_t updateMassFlowsFreq = _config.getSolutionOutputFrequency();                      
    size_t monitorOutputFreq = _config.getSolutionOutputFrequency();
    size_t solutionOutputFreq = _config.getSolutionOutputFrequency();                       
    bool turboOutput = _config.saveTurboOutput();                                           
    bool monitorPointsActive = _config.isMonitorPointsActive();                             
    bool exitLoop = false;                                                                  
    bool steadySimulation = _config.isSimulationSteady();                                   
    if (monitorPointsActive) initializeMonitorPoints();                                     

    // place holder for the solution
    preprocessSolution(_conservativeSolution, false);
    FlowSolution solutionTmp(_conservativeSolution);                              
    std::map<SolutionName, Matrix3D<Vector3D>> solutionGradTmp = _solutionGrad;               
    
    // explict time-stepping
    for (size_t it=1; it<=nIterMax; it++){        
        updateMassFlows(solutionTmp);
        
        if (turboOutput) updateTurboPerformance(solutionTmp);                               
        if (monitorPointsActive) updateMonitorPoints(solutionTmp);                          

        computeTimestepArray(solutionTmp, timestep);                                        
        
        // runge-kutta steps
        preprocessSolution(solutionTmp);
        computeSolutionGradient(solutionTmp, solutionGradTmp);
        updateTurbulenceSolution(solutionTmp, solutionGradTmp, 1.0, timestep);
        for (const auto &integrationCoeff: timeIntegrationCoeffs){
            computeSolutionGradient(solutionTmp, solutionGradTmp);
            computeResiduals(solutionTmp, solutionGradTmp, it, _currentTime, timestep, residuals);
            updateSolution(_conservativeSolution, solutionTmp, residuals, integrationCoeff, timestep);   
            enforcePeriodicityOnSolution(solutionTmp);
        }

        // update the solution and prepare for next iteration
        _conservativeSolution = solutionTmp;
        
        // update the physical time
        _currentTime += timestep.min();
        _time.push_back(_currentTime);
        
        // print information on screen
        printInfoResiduals(residuals, it);
        if (it%updateMassFlowsFreq == 0) {
            printCheckOfMassFlowConservation(it);
            if (turboOutput) printTurboPerformance(it);
        }

        // check the convergence process
        checkConvergence(exitLoop, steadySimulation, it); 
        if (exitLoop && steadySimulation) {
            _output->writeSolution(it);
            writeLogResidualsToCsvFile();
            if (turboOutput) writeTurboPerformanceToCsvFile();
            if (monitorPointsActive) writeMonitorPointsToCsvFile();
            if (_isGreitzerModelingActive) writeGreitzerDynamicsToCsvFile();
            break;
        }

        // write volume output file
        if (it%solutionOutputFreq == 0 || it == nIterMax) {
            _output->writeSolution(it);
        } 

        // write additional text files
        if (it%monitorOutputFreq == 0 || it == nIterMax || _logResiduals.size() >= _historyBufferSize) {
            writeLogResidualsToCsvFile();
            if (turboOutput) writeTurboPerformanceToCsvFile();
            if (monitorPointsActive) writeMonitorPointsToCsvFile();
            if (_isGreitzerModelingActive) writeGreitzerDynamicsToCsvFile();
        } 
        
    }
}

void Solver::printInfoResiduals(FlowSolution &residuals, size_t it) {
    if (it == 1) {printLogResidualsHeader();}
    auto logRes = computeLogResidualNorm(residuals);
    printLogResiduals(logRes, it);
    _logResiduals.push_back(logRes);
    if (!_hasInitialLogResiduals) {
        _initialLogResiduals = logRes;
        _hasInitialLogResiduals = true;
    }
}


void Solver::printCheckOfMassFlowConservation(size_t it) const {
    std::cout << "\nMASS FLOWS CHECK [kg/s]:\n";
    std::cout << "I_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::I_START) << std::endl;
    std::cout << "I_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::I_END) << std::endl;
    std::cout << "J_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::J_START) << std::endl;
    std::cout << "J_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::J_END) << std::endl;
    std::cout << "K_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::K_START) << std::endl;
    std::cout << "K_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::K_END) << std::endl << std::endl;
}

void Solver::printTurboPerformance(size_t it) const {
    std::cout << "\nTURBOMACHINERY PERFORMANCE:\n";
    std::cout << "Mass Flow [kg/s]: " 
              << std::setprecision(6) 
              << _turboPerformance.at(TurboPerformance::MASS_FLOW).back() 
              << std::endl;
    std::cout << "Total Pressure Ratio [-]: " 
              << std::setprecision(6) 
              << _turboPerformance.at(TurboPerformance::TOTAL_PRESSURE_RATIO).back() 
              << std::endl;
    std::cout << "Total Temperature Ratio [-]: " 
              << std::setprecision(6) 
              << _turboPerformance.at(TurboPerformance::TOTAL_TEMPERATURE_RATIO).back() 
              << std::endl;
    std::cout << "Total Efficiency [-]: " 
              << std::setprecision(6) 
              << _turboPerformance.at(TurboPerformance::TOTAL_EFFICIENCY).back() 
              << std::endl; 
    std::cout<< std::endl;
}

void Solver::printLogResiduals(const StateVector &logRes, unsigned long int it) const {
    int col_width = 14;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "|" << std::setw(col_width) << std::setfill(' ') << std::left << it << "|"
              << std::setw(col_width) << std::left << _currentTime*1E6 << "|"
              << std::setw(col_width) << std::right << logRes[0] << "|"
              << std::setw(col_width) << std::right << logRes[1] << "|"
              << std::setw(col_width) << std::right << logRes[2] << "|"
              << std::setw(col_width) << std::right << logRes[3] << "|"
              << std::setw(col_width) << std::right << logRes[4] << "|"
              << std::endl;
}



StateVector Solver::computeLogResidualNorm(const FlowSolution &residuals) const {
    StateVector logResidualNorm{};
    constexpr FloatType minDouble = std::numeric_limits<FloatType>::min();

    for (int i = 0; i < 5; i++) {
        auto residualNorm = residuals.norm(i);
        if (residualNorm >= minDouble) {
            logResidualNorm[i] = std::log10(residualNorm / (_nPointsI * _nPointsJ * _nPointsK));
        } else {
            logResidualNorm[i] = -16.0;
        }
    }
    return logResidualNorm;
}

void Solver::printLogResidualsHeader() const {
    int col_width = 14;
    
    // Print the first separator (top border)
    std::cout << "|" << std::setw(col_width * 7 + 6) << std::setfill('-') << "" << "|" << std::endl;

    // Print the column headers
    std::cout << "|"
                << std::setw(col_width) << std::setfill(' ') << std::left << "Iteration" << "|"
                << std::setw(col_width) << std::left << "Time[μs]" << "|"
                << std::setw(col_width) << std::right << "rms[Rho]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoU]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoV]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoW]" << "|"
                << std::setw(col_width) << std::right << "rms[RhoE]" << "|"
                << std::endl;

    // Print the last separator (bottom border)
    std::cout << "|" << std::setw(col_width * 7 + 6) << std::setfill('-') << "" << "|" << std::endl;
}

void Solver::computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep){
    FloatType cflMax = _config.getCFL();
    bool isViscosityActive = _config.isViscosityActive();

    if (_config.getTimeStepMethod() == TimeStepMethod::FIXED){
        FloatType dtMin = _config.getFixedTimeStep();
        #pragma omp parallel for collapse(2) schedule(static) if(_nPointsI * _nPointsJ >= 64)
        for (size_t i=0; i<_nPointsI; i++) {
            for (size_t j=0; j<_nPointsJ; j++){
                for (size_t k=0; k<_nPointsK; k++){
                    timestep(i,j,k) = dtMin;
                }
            }
        }
        return;
    }

    #pragma omp parallel for collapse(2) schedule(static) if(_nPointsI * _nPointsJ >= 64)
    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++) {
            for (size_t k=0; k<_nPointsK; k++) {

                Vector3D iEdge, jEdge, kEdge;
                _mesh.getElementEdges(i,j,k,iEdge,jEdge,kEdge);

                FloatType dsI = iEdge.magnitude();
                FloatType dsJ = jEdge.magnitude();
                FloatType dsK = kEdge.magnitude();

                Vector3D iDir = iEdge / dsI;
                Vector3D jDir = jEdge / dsJ;
                Vector3D kDir = kEdge / dsK;

                StateVector conservative = solution.at(i,j,k);
                StateVector primitive = getPrimitiveVariablesFromConservative(conservative);

                Vector3D velocity(primitive[1], primitive[2], primitive[3]);

                FloatType rho = primitive[0];
                FloatType temperature = _fluid->computeTemperature_rho_u_et(primitive[0], velocity, primitive[4]);

                FloatType a = _fluid->computeSoundSpeed_rho_u_et(
                                primitive[0],
                                velocity,
                                primitive[4]
                            );
                
                FloatType mu, mut;
                if (isViscosityActive){
                    mu = _fluid->computeMolecularDynamicViscosity(temperature);
                    mut = _turbulenceModel->getEddyViscosity(rho, i, j, k);
                }
                else {
                    mu = 0.0;
                    mut = 0.0;
                }
                
                FloatType nu  = mu  / rho;
                FloatType nut = mut / rho;

                FloatType ui = std::abs(velocity.dot(iDir));
                FloatType uj = std::abs(velocity.dot(jDir));
                FloatType uk = std::abs(velocity.dot(kDir));

                FloatType lambdaConv = (ui + a)/dsI;
                FloatType invH2 = 1.0/(dsI*dsI);

                if (_topology != Topology::ONE_DIMENSIONAL){
                    lambdaConv += (uj + a)/dsJ;
                    invH2 += 1.0/(dsJ*dsJ);
                }

                if (_topology == Topology::THREE_DIMENSIONAL){
                    lambdaConv += (uk + a)/dsK;
                    invH2 += 1.0/(dsK*dsK);
                }

                FloatType lambdaVisc =
                    4.0 * (nu + nut) * invH2;

                timestep(i,j,k) =
                    cflMax /
                    (lambdaConv + lambdaVisc);
            }
        }
    }

    // if time step method is global, the transient is coherent, and dt must be the same for every cell
    if (_config.getTimeStepMethod() == TimeStepMethod::GLOBAL){
        FloatType dtMin = timestep.min();
        #pragma omp parallel for collapse(2) schedule(static)
        for (size_t i=0; i<_nPointsI; i++) {
            for (size_t j=0; j<_nPointsJ; j++){
                for (size_t k=0; k<_nPointsK; k++){
                    timestep(i,j,k) = dtMin;
                }
            }
        }
    }
}


void Solver::updateMassFlows(const FlowSolution&solution){
    std::array<BoundaryIndex, 6> bcIndices {BoundaryIndex::I_START,
                                              BoundaryIndex::I_END,
                                              BoundaryIndex::J_START,
                                              BoundaryIndex::J_END,
                                              BoundaryIndex::K_START,
                                              BoundaryIndex::K_END};
    
    for (auto& bcIndex: bcIndices){
        Matrix2D<Vector3D> surface = _mesh.getMeshBoundary(bcIndex);
        Matrix2D<FloatType> rhoUX = solution._rhoU.getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUV = solution._rhoV.getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUW = solution._rhoW.getBoundarySlice(bcIndex);
        _massFlows[bcIndex] = computeSurfaceIntegral(surface, rhoUX, rhoUV, rhoUW);
    }
    
}


void Solver::updateTurboPerformance(const FlowSolution&solution){
    
    FloatType massFlow = 0.5 * (_massFlows[BoundaryIndex::I_START] + _massFlows[BoundaryIndex::I_END]);
    if (_config.getTopology() == Topology::AXISYMMETRIC){
        massFlow *= 2.0 * M_PI / _mesh.getWedgeAngle();
    }
    else {
        FloatType periodicAngle = _periodicityAngleDeg;
        if (periodicAngle != 0.0) {
            massFlow *= 360.0 / periodicAngle;
        }
    }
    _turboPerformance[TurboPerformance::MASS_FLOW].push_back(massFlow);
    
    std::array<BoundaryIndex, 2> bcIndices {BoundaryIndex::I_START, BoundaryIndex::I_END};
    std::vector<FloatType> totalPressure;
    std::vector<FloatType> totalTemperature;
    for (auto& bcIndex: bcIndices){
        Matrix2D<Vector3D> surface = _mesh.getMeshBoundary(bcIndex);

        size_t nj = surface.sizeI();
        size_t nk = surface.sizeJ();
        Matrix2D<FloatType> rhoUxPt(nj, nk);
        Matrix2D<FloatType> rhoUyPt(nj, nk);
        Matrix2D<FloatType> rhoUzPt(nj, nk);
        Matrix2D<FloatType> rhoUxTt(nj, nk);
        Matrix2D<FloatType> rhoUyTt(nj, nk);
        Matrix2D<FloatType> rhoUzTt(nj, nk);

        for (size_t j=0; j<nj; j++){
            for (size_t k=0; k<nk; k++){
                StateVector primitive;
                if (bcIndex == BoundaryIndex::I_START){
                    primitive = getPrimitiveVariablesFromConservative(solution.at(0,j,k));
                }
                else{
                    primitive = getPrimitiveVariablesFromConservative(solution.at(_nPointsI-1,j,k));
                }
                FloatType rho = primitive[0];
                FloatType ux = primitive[1];
                FloatType uy = primitive[2];
                FloatType uz = primitive[3];
                FloatType et = primitive[4];
                FloatType totalPressure = _fluid->computeTotalPressure_rho_u_et(rho, {ux,uy,uz}, et);
                FloatType totalTemperature = _fluid->computeTotalTemperature_rho_u_et(rho, {ux,uy,uz}, et);

                rhoUxPt(j,k) = rho * ux * totalPressure;
                rhoUyPt(j,k) = rho * uy * totalPressure;
                rhoUzPt(j,k) = rho * uz * totalPressure;
                rhoUxTt(j,k) = rho * ux * totalTemperature;
                rhoUyTt(j,k) = rho * uy * totalTemperature;
                rhoUzTt(j,k) = rho * uz * totalTemperature;
            }
        }
        totalPressure.push_back(computeSurfaceIntegral(surface, rhoUxPt, rhoUyPt, rhoUzPt) / _massFlows[bcIndex]);
        totalTemperature.push_back(computeSurfaceIntegral(surface, rhoUxTt, rhoUyTt, rhoUzTt) / _massFlows[bcIndex]);
    }
    
    FloatType pressureRatio = totalPressure.at(1) / totalPressure.at(0);
    FloatType temperatureRatio = totalTemperature.at(1) / totalTemperature.at(0);
    FloatType efficiency = _fluid->computeTotalEfficiency_PRtt_TRt(pressureRatio, temperatureRatio);

    _turboPerformance[TurboPerformance::TOTAL_PRESSURE_RATIO].push_back(pressureRatio);
    _turboPerformance[TurboPerformance::TOTAL_TEMPERATURE_RATIO].push_back(temperatureRatio);
    _turboPerformance[TurboPerformance::TOTAL_EFFICIENCY].push_back(efficiency);
    
}

void Solver::computeResiduals(
    FlowSolution& solution, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad, 
    const size_t iterationCounter, 
    const FloatType timePhysical, 
    Matrix3D<FloatType> &timestep,
    FlowSolution &residuals) {
        
    residuals.setToZero(); 
    
    // advection fluxes
    computeAdvectionFluxResiduals(FluxDirection::I, solution, iterationCounter, residuals);
    if (_topology!=Topology::ONE_DIMENSIONAL){
        computeAdvectionFluxResiduals(FluxDirection::J, solution, iterationCounter, residuals);
    }
    if (_topology==Topology::THREE_DIMENSIONAL){
        computeAdvectionFluxResiduals(FluxDirection::K, solution, iterationCounter, residuals);
    }

    // viscous fluxes
    if (_config.isViscosityActive()){
        computeViscousFluxResiduals(FluxDirection::I, solution, solutionGrad, iterationCounter, residuals);
        computeViscousFluxResiduals(FluxDirection::J, solution, solutionGrad, iterationCounter, residuals);
        if (_topology==Topology::THREE_DIMENSIONAL){
            computeViscousFluxResiduals(FluxDirection::K, solution, solutionGrad, iterationCounter, residuals);
        }
    }

    // source terms
    computeSourceResiduals(
        solution, 
        solutionGrad, 
        iterationCounter, 
        residuals, 
        _inviscidForce, 
        _viscousForce, 
        _deviationAngle, 
        timePhysical,
        timestep);

    // periodicity enforcement on residuals
    if (_mesh.isPeriodicityActive()){
        enforcePeriodicityOnResiduals(residuals, _periodicityAngleRad);
    }

    // no-slip walls enforcement on residuals
    if (_config.isViscosityActive()){
        enforceNoSlipWallsOnResiduals(residuals);
    }

}


void Solver::enforcePeriodicityOnResiduals(FlowSolution& residuals, FloatType& angleRad) const {
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            StateVector R1 = residuals.at(i, j, 0);
            StateVector R2 = residuals.at(i, j, _nPointsK - 1);

            // Rotate the second in place of a common frame (e.g. frame of the "first" side)
            StateVector R2_frame1 = rotateStateVectorAlongXAxis(R2, -angleRad);

            // Combine residuals to get an average, in the frame of the "first" side
            StateVector R1_avg = (R1 + R2_frame1) * 0.5;

            // Symmetrize: give each half (ensures conservation)
            residuals.set(i, j, 0, R1_avg);
            residuals.set(i, j, _nPointsK - 1, rotateStateVectorAlongXAxis(R1_avg, +angleRad));
        }
    }
}

void Solver::enforceNoSlipWallsOnResiduals(FlowSolution& residuals) const {
    Vector3D zeroWallEffect{0.0, 0.0, 0.0}; 
    for (auto& bc : _boundaries){
        if (bc.type == BoundaryType::NO_SLIP_WALL){
            setMomentumSolutionOnViscousWalls(residuals, bc, zeroWallEffect);
        }
    }
}


template <typename AdvectionScheme>
void Solver::computeAdvectionFluxResidualsImpl(
    const AdvectionScheme& advection,
    FluxDirection direction, 
    const FlowSolution& solution, 
    size_t itCounter, 
    FlowSolution &residuals) const {

    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    const Matrix3D<Vector3D>& midPoints = _mesh.getMidPoints(direction);
    const Matrix3D<std::shared_ptr<BoundaryBase>>& boundaryConditionMap = getBoundaryConditionsMap(direction);
    
    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();

    auto processFace = [&](size_t iFace, size_t jFace, size_t kFace) {
        size_t dirFace = 0;
        size_t stopFace = 0;
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

        if (direction == FluxDirection::K && _isBfmActive) {
            // don't compute the flux in the circumferential direction if the blade is present
            FloatType bladeIsPresent = _mesh.getInputFields(InputField::BLADE_PRESENT, iFace, jFace, 0);
            if (bladeIsPresent > 0.0) {
                return; 
            }
        }
        
        if (dirFace == 0) { // starting boundary fluxes
            StateVector Uinternal = solution.at(iFace, jFace, kFace);
            Vector3D surface = -surfaces(iFace, jFace, kFace); // outward point
            Vector3D midPoint = midPoints(iFace, jFace, kFace);
            StateVector flux = boundaryConditionMap(iFace, jFace, kFace)->computeBoundaryFlux(
                Uinternal, 
                surface, 
                midPoint, 
                {iFace, jFace, kFace}, 
                solution, 
                itCounter);
            residuals.add(iFace, jFace, kFace, flux * surface.magnitude());
        }
        else if (dirFace == stopFace) { // ending boundary fluxes
            StateVector Uinternal = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
            Vector3D surface = surfaces(iFace, jFace, kFace); // outward pointing
            Vector3D midPoint = midPoints(iFace, jFace, kFace);
            size_t ii = iFace-1*stepMask[0];
            size_t jj = jFace-1*stepMask[1];
            size_t kk = kFace-1*stepMask[2];
            switch (direction)
            {                    
            case (FluxDirection::I):
                ii = 1;
                break;
            case (FluxDirection::J):
                jj = 1;
                break;
            case (FluxDirection::K):
                kk = 1;
                break;
            default:
                throw std::runtime_error("Invalid FluxDirection.");
            }
            StateVector flux = boundaryConditionMap(ii, jj, kk)->computeBoundaryFlux(
                Uinternal, 
                surface, 
                midPoint, 
                {iFace, jFace, kFace}, 
                solution, 
                itCounter);
            residuals.add(
                iFace-1*stepMask[0], 
                jFace-1*stepMask[1], 
                kFace-1*stepMask[2], 
                flux * surface.magnitude());
        } 
        else { // flux across internal faces
            StateVector Uleft = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
            StateVector Uright = solution.at(iFace, jFace, kFace);
            StateVector Uleftleft{}, Urightright{};

            // get the extended stencil of values used for muscl reconstruction
            if (direction==FluxDirection::K && _mesh.isPeriodicityActive()){ 
                if (dirFace==1){                                                    // first internal element
                    Uleftleft = solution.at(iFace , jFace , nk-3);                  // last internal element
                    Uleftleft = rotateStateVectorAlongXAxis(Uleftleft, -_mesh.getPeriodicityAngleRad());
                    Urightright = solution.at(iFace , jFace, kFace+1);              // second internal element
                }
                else if (dirFace==stopFace-1){ 
                    Uleftleft = solution.at(iFace, jFace, kFace-2);             // last-1 internal element
                    Urightright = solution.at(iFace, jFace, 1);                 // first internal element
                    Urightright = rotateStateVectorAlongXAxis(Urightright, _mesh.getPeriodicityAngleRad());
                }
                else {
                    Uleftleft = solution.at(iFace-2*stepMask[0], jFace-2*stepMask[1], kFace-2*stepMask[2]);
                    Urightright = solution.at(iFace+1*stepMask[0], jFace+1*stepMask[1], kFace+1*stepMask[2]);
                }
            }
            else { // normal topology -> plain extrapolation for elements close to boundaries
                if (dirFace==1){ 
                    Uleftleft = Uleft - (Uright - Uleft); 
                    Urightright = solution.at(iFace+1*stepMask[0], jFace+1*stepMask[1], kFace+1*stepMask[2]);
                }
                else if (dirFace==stopFace-1){ 
                    Uleftleft = solution.at(iFace-2*stepMask[0], jFace-2*stepMask[1], kFace-2*stepMask[2]);
                    Urightright = Uright + (Uright - Uleft); 
                }
                else {
                    Uleftleft = solution.at(iFace-2*stepMask[0], jFace-2*stepMask[1], kFace-2*stepMask[2]);
                    Urightright = solution.at(iFace+1*stepMask[0], jFace+1*stepMask[1], kFace+1*stepMask[2]);
                }
            }

            Vector3D surface = surfaces(iFace, jFace, kFace); // outward pointing with respect to the left cell
            StateVector flux = advection.computeFlux(Uleftleft, Uleft, Uright, Urightright, surface);
            residuals.add(
                iFace-1*stepMask[0], 
                jFace-1*stepMask[1], 
                kFace-1*stepMask[2], 
                flux * surface.magnitude());
            residuals.subtract(
                iFace, 
                jFace, 
                kFace, 
                flux * surface.magnitude());
        }
    };

    if (direction == FluxDirection::I) {
        #pragma omp parallel for collapse(2) schedule(static) if(nj * nk >= 64)
        for (size_t jFace = 0; jFace < nj; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                for (size_t iFace = 0; iFace < ni; ++iFace) {
                    processFace(iFace, jFace, kFace);
                }
            }
        }
    }
    else if (direction == FluxDirection::J) {
        #pragma omp parallel for collapse(2) schedule(static) if(ni * nk >= 64)
        for (size_t iFace = 0; iFace < ni; ++iFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                for (size_t jFace = 0; jFace < nj; ++jFace) {
                    processFace(iFace, jFace, kFace);
                }
            }
        }
    }
    else { // FluxDirection::K
        #pragma omp parallel for collapse(2) schedule(static) if(ni * nj >= 64)
        for (size_t iFace = 0; iFace < ni; ++iFace) {
            for (size_t jFace = 0; jFace < nj; ++jFace) {
                for (size_t kFace = 0; kFace < nk; ++kFace) {
                    processFace(iFace, jFace, kFace);
                }
            }
        }
    }
}

void Solver::computeAdvectionFluxResiduals(
    FluxDirection direction, 
    const FlowSolution& solution, 
    size_t itCounter, 
    FlowSolution &residuals) const {

    if (const auto* roe = dynamic_cast<const AdvectionRoe*>(_advection.get())) {
        computeAdvectionFluxResidualsImpl(*roe, direction, solution, itCounter, residuals);
    }
    else if (const auto* jst = dynamic_cast<const AdvectionJst*>(_advection.get())) {
        computeAdvectionFluxResidualsImpl(*jst, direction, solution, itCounter, residuals);
    }
    else {
        computeAdvectionFluxResidualsImpl(*_advection, direction, solution, itCounter, residuals);
    }
}


void Solver::computeViscousFluxResiduals(
    FluxDirection direction, 
    const FlowSolution& solution, 
    const std::map<SolutionName, Matrix3D<Vector3D>>& gradients, 
    size_t itCounter, 
    FlowSolution &residuals) const {

    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    
    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();
    
    auto processViscousFace = [&](size_t iFace, size_t jFace, size_t kFace) {
        size_t dirFace = 0;
        size_t stopFace = 0;
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
        
        if (dirFace==0 || dirFace==stopFace){ // skip boundary faces, no viscous flux contribution
            return;
        }

        size_t il = iFace - stepMask[0];
        size_t jl = jFace - stepMask[1];
        size_t kl = kFace - stepMask[2];
        size_t ir = iFace;
        size_t jr = jFace;
        size_t kr = kFace;
        Vector3D surface = surfaces(ir, jr, kr); // pointing from left to right cell
        
        StateVector Uleft = solution.at(il, jl, kl);
        StateVector Uright = solution.at(ir, jr, kr);
        StateVector Uavg = (Uleft + Uright) * 0.5;

        Vector3D uxGrad = (
            gradients.at(SolutionName::VELOCITY_X)(il, jl, kl) + 
            gradients.at(SolutionName::VELOCITY_X)(ir, jr, kr)
            ) * 0.5;
        
        Vector3D uyGrad = (
            gradients.at(SolutionName::VELOCITY_Y)(il, jl, kl) + 
            gradients.at(SolutionName::VELOCITY_Y)(ir, jr, kr)
            ) * 0.5;

        Vector3D uzGrad = (
            gradients.at(SolutionName::VELOCITY_Z)(il, jl, kl) + 
            gradients.at(SolutionName::VELOCITY_Z)(ir, jr, kr)
            ) * 0.5;
        
        Vector3D tempGrad = (
            gradients.at(SolutionName::TEMPERATURE)(il, jl, kl) + 
            gradients.at(SolutionName::TEMPERATURE)(ir, jr, kr)
            ) * 0.5;
        
        FloatType muEddy = (
                _turbulenceModel->getEddyViscosity(solution.at(il, jl, kl)[0], il, jl, kl) +
                _turbulenceModel->getEddyViscosity(solution.at(ir, jr, kr)[0], ir, jr, kr)
                ) * 0.5;
        
        StateVector flux = computeViscousFlux(Uavg, uxGrad, uyGrad, uzGrad, tempGrad, surface, muEddy);
        residuals.add(il, jl, kl, flux * surface.magnitude() * (-1.0));
        residuals.subtract(ir, jr, kr, flux * surface.magnitude() * (-1.0));
    };

    if (direction == FluxDirection::I) {
        #pragma omp parallel for collapse(2) schedule(static) if(nj * nk >= 64)
        for (size_t jFace = 0; jFace < nj; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                for (size_t iFace = 0; iFace < ni; ++iFace) {
                    processViscousFace(iFace, jFace, kFace);
                }
            }
        }
    }
    else if (direction == FluxDirection::J) {
        #pragma omp parallel for collapse(2) schedule(static) if(ni * nk >= 64)
        for (size_t iFace = 0; iFace < ni; ++iFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
                for (size_t jFace = 0; jFace < nj; ++jFace) {
                    processViscousFace(iFace, jFace, kFace);
                }
            }
        }
    }
    else { // FluxDirection::K
        #pragma omp parallel for collapse(2) schedule(static) if(ni * nj >= 64)
        for (size_t iFace = 0; iFace < ni; ++iFace) {
            for (size_t jFace = 0; jFace < nj; ++jFace) {
                for (size_t kFace = 0; kFace < nk; ++kFace) {
                    processViscousFace(iFace, jFace, kFace);
                }
            }
        }
    }
}


StateVector Solver::computeViscousFlux(
    const StateVector& conservative, 
    const Vector3D& velXGrad, 
    const Vector3D& velYGrad, 
    const Vector3D& velZGrad, 
    const Vector3D& tempGrad, 
    const Vector3D& surface,
    const FloatType& eddyViscosity) const{
    
    StateVector primitive = getPrimitiveVariablesFromConservative(conservative);

    Vector3D vel = Vector3D(primitive[1], primitive[2], primitive[3]);
    
    FloatType temperature = _fluid->computeTemperature_rho_u_et(
        primitive[0], 
        vel, 
        primitive[4]);
    
    // flow quantities
    FloatType muL = _fluid->computeMolecularDynamicViscosity(temperature);
    FloatType muEddy = eddyViscosity;

    FloatType kappaL = _fluid->computeThermalConductivity(muL);
    FloatType cp = _config.getFluidHeatCapacity();
    FloatType Prt = _config.getTurbulentPrandtlNumber();
    FloatType kappaEddy = _turbulenceModel->getEddyThermalConductivity(muEddy, cp, Prt);

    // total flow quantities
    FloatType muTotal = muL + muEddy;
    FloatType kappaTotal = kappaL + kappaEddy;
    FloatType secondaryViscosity = -2.0 / 3.0 * muTotal;
    
    ViscousStressTensor tau = computeViscousStressTensor(muTotal, secondaryViscosity, velXGrad, velYGrad, velZGrad);
    Vector3D tauX = Vector3D(tau.xx, tau.xy, tau.xz);
    Vector3D tauY = Vector3D(tau.xy, tau.yy, tau.yz);
    Vector3D tauZ = Vector3D(tau.xz, tau.yz, tau.zz);

    // theta terms (Blazek book pag 17)
    FloatType thetaX = tauX.dot(vel) + kappaTotal * tempGrad.x();
    FloatType thetaY = tauY.dot(vel) + kappaTotal * tempGrad.y();
    FloatType thetaZ = tauZ.dot(vel) + kappaTotal * tempGrad.z();

    StateVector flux({0.0, 0.0, 0.0, 0.0, 0.0});
    Vector3D surfDir = surface / surface.magnitude();
    
    flux[0] = 0.0;
    flux[1] = tauX.dot(surfDir);
    flux[2] = tauY.dot(surfDir);
    flux[3] = tauZ.dot(surfDir);
    flux[4] = thetaX*surfDir.x() + thetaY*surfDir.y() + thetaZ*surfDir.z();

    return flux; 
}



void Solver::updateSolution(
    const FlowSolution &solOld, 
    FlowSolution &solNew, 
    const FlowSolution &residuals, 
    const FloatType &integrationCoeff, 
    const Matrix3D<FloatType> &dt){
    
    // U^{n+1} = U^n - dt / V * R
    const size_t totalCells = _nPointsI * _nPointsJ * _nPointsK;
    const FloatType* d_dt = dt.data();
    const FloatType* d_vol = _mesh.getVolumes().data();

    #pragma omp parallel for schedule(static) if(totalCells >= 64)
    for (size_t idx = 0; idx < totalCells; ++idx) {
        const auto U = solOld.at(idx);
        const auto R = residuals.at(idx);
        const FloatType factor = integrationCoeff * d_dt[idx] / d_vol[idx];
        solNew.set(idx, U - R * factor);
    }
}


void Solver::enforcePeriodicityOnSolution(FlowSolution &solNew){
    if (!_mesh.isPeriodicityActive()) return;

    const auto angle = _periodicityAngleRad;
    const auto inverseAngle = -angle;

    StateVector U1, U2, Uavg;
    for (size_t i = 0; i < _nPointsI; ++i) {
        for (size_t j = 0; j < _nPointsJ; ++j) {

            U1 = solNew.at(i, j, 0);
            U2 = solNew.at(i, j, _nPointsK - 1);

            Uavg = (U1 + rotateStateVectorAlongXAxis(U2, inverseAngle)) * 0.5;

            solNew.set(i, j, 0, Uavg);
            solNew.set(i, j, _nPointsK - 1, rotateStateVectorAlongXAxis(Uavg, angle));
        }
    }
}



void Solver::writeLogResidualsToCsvFile() {

    if (_logResiduals.empty()) return;

    std::string filename = "residuals.csv";
    std::ofstream file;
    if (_isFirstResidualWrite) {
        file.open(filename, std::ios::out);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open log residuals file: " << filename << std::endl;
            return;
        }
        file << "Rho,RhoU,RhoV,RhoW,RhoE\n";
        _isFirstResidualWrite = false;
    } else {
        file.open(filename, std::ios::app);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open log residuals file: " << filename << std::endl;
            return;
        }
    }

    size_t size = _logResiduals.size();
    for (size_t i = 0; i < size; i++) {
        file <<
             _logResiduals[i][0] << "," << 
             _logResiduals[i][1] << "," << 
             _logResiduals[i][2] << "," <<
             _logResiduals[i][3] << "," <<
             _logResiduals[i][4] << "\n";
    }

    file.close();
    _logResiduals.clear();
    std::cout << std::endl;
    std::cout << "Log residuals appended to " << filename << std::endl;
    std::cout << std::endl;
}

void Solver::writeTurboPerformanceToCsvFile() {

    if (_turboPerformance.empty() || _turboPerformance.at(TurboPerformance::MASS_FLOW).empty()) return;

    std::string filename = "turbo.csv";
    std::ofstream file;
    if (_isFirstTurboWrite) {
        file.open(filename, std::ios::out);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open turbo performance file: " << filename << std::endl;
            return;
        }
        file << "Time[μs],Massflow[kg/s],PRtt,TRtt,ETAtt\n";
        _isFirstTurboWrite = false;
    } else {
        file.open(filename, std::ios::app);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open turbo performance file: " << filename << std::endl;
            return;
        }
    }

    size_t size = _turboPerformance.at(TurboPerformance::MASS_FLOW).size();
    for (size_t i = 0; i < size; i++) {
        FloatType timeVal = (i < _time.size()) ? _time.at(i) : _currentTime;
        file << timeVal*1E6 << ","
             << _turboPerformance.at(TurboPerformance::MASS_FLOW)[i] << "," 
             << _turboPerformance.at(TurboPerformance::TOTAL_PRESSURE_RATIO)[i] << "," 
             << _turboPerformance.at(TurboPerformance::TOTAL_TEMPERATURE_RATIO)[i] << "," 
             << _turboPerformance.at(TurboPerformance::TOTAL_EFFICIENCY)[i] << "\n"; 
    }

    file.close();
    for (auto& [key, vec] : _turboPerformance) {
        vec.clear();
    }
    _time.clear();
    std::cout << std::endl;
    std::cout << "Turbo performance appended to " << filename << std::endl;
    std::cout << std::endl;
}


void Solver::writeGreitzerDynamicsToCsvFile() {

    if (!_greitzerModel || _greitzerModel->getSize() == 0) return;

    std::string filename = "greitzer_dynamics.csv";
    std::ofstream file;
    if (_isFirstGreitzerWrite) {
        file.open(filename, std::ios::out);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open greitzer dynamics file: " << filename << std::endl;
            return;
        }
        file << "Time[s],PlenumPressure[Pa],PlenumInletMassflow[kg/s],PlenumOutletMassflow[kg/s]\n";
        _isFirstGreitzerWrite = false;
    } else {
        file.open(filename, std::ios::app);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open greitzer dynamics file: " << filename << std::endl;
            return;
        }
    }

    size_t size = _greitzerModel->getSize();
    for (size_t i = 0; i < size; i++) {
        file << _greitzerModel->getTime(i) << ","
             << _greitzerModel->getPlenumPressure(i) << ","
             << _greitzerModel->getPlenumInletMassflow(i) << ","
             << _greitzerModel->getPlenumOutletMassflow(i) << "\n";
    }

    file.close();
    _greitzerModel->clearBuffer();
    std::cout << std::endl;
    std::cout << "Greitzer dynamics appended to " << filename << std::endl;
    std::cout << std::endl;
}



void Solver::writeMonitorPointsToCsvFile() {

    if (_monitorPoints.empty()) return;

    std::string folder = "Monitor_Points";
    std::filesystem::create_directories(folder); // Ensure the folder exists

    for (size_t iPoint = 0; iPoint < _monitorPoints.size(); iPoint++) {
        if (_monitorPoints[iPoint].empty() || _monitorPoints[iPoint].at(MonitorOutputField::TIME).empty()) continue;

        std::string filename = folder + "/Monitor_Point_" + std::to_string(iPoint) + ".csv";
        std::ofstream file;
        if (_isFirstMonitorPointsWrite) {
            file.open(filename, std::ios::out);
            if (!file.is_open()) {
                std::cerr << "Error: Could not open monitor point file: " << filename << std::endl;
                return;
            }
            file << "Time[s],Pressure[Pa],Velocity_X[m/s],Velocity_Y[m/s],Velocity_Z[m/s]\n";
        } else {
            file.open(filename, std::ios::app);
            if (!file.is_open()) {
                std::cerr << "Error: Could not open monitor point file: " << filename << std::endl;
                return;
            }
        }

        file << std::scientific << std::setprecision(6);
        size_t size = _monitorPoints[iPoint].at(MonitorOutputField::TIME).size();
        for (size_t i = 0; i < size; i++) {
            file    << _monitorPoints[iPoint].at(MonitorOutputField::TIME)[i] << "," 
                    << _monitorPoints[iPoint].at(MonitorOutputField::PRESSURE)[i] << "," 
                    << _monitorPoints[iPoint].at(MonitorOutputField::VELOCITY_X)[i] << "," 
                    << _monitorPoints[iPoint].at(MonitorOutputField::VELOCITY_Y)[i] << "," 
                    << _monitorPoints[iPoint].at(MonitorOutputField::VELOCITY_Z)[i] << "\n"; 
        }
        file.close();

        for (auto& [field, vec] : _monitorPoints[iPoint]) {
            vec.clear();
        }
    }

    _isFirstMonitorPointsWrite = false;
    std::cout << std::endl;
    std::cout << "Written monitor points to " << folder  << std::endl;
    std::cout << std::endl;

}


void Solver::updateRadialProfiles(FlowSolution &solution){
    StateVector conservative, primitive;
    Vector3D velocityCart, velocityCyl;
    if (_isGreitzerModelingActive){
        FloatType mflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
        _hubStaticPressure = _greitzerModel->computePlenumPressure(mflow);
    }
    
    for (auto& radialProfile : _radialEquilibriumProfiles){
        std::vector<FloatType> densityProfile(radialProfile.pressure.size());
        std::vector<FloatType> velTangProfile(radialProfile.pressure.size());
        FloatType theta;
        size_t i = radialProfile.boundary.i_min-1;
        size_t k = 0;
        for (size_t j = radialProfile.boundary.j_min; j < radialProfile.boundary.j_max; j++) {
            conservative = solution.at(i, j, 0);
            primitive = getPrimitiveVariablesFromConservative(conservative);    
            velocityCart(0) = primitive[1];
            velocityCart(1) = primitive[2];
            velocityCart(2) = primitive[3];
            theta = _mesh.getTheta(i, j, 0);
            velocityCyl = computeCylindricalComponentsFromCartesian(velocityCart, theta);
            
            densityProfile[k] = primitive[0];
            velTangProfile[k] = std::abs(velocityCyl.z());
            k++;
        }
        if (_isGreitzerModelingActive){
            // _hubStaticPressure already updated from Greitzer model
        }
        else if (radialProfile.boundary.type == BoundaryType::THROTTLE){
            FloatType kt = radialProfile.boundary.values[0];
            FloatType mflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
            FloatType Pt_in=0.0;
            for (auto& boundary: _boundaries){
                if (boundary.type == BoundaryType::INLET){
                    Pt_in = boundary.values[0];
                    break;
                }
            }
            _hubStaticPressure = Pt_in + kt * mflow*mflow;
        }
        else if (radialProfile.boundary.type == BoundaryType::RADIAL_EQUILIBRIUM){  
            _hubStaticPressure = radialProfile.pressure[0];
        }

        integrateRadialEquilibrium(
            densityProfile, 
            velTangProfile, 
            radialProfile.radius, 
            radialProfile.pressure,
            _hubStaticPressure);
        
    }
}


void Solver::computeSourceResiduals(
    FlowSolution& solution, 
    const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad, 
    const size_t itCounter, 
    FlowSolution &residuals, Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    Matrix3D<FloatType> &deviationAngle, 
    FloatType timePhysical,
    Matrix3D<FloatType> &timestep) {
    
    if (_topology == Topology::THREE_DIMENSIONAL && _isBfmActive == false) {
        return;
    }
    
    StateVector primitive, conservative, bfmSource, sourceGeometrical;
    FloatType omega, radius, theta, volume, pressure;
    Vector3D densityGrad, velXGrad, velYGrad, velZGrad, totEnergyGrad;
    Vector3D blockageGradient, velocityCart, velocityCyl, sourceCyl, sourceCart;
    StateVector gongSource;
    bool geometricSourceFlag{false};
    bool gongSourceFlag{false};
    const auto& surfacesK = _mesh.getSurfacesK();
    const FloatType wedgeAngle = _mesh.getWedgeAngle();

    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {

                conservative = solution.at(i, j, k);
                primitive = getPrimitiveVariablesFromConservative(conservative);
                volume = _mesh.getVolume(i, j, k);
                pressure = _fluid->computePressure_rho_u_et(
                    primitive[0], 
                    {primitive[1], primitive[2], primitive[3]}, 
                    primitive[4]);
                radius = _mesh.getRadius(i, j, k);
                theta = _mesh.getTheta(i, j, k);

                understandWhatSourcesAreNeeded(i, j, k, geometricSourceFlag, gongSourceFlag);
        
                if (geometricSourceFlag) {
                    velocityCart = {primitive[1], primitive[2], primitive[3]};
                    velocityCyl = computeCylindricalComponentsFromCartesian(velocityCart, theta);

                    if (_topology == Topology::AXISYMMETRIC) {
                        // By Pappus's centroid theorem, the volume of a wedge cell is: V = r_centroid * A_K * deltaTheta.
                        // Therefore, the volume-integrated geometric source terms:
                        //   int (S / r) dV = S * A_K * deltaTheta
                        // which cancels out 'r' analytically and avoids coordinate singularity near the axis (r -> 0).
                        FloatType areaK = surfacesK(i, j, 0).magnitude();
                        FloatType geomFactor = areaK * wedgeAngle;

                        sourceCyl.x() = 0.0;
                        sourceCyl.y() = (+ primitive[0] * velocityCyl.z() * velocityCyl.z() + pressure) * geomFactor;
                        sourceCyl.z() = - primitive[0] * velocityCyl.y() * velocityCyl.z() * geomFactor;

                        sourceCart = computeCartesianComponentsFromCylindrical(sourceCyl, theta);
                        sourceGeometrical[0] = 0.0;
                        sourceGeometrical[1] = sourceCart.x();
                        sourceGeometrical[2] = sourceCart.y();
                        sourceGeometrical[3] = sourceCart.z();
                        sourceGeometrical[4] = 0.0;

                        residuals.subtract(i, j, k, sourceGeometrical);
                    } else {
                        sourceCyl.x() = 0.0; 
                        sourceCyl.y() = (+ primitive[0] * velocityCyl.z() * velocityCyl.z() + pressure) / radius; 
                        sourceCyl.z() = - primitive[0] * velocityCyl.y() * velocityCyl.z() / radius; 

                        sourceCart = computeCartesianComponentsFromCylindrical(sourceCyl, theta);
                        sourceGeometrical[0] = 0.0;
                        sourceGeometrical[1] = sourceCart.x();
                        sourceGeometrical[2] = sourceCart.y();
                        sourceGeometrical[3] = sourceCart.z();
                        sourceGeometrical[4] = 0.0;
                        
                        residuals.subtract(i, j, k, sourceGeometrical * volume);
                    }
                }

                if (gongSourceFlag){
                    omega = _mesh.getInputFields(InputField::RPM, i, j, k) * 2 * M_PI / 60;
                    FloatType scalingFactor = _config.getRotationalSpeedScalingFactor(timePhysical);
                    omega *= scalingFactor;
                    gongSource = computeGongSource(radius, theta, omega, i, j, k, volume, solution);
                    residuals.subtract(i, j, k, gongSource);
                }

                if (_isBfmActive){
                    blockageGradient = _mesh.getInputFieldsGradient(InputField::BLOCKAGE, i, j, k); 
                    if (blockageGradient.magnitude() > 1E-10){                     
                        bfmSource = _bfmSource->computeTotalSource(
                            i, j, k, 
                            primitive, 
                            inviscidForce, 
                            viscousForce, 
                            deviationAngle, 
                            timePhysical, 
                            solution,
                            timestep,
                            timePhysical);
                        residuals.subtract(i, j, k, bfmSource);
                    }
                }
            }
        }
    }
}

void Solver::understandWhatSourcesAreNeeded(
    size_t i, size_t j, size_t k, 
    bool &geometricSourceFlag, 
    bool &gongSourceFlag) const {

    if (_topology == Topology::THREE_DIMENSIONAL && _isBfmActive) {
        FloatType bladePresent = _mesh.getInputFields(InputField::BLADE_PRESENT, i, j, k);
        if (bladePresent > 0.0) {
            geometricSourceFlag = true;
            gongSourceFlag = true;
        }
        else {
            geometricSourceFlag = false;
            gongSourceFlag = false;
        }
    }
    else if (_topology == Topology::AXISYMMETRIC){
        geometricSourceFlag = true;
        gongSourceFlag = false;
    }
    else {
        geometricSourceFlag = false;
        gongSourceFlag = false;
    }

}

void Solver::initializeMonitorPoints(){

    size_t seedI = _config.getMonitorPointsCoordsI();
    size_t seedJ = _config.getMonitorPointsCoordsJ();
    _monitorPointsIdxI.push_back(seedI);
    _monitorPointsIdxJ.push_back(seedJ);
    _monitorPointsIdxK.push_back(0);

    if (_topology == Topology::THREE_DIMENSIONAL){
        size_t circumferentialPoints = _config.getCircumferentialNumberMonitorPoints();
        size_t deltaK = _nPointsK / circumferentialPoints;
        int leftOver = _nPointsK%circumferentialPoints;

        if (leftOver > 0){
            std::cerr << "Error: Monitor points not evenly distributed" 
                      << "Please choose a divisor of the total circumferential points" << std::endl;
            return;
        }

        FloatType deltaAngle = _mesh.getPeriodicityAngleDeg();
        if (deltaAngle==0.0) {
            deltaAngle = 360.0;
        }
        std::cout << "The angle between monitor points is: " 
                  << deltaAngle / (circumferentialPoints) 
                  << " degrees" << std::endl;

        for (size_t k = 1; k < circumferentialPoints; k++){
            _monitorPointsIdxI.push_back(seedI);
            _monitorPointsIdxJ.push_back(seedJ);
            _monitorPointsIdxK.push_back(k * deltaK);
        }
    }

    _numberMonitorPoints = _monitorPointsIdxI.size();
    _monitorPoints.resize(_numberMonitorPoints);

}

void Solver::updateMonitorPoints(const FlowSolution &solution){
    for (unsigned int i = 0; i < _numberMonitorPoints; i++){
        size_t idxI = _monitorPointsIdxI[i];
        size_t idxJ = _monitorPointsIdxJ[i];
        size_t idxK = _monitorPointsIdxK[i];
        StateVector conservative = solution.at(idxI, idxJ, idxK);
        StateVector primitive = getPrimitiveVariablesFromConservative(conservative);
        FloatType pressure = _fluid->computePressure_rho_u_et(
            primitive[0], 
            {primitive[1], primitive[2], primitive[3]}, 
            primitive[4]);

        _monitorPoints[i][MonitorOutputField::PRESSURE].push_back(pressure);
        _monitorPoints[i][MonitorOutputField::VELOCITY_X].push_back(primitive[1]);
        _monitorPoints[i][MonitorOutputField::VELOCITY_Y].push_back(primitive[2]);
        _monitorPoints[i][MonitorOutputField::VELOCITY_Z].push_back(primitive[3]);
        _monitorPoints[i][MonitorOutputField::TIME].push_back(_currentTime);
    }
}

void Solver::checkConvergence(bool &exitLoop, bool &isSteady, size_t it) const {
    if (!isSteady || _logResiduals.empty()) return;

    StateVector current = _logResiduals.back();

    if (current[0] < _initialLogResiduals[0] - _residualsDropConvergence &&
        current[1] < _initialLogResiduals[1] - _residualsDropConvergence &&
        current[2] < _initialLogResiduals[2] - _residualsDropConvergence &&
        current[3] < _initialLogResiduals[3] - _residualsDropConvergence &&
        current[4] < _initialLogResiduals[4] - _residualsDropConvergence) {
        std::cout << "\nConvergence reached at iteration " << it << std::endl;
        std::cout << std::endl;
        exitLoop = true;
    } 
}


void Solver::computeGradientOfField(const Matrix3D<FloatType> &var, Matrix3D<Vector3D> &grad) const {
    computeGradientGreenGauss(  
        _mesh.getSurfacesI(), 
        _mesh.getSurfacesJ(), 
        _mesh.getSurfacesK(), 
        _mesh.getMidPointsI(), 
        _mesh.getMidPointsJ(), 
        _mesh.getMidPointsK(), 
        _mesh.getVertices(), 
        _mesh.getVolumes(), 
        var, 
        grad);
}


void Solver::computeSolutionGradient(FlowSolution &sol, std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad){
    if (!_config.isViscosityActive()){
        return;
    }

    Matrix3D<FloatType> rho = sol.getDensity();
    Matrix3D<FloatType> ux = sol.getVelocityX();
    Matrix3D<FloatType> uy = sol.getVelocityY();
    Matrix3D<FloatType> uz = sol.getVelocityZ();
    Matrix3D<FloatType> et = sol.getTotalEnergy();
    Matrix3D<FloatType> temperature = _fluid->computeTemperature_conservative(rho, ux, uy, uz, et);

    computeGradientsGreenGauss4(
        _mesh.getSurfacesI(), 
        _mesh.getSurfacesJ(), 
        _mesh.getSurfacesK(), 
        _mesh.getMidPointsI(), 
        _mesh.getMidPointsJ(), 
        _mesh.getMidPointsK(), 
        _mesh.getVertices(), 
        _mesh.getVolumes(), 
        ux, uy, uz, temperature,
        solutionGrad[SolutionName::VELOCITY_X],
        solutionGrad[SolutionName::VELOCITY_Y],
        solutionGrad[SolutionName::VELOCITY_Z],
        solutionGrad[SolutionName::TEMPERATURE]);
}

void Solver::updateTurbulenceSolution(
    FlowSolution &sol, std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad,
    const FloatType &integrationCoeff, 
    const Matrix3D<FloatType> &dt){
    
    if (!_config.isTurbulenceActive()){
        return;
    }
    
    _turbulenceModel->solve(sol, solutionGrad, dt*integrationCoeff);  
}


StateVector Solver::computeGongSource(
    const FloatType& radius,
    const FloatType& theta,
    const FloatType& omega,
    const size_t i,
    const size_t j,
    const size_t k,
    const FloatType& volume,
    const FlowSolution& solution) const{

    if (std::abs(omega)<1E-3){
        // stator blades -> no source term 
        return StateVector({0,0,0,0,0});
    }

    if (_periodicityAngleDeg != 0){
        std::cerr << "Error: Periodicity angle different from full annulus is not supported for Gong BF formulation\n";
    }

    StateVector U0, U1, U2;
    if (omega > 0) {
        if (k==0){ // first point
            U0 = solution.at(i, j, k);
            U1 = solution.at(i, j, _nPointsK-2);
            U2 = solution.at(i, j, _nPointsK-3);
        }
        else if (k==1){ // second point
            U0 = solution.at(i, j, k);
            U1 = solution.at(i, j, 0);
            U2 = solution.at(i, j, _nPointsK-2);
        }
        else { // internals
            U0 = solution.at(i, j, k);
            U1 = solution.at(i, j, k-1);
            U2 = solution.at(i, j, k-2);
        }
    } else {
        if (k==_nPointsK-1){ // last point
            U0 = solution.at(i, j, 0);
            U1 = solution.at(i, j, 1);
            U2 = solution.at(i, j, 2);
        }
        else if (k==_nPointsK-2){ // second to last
            U0 = solution.at(i, j, k);
            U1 = solution.at(i, j, 0);
            U2 = solution.at(i, j, 1);
        }
        else {
            U0 = solution.at(i, j, k);
            U1 = solution.at(i, j, k+1);
            U2 = solution.at(i, j, k+2);
        }
    }

    FloatType dTheta = (2 * M_PI / static_cast<FloatType>(_nPointsK-1));
    
    StateVector dU_dTheta({0,0,0,0,0});
    if (omega > 0.0) {
        // backward stencil → backward derivative
        dU_dTheta = ( U0*3.0 - U1*4.0 + U2 ) / (2*dTheta);
    } else {
        // forward stencil → forward derivative
        dU_dTheta = ( U0*(-3.0) + U1*4.0 - U2 ) / (2*dTheta);
    }

    StateVector source({0,0,0,0,0});
    for (size_t i = 0; i < 5; i++) {
        source[i] = -dU_dTheta[i] * omega;
    }

    return source*volume;
}

void Solver::setMomentumSolutionOnViscousWalls(
    FlowSolution &sol, 
    const Boundary &boundary,
    const Vector3D& wallVelocity) const{
    
    size_t iStart, iLast, jStart, jLast, kStart, kLast;
    iStart = boundary.i_min;
    iLast = boundary.i_max;
    jStart = boundary.j_min;
    jLast = boundary.j_max;
    kStart = boundary.k_min;
    kLast = boundary.k_max;

    auto applyWallCondition = [&](size_t i, size_t j, size_t k) {
    FloatType density = sol.at(i, j, k)[0];
    for (int eq = 1; eq <= 3; ++eq)
        sol.set(i, j, k, eq, density * wallVelocity(eq - 1));
    };

    if (iStart == iLast) {
        size_t i = (iStart == 0) ? 0 : iLast - 1;
        for (size_t j = jStart; j < jLast; ++j)
            for (size_t k = kStart; k < kLast; ++k)
                applyWallCondition(i, j, k);
    }
    else if (jStart == jLast) {
        size_t j = (jStart == 0) ? 0 : jLast - 1;
        for (size_t i = iStart; i < iLast; ++i)
            for (size_t k = kStart; k < kLast; ++k)
                applyWallCondition(i, j, k);
    }
    else if (kStart == kLast) {
        size_t k = (kStart == 0) ? 0 : kLast - 1;
        for (size_t i = iStart; i < iLast; ++i)
            for (size_t j = jStart; j < jLast; ++j)
                applyWallCondition(i, j, k);
    }


}


void Solver::preprocessSolution(FlowSolution &sol, bool updateRadialProf) {
    
    if (updateRadialProf) {
        updateRadialProfiles(sol);
    };

    // if there are no-slip walls impose zero velocity
    Vector3D wallVel(0.0, 0.0, 0.0);
    if (_config.isViscosityActive()){
        for (auto& bc : _boundaries){
            if (bc.type == BoundaryType::NO_SLIP_WALL){
                wallVel = {bc.values[0], bc.values[1], bc.values[2]};
                setMomentumSolutionOnViscousWalls(sol, bc, wallVel);
            }
        }
    }
    
}