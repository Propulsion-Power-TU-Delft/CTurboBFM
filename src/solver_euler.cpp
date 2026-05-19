#include "solver_euler.hpp"
#include "math_utils.hpp"
#include "types.hpp"
#include <iostream> // Optional, for logging/debugging
#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

SolverEuler::SolverEuler(Config& config, Mesh& mesh)
    : SolverBase(config, mesh)  
{
    _isGreitzerModelingActive = _config.enableGreitzerModeling();
    if (_isGreitzerModelingActive) {
        _greitzerModel = std::make_unique<GreitzerModel>(_config, *_fluid);
    }
    
    initializeSolutionArrays();
    
    _output = std::make_unique<OutputCSV>(
        _config, 
        _mesh, 
        _conservativeSolution, 
        *_fluid, _inviscidForce, 
        _viscousForce, 
        _deviationAngle);

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

    _isTurbulenceActive = _config.isTurbulenceActive();
    TurbulenceModel turbulenceModel = _config.getTurbulenceModel();
    if (turbulenceModel == TurbulenceModel::SPALART_ALLMARAS) {
        _turbulenceModel = std::make_unique<TurbulenceModelSA>();
    }
    else if (turbulenceModel == TurbulenceModel::NONE) {
        _turbulenceModel = std::make_unique<TurbulenceModelBase>();
    }
    else {
        throw std::runtime_error("Unsupported turbulence model selected.");
    }

}


void SolverEuler::initializeSolutionArrays(){
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

    _radialProfilePressure.resize(_nPointsJ);
    _radialProfileRadialCoords.resize(_nPointsJ);
    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    for (size_t j = 0; j < nj; j++) {
        _radialProfileRadialCoords[j] = std::sqrt(_mesh.getVertex(ni-1,j,0).y() * _mesh.getVertex(ni-1,j,0).y() + 
                                            _mesh.getVertex(ni-1,j,0).z() * _mesh.getVertex(ni-1,j,0).z());
    }


    if (_isGreitzerModelingActive){
        updateMassFlows(_conservativeSolution);
        updateTurboPerformance(_conservativeSolution);
        FloatType massflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
        StateVector conservativeHubOutlet = _conservativeSolution.at(_nPointsI-1, 0, 0);
        StateVector primitiveHubOutlet = getPrimitiveVariablesFromConservative(conservativeHubOutlet);
        FloatType outletPressureHub = _fluid->computePressure_primitive(primitiveHubOutlet);
        _greitzerModel->initializeState(outletPressureHub, massflow, massflow);
    }

}


void SolverEuler::initializeSolutionFromScratch(){
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

void SolverEuler::initializeSolutionFromRestart(){
    std::string restartFileName = _config.getRestartFilepath();
    
    size_t NI=0, NJ=0, NK=0;    
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
        inputDensity, 
        inputVelX, 
        inputVelY, 
        inputVelZ, 
        inputTemperature, 
        inputForceViscous, 
        inputForceInviscid);

    if (NI != _nPointsI || NJ != _nPointsJ || NK != _nPointsK) {
        if (NI == _nPointsI && NJ == _nPointsJ && _config.getRestartType()=="axisymmetric") {
            axisymmetricRestart(inputDensity, inputVelX, inputVelY, inputVelZ, inputTemperature);
        }
        else if (NI == _nPointsI && NJ == _nPointsJ && _config.getRestartType()!="axisymmetric") {
            std::cerr << "Restart file dimensions (I,J) coincides with solver dimensions, but K does not." << 
                        "If you want to restart in axisymmetric mode, specify RESTART_TYPE=axisymmetric\n";
            exit(1);
        }
        else {
            std::cerr << "Restart file dimensions do not match solver dimensions.\n";
            exit(1);
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


void SolverEuler::standardRestart(
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


void SolverEuler::axisymmetricRestart(
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


void SolverEuler::readRestartFile(
    const std::string &restartFileName, 
    size_t &NI, 
    size_t &NJ, 
    size_t &NK,
    Matrix3D<FloatType> &inputDensity, 
    Matrix3D<FloatType> &inputVelX, 
    Matrix3D<FloatType> &inputVelY, 
    Matrix3D<FloatType> &inputVelZ, 
    Matrix3D<FloatType> &inputTemperature, 
    Matrix3D<Vector3D> &inputForceViscous, 
    Matrix3D<Vector3D> &inputForceInviscid) {

    std::ifstream file(restartFileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open file.\n";
    }

    std::string line;

    // Read NI, NJ, NK
    std::getline(file, line); NI = std::stoi(line.substr(line.find('=') + 1));
    std::getline(file, line); NJ = std::stoi(line.substr(line.find('=') + 1));
    std::getline(file, line); NK = std::stoi(line.substr(line.find('=') + 1));

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
        columnIndex[column] = idx++;
    }

    // Get indexes of the fields of interest (primary fields should always be available in the restart)
    size_t iDensity   = columnIndex.at("Density");
    size_t iVelX      = columnIndex.at("Velocity X");
    size_t iVelY      = columnIndex.at("Velocity Y");
    size_t iVelZ      = columnIndex.at("Velocity Z");
    size_t iTotEnergy = columnIndex.at("Total Energy");

    // these could also not be there, not a problem
    size_t iForceInviscidX{1000};
    size_t iForceInviscidY{1000};
    size_t iForceInviscidZ{1000};
    size_t iForceViscousX{1000};
    size_t iForceViscousY{1000};
    size_t iForceViscousZ{1000};
    if (isBfmSimulation){
        iForceInviscidX = columnIndex["Inviscid Body Force X"];
        iForceInviscidY = columnIndex["Inviscid Body Force Y"];
        iForceInviscidZ = columnIndex["Inviscid Body Force Z"];
        iForceViscousX = columnIndex["Viscous Body Force X"];
        iForceViscousY = columnIndex["Viscous Body Force Y"];
        iForceViscousZ = columnIndex["Viscous Body Force Z"];
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
            inputForceInviscid(i,j,k).x() = (iForceInviscidX != 0) ? row[iForceInviscidX] : 0.0;
            inputForceInviscid(i,j,k).y() = (iForceInviscidY != 0) ? row[iForceInviscidY] : 0.0;
            inputForceInviscid(i,j,k).z() = (iForceInviscidZ != 0) ? row[iForceInviscidZ] : 0.0;

            inputForceViscous(i,j,k).x()  = (iForceViscousX != 0) ? row[iForceViscousX] : 0.0;
            inputForceViscous(i,j,k).y()  = (iForceViscousY != 0) ? row[iForceViscousY] : 0.0;
            inputForceViscous(i,j,k).z()  = (iForceViscousZ != 0) ? row[iForceViscousZ] : 0.0;
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


void SolverEuler::solve(){
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
        for (const auto &integrationCoeff: timeIntegrationCoeffs){
            computeSolutionGradient(solutionTmp, solutionGradTmp);
            computeResiduals(solutionTmp, solutionGradTmp, it, _time.back(), timestep, residuals);
            updateSolution(_conservativeSolution, solutionTmp, residuals, integrationCoeff, timestep);   
            enforcePeriodicityOnSolution(solutionTmp);
        }

        // update the solution and prepare for next iteration
        _conservativeSolution = solutionTmp;
        
        // update the physical time
        _time.push_back(_time.back() + timestep.min());
        
        // print information on screen
        printInfoResiduals(residuals, it);
        if (it%updateMassFlowsFreq == 0) {
            printCheckOfMassFlowConservation(it);
            if (turboOutput) printTurboPerformance(it);
        }

        // check the convergence process
        checkConvergence(exitLoop, steadySimulation); 
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
            _output->writeSolution(it, false);
        } 

        // write additional text files
        if (it%monitorOutputFreq == 0) {
            writeLogResidualsToCsvFile();
            if (turboOutput) writeTurboPerformanceToCsvFile();
            if (monitorPointsActive) writeMonitorPointsToCsvFile();
            if (_isGreitzerModelingActive) writeGreitzerDynamicsToCsvFile();
        } 
        
    }
}

void SolverEuler::printInfoResiduals(FlowSolution &residuals, size_t it) {
    if (it == 1) {printLogResidualsHeader();}
    auto logRes = computeLogResidualNorm(residuals);
    printLogResiduals(logRes, it);
    _logResiduals.push_back(logRes);
}


void SolverEuler::printCheckOfMassFlowConservation(size_t it) const {
    std::cout << "\nMASS FLOWS CHECK [kg/s]:\n";
    std::cout << "I_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::I_START) << std::endl;
    std::cout << "I_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::I_END) << std::endl;
    std::cout << "J_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::J_START) << std::endl;
    std::cout << "J_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::J_END) << std::endl;
    std::cout << "K_START: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::K_START) << std::endl;
    std::cout << "K_END: " << std::setprecision(6) << _massFlows.at(BoundaryIndex::K_END) << std::endl << std::endl;
}

void SolverEuler::printTurboPerformance(size_t it) const {
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

void SolverEuler::printLogResiduals(const StateVector &logRes, unsigned long int it) const {
    int col_width = 14;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "|" << std::setw(col_width) << std::setfill(' ') << std::left << it << "|"
              << std::setw(col_width) << std::left << _time.back()*1E6 << "|"
              << std::setw(col_width) << std::right << logRes[0] << "|"
              << std::setw(col_width) << std::right << logRes[1] << "|"
              << std::setw(col_width) << std::right << logRes[2] << "|"
              << std::setw(col_width) << std::right << logRes[3] << "|"
              << std::setw(col_width) << std::right << logRes[4] << "|"
              << std::endl;
}



StateVector SolverEuler::computeLogResidualNorm(const FlowSolution &residuals) const {
    StateVector logResidualNorm{};
    constexpr double minDouble = std::numeric_limits<FloatType>::min();

    for (int i = 0; i < 5; i++) {
        auto residualNorm = residuals.norm(i);
        if (residualNorm >= minDouble) {
            logResidualNorm[i] = std::log10(residualNorm / (_nPointsI * _nPointsJ * _nPointsK));
        } else {
            logResidualNorm[i] = 0.0;
        }
    }
    return logResidualNorm;
}

void SolverEuler::printLogResidualsHeader() const {
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

void SolverEuler::computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep){
    FloatType cflMax = _config.getCFL();
    Vector3D iEdge, jEdge, kEdge;
    Vector3D iDir, jDir, kDir;
    Vector3D velocity;
    StateVector primitive, conservative;
    FloatType soundSpeed;
    std::array<FloatType, 3> dtEdge;
    FloatType dtMin;

    if (_config.getTimeStepMethod() == TimeStepMethod::FIXED){
        dtMin = _config.getFixedTimeStep();
        for (size_t i=0; i<_nPointsI; i++) {
            for (size_t j=0; j<_nPointsJ; j++){
                for (size_t k=0; k<_nPointsK; k++){
                    timestep(i,j,k) = dtMin;
                }
            }
        }
        return;
    }

    for (size_t i=0; i<_nPointsI; i++) {
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                _mesh.getElementEdges(i, j, k, iEdge, jEdge, kEdge);
                iDir = iEdge / iEdge.magnitude();
                jDir = jEdge / jEdge.magnitude();
                kDir = kEdge / kEdge.magnitude();
                conservative = solution.at(i,j,k);
                primitive = getPrimitiveVariablesFromConservative(conservative);
                velocity(0) = primitive[1];
                velocity(1) = primitive[2];
                velocity(2) = primitive[3];
                soundSpeed = _fluid->computeSoundSpeed_rho_u_et(primitive[0], velocity, primitive[4]);

                dtEdge[0] = iEdge.magnitude() / (std::abs(velocity.dot(iDir))+soundSpeed);
                dtEdge[1] = jEdge.magnitude() / (std::abs(velocity.dot(jDir))+soundSpeed);
                dtEdge[2] = kEdge.magnitude() / (std::abs(velocity.dot(kDir))+soundSpeed);

                dtMin = *std::min_element(dtEdge.begin(), dtEdge.end());
                timestep(i,j,k) = dtMin*cflMax;
            }
        }
    }

    // if time step method is global, the transient is coherent, and dt must be the same for every cell
    if (_config.getTimeStepMethod() == TimeStepMethod::GLOBAL){
        dtMin = timestep.min();
        for (size_t i=0; i<_nPointsI; i++) {
            for (size_t j=0; j<_nPointsJ; j++){
                for (size_t k=0; k<_nPointsK; k++){
                    timestep(i,j,k) = dtMin;
                }
            }
        }
    }
}


void SolverEuler::updateMassFlows(const FlowSolution&solution){
    std::array<BoundaryIndex, 6> bcIndices {BoundaryIndex::I_START,
                                              BoundaryIndex::I_END,
                                              BoundaryIndex::J_START,
                                              BoundaryIndex::J_END,
                                              BoundaryIndex::K_START,
                                              BoundaryIndex::K_END};
    
    for (auto& bcIndex: bcIndices){
        Matrix2D<Vector3D> surface = _mesh.getMeshBoundary(bcIndex);
        Matrix2D<FloatType> rhoUX = (_conservativeSolution._rhoU).getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUV = (_conservativeSolution._rhoV).getBoundarySlice(bcIndex);
        Matrix2D<FloatType> rhoUW = (_conservativeSolution._rhoW).getBoundarySlice(bcIndex);
        _massFlows[bcIndex] = computeSurfaceIntegral(surface, rhoUX, rhoUV, rhoUW);
    }
    
}


void SolverEuler::updateTurboPerformance(const FlowSolution&solution){
    
    FloatType massFlow = 0.5 * (_massFlows[BoundaryIndex::I_START] + _massFlows[BoundaryIndex::I_END]);
    if (_config.getTopology() == Topology::AXISYMMETRIC){
        massFlow *= 2.0 * M_PI / _mesh.getWedgeAngle();
    }
    else {
        FloatType periodicAngle = _config.getPeriodicityAngleDeg();
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

void SolverEuler::computeResiduals(
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
        FloatType angle = _mesh.getPeriodicityAngleRad();
        enforcePeriodicityOnResiduals(residuals, angle);
    }

    // no-slip walls enforcement on residuals
    if (_config.isViscosityActive()){
        enforceNoSlipWallsOnResiduals(residuals);
    }

}


void SolverEuler::enforcePeriodicityOnResiduals(FlowSolution& residuals, FloatType& angleRad) const {
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

void SolverEuler::enforceNoSlipWallsOnResiduals(FlowSolution& residuals) const {
    Vector3D zeroWallEffect{0.0, 0.0, 0.0}; 
    for (auto& bc : _boundaryTypes){
        if (bc.second == BoundaryType::NO_SLIP_WALL){
            setMomentumSolutionOnViscousWalls(residuals, bc.first, zeroWallEffect);
        }
    }
}


void SolverEuler::computeAdvectionFluxResiduals(
    FluxDirection direction, 
    const FlowSolution& solution, 
    size_t itCounter, 
    FlowSolution &residuals) const {

    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    const Matrix3D<Vector3D>& midPoints = _mesh.getMidPoints(direction);
    
    BoundaryIndex boundaryStart, boundaryEnd;
    if (direction==FluxDirection::I){
        boundaryStart = BoundaryIndex::I_START;
        boundaryEnd = BoundaryIndex::I_END;
    }
    else if (direction==FluxDirection::J){
        boundaryStart = BoundaryIndex::J_START;
        boundaryEnd = BoundaryIndex::J_END;
    }
    else {
        boundaryStart = BoundaryIndex::K_START;
        boundaryEnd = BoundaryIndex::K_END;
    }
    
    StateVector Uinternal{}, Uleft{}, Uright{}, Uleftleft {}, Urightright {}, flux {};
    Vector3D surface {}, midPoint {};

    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();
    for (size_t iFace = 0; iFace < ni; ++iFace) {
        for (size_t jFace = 0; jFace < nj; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
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
                        continue; 
                    }
                }
                
                
                if (dirFace == 0) { // starting boundary fluxes
                    Uinternal = solution.at(iFace, jFace, kFace);
                    surface = -surfaces(iFace, jFace, kFace);
                    midPoint = midPoints(iFace, jFace, kFace);
                    flux = _boundaryConditions.at(boundaryStart)->computeBoundaryFlux(
                        Uinternal, 
                        surface, 
                        midPoint, 
                        {iFace, jFace, kFace}, 
                        solution, 
                        itCounter);
                    residuals.add(iFace, jFace, kFace, flux * surface.magnitude());
                } 
                else if (dirFace == stopFace) { // ending boundary fluxes
                    Uinternal = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    surface = surfaces(iFace, jFace, kFace);
                    midPoint = midPoints(iFace, jFace, kFace);
                    flux = _boundaryConditions.at(boundaryEnd)->computeBoundaryFlux(
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
                    Uleft = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    Uright = solution.at(iFace, jFace, kFace);

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
        
                    surface = surfaces(iFace, jFace, kFace);
                    flux = _advection->computeFlux(Uleftleft, Uleft, Uright, Urightright, surface);
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

            
            }
        }
    }
    
}


void SolverEuler::computeViscousFluxResiduals(
    FluxDirection direction, 
    const FlowSolution& solution, 
    const std::map<SolutionName, Matrix3D<Vector3D>>& gradients, 
    size_t itCounter, 
    FlowSolution &residuals) const {

    const auto stepMask = getStepMask(direction);
    const Matrix3D<Vector3D>& surfaces = _mesh.getSurfaces(direction);
    
    StateVector Uleft{}, Uright{}, Uavg{}, flux {};
    Vector3D surface {};
    Vector3D uxGrad, uyGrad, uzGrad, tempGrad;

    size_t ni = surfaces.sizeI(); 
    size_t nj = surfaces.sizeJ(); 
    size_t nk = surfaces.sizeK();
    
    // The viscous fluxes have no contribution here to the boundary elements, since their effect is already
    // treated with the boundary conditions in the advection part. Extra things (such as no-slip walls) are
    // also treated elsewhere --> so only internal fluxes are computed here.
    for (size_t iFace = 1; iFace < ni-1; ++iFace) {
        for (size_t jFace = 1; jFace < nj-1; ++jFace) {
            for (size_t kFace = 0; kFace < nk; ++kFace) {
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
                
                if (dirFace == 0) {
                    continue;
                } 
                else if (dirFace == stopFace) { 
                    continue;
                } 
                else {
                    Uleft = solution.at(iFace-1*stepMask[0], jFace-1*stepMask[1], kFace-1*stepMask[2]);
                    Uright = solution.at(iFace, jFace, kFace);
                    Uavg = (Uleft + Uright) * 0.5;

                    uxGrad = (gradients.at(SolutionName::VELOCITY_X)(
                        iFace-1*stepMask[0], 
                        jFace-1*stepMask[1], 
                        kFace-1*stepMask[2]) + 
                        gradients.at(SolutionName::VELOCITY_X)(iFace, jFace, kFace)) * 0.5;
                    
                    uyGrad = (gradients.at(SolutionName::VELOCITY_Y)(
                        iFace-1*stepMask[0], 
                        jFace-1*stepMask[1], 
                        kFace-1*stepMask[2]) + 
                        gradients.at(SolutionName::VELOCITY_Y)(iFace, jFace, kFace)) * 0.5;

                    uzGrad = (gradients.at(SolutionName::VELOCITY_Z)(
                        iFace-1*stepMask[0], 
                        jFace-1*stepMask[1], 
                        kFace-1*stepMask[2]) + 
                        gradients.at(SolutionName::VELOCITY_Z)(iFace, jFace, kFace)) * 0.5;
                    
                    tempGrad = (gradients.at(SolutionName::TEMPERATURE)(
                        iFace-1*stepMask[0], 
                        jFace-1*stepMask[1], 
                        kFace-1*stepMask[2]) + 
                        gradients.at(SolutionName::TEMPERATURE)(iFace, jFace, kFace)) * 0.5;
                    
                    surface = surfaces(iFace, jFace, kFace);
                    flux = computeViscousFlux(Uavg, uxGrad, uyGrad, uzGrad, tempGrad, surface);
                    
                    residuals.add(
                        iFace-1*stepMask[0], 
                        jFace-1*stepMask[1], 
                        kFace-1*stepMask[2], 
                        flux * surface.magnitude());
                    residuals.subtract(iFace, 
                        jFace, 
                        kFace, 
                        flux * surface.magnitude());
                
                }
            
            }
        }
    }
    
}


StateVector SolverEuler::computeViscousFlux(
    const StateVector& conservative, 
    const Vector3D& velXGrad, 
    const Vector3D& velYGrad, 
    const Vector3D& velZGrad, 
    const Vector3D& tempGrad, 
    const Vector3D& surface) const{
    
    StateVector primitive = getPrimitiveVariablesFromConservative(conservative);
    Vector3D vel = Vector3D(primitive[1], primitive[2], primitive[3]);
    
    // fluid quantities
    FloatType nu = _config.getFluidKinematicViscosity();
    FloatType mu = nu * primitive[0]; 
    FloatType lmbda = -2.0 / 3.0 * mu;         
    FloatType cp = _config.getFluidHeatCapacity(); 
    FloatType Pr = _config.getFluidPrandtlNumber();             
    FloatType kappa = cp * mu / Pr;   
    
    // viscous stresses
    FloatType tauxx, tauyy, tauzz, tauxy, tauxz, tauyz;

    FloatType divVel = velXGrad.x() + velYGrad.y() + velZGrad.z();
    
    // viscous stresses
    tauxx = lmbda * divVel + 2.0 * mu * velXGrad.x();
    tauyy = lmbda * divVel + 2.0 * mu * velYGrad.y();
    tauzz = lmbda * divVel + 2.0 * mu * velZGrad.z();

    tauxy = mu * (velXGrad.y() + velYGrad.x());
    tauxz = mu * (velXGrad.z() + velZGrad.x());
    tauyz = mu * (velYGrad.z() + velZGrad.y());

    Vector3D tauX = Vector3D(tauxx, tauxy, tauxz);
    Vector3D tauY = Vector3D(tauxy, tauyy, tauyz);
    Vector3D tauZ = Vector3D(tauxz, tauyz, tauzz);

    // theta terms (Blazek pag 17)
    FloatType thetaX = tauX.dot(vel) + kappa * tempGrad.x();
    FloatType thetaY = tauY.dot(vel) + kappa * tempGrad.y();
    FloatType thetaZ = tauZ.dot(vel) + kappa * tempGrad.z();

    StateVector flux({0.0, 0.0, 0.0, 0.0, 0.0});
    Vector3D surfDir = surface / surface.magnitude();
    
    flux[0] = 0.0;
    flux[1] = tauX.dot(surfDir);
    flux[2] = tauY.dot(surfDir);
    flux[3] = tauZ.dot(surfDir);
    flux[4] = thetaX*surfDir.x() + thetaY*surfDir.y() + thetaZ*surfDir.z();

    return flux*(-1.0); // minus due to viscous flux positive considered on the left hand side of equations

}



void SolverEuler::updateSolution(
    const FlowSolution &solOld, 
    FlowSolution &solNew, 
    const FlowSolution &residuals, 
    const FloatType &integrationCoeff, 
    const Matrix3D<FloatType> &dt){
    
    // U^{n+1} = U^n - dt / V * R
    for (size_t i = 0; i < _nPointsI; ++i) {
        for (size_t j = 0; j < _nPointsJ; ++j) {
            for (size_t k = 0; k < _nPointsK; ++k) {

                const auto U  = solOld.at(i,j,k);
                const auto R  = residuals.at(i,j,k);
                const auto dti = dt(i,j,k);
                const auto vol = _mesh.getVolume(i,j,k);

                solNew.set(i,j,k, U - R * (integrationCoeff * dti / vol));
            }
        }
    }
}


void SolverEuler::enforcePeriodicityOnSolution(FlowSolution &solNew){
    if (!_mesh.isPeriodicityActive()) return;

    const auto angle = _mesh.getPeriodicityAngleRad();
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



void SolverEuler::writeLogResidualsToCsvFile() const {

    std::string filename = "residuals.csv";
    std::ofstream file(filename); // open in truncate (default) mode

    if (!file.is_open()) {
        std::cerr << "Error: Could not open log residuals file: " << filename << std::endl;
        return;
    }

    // Write header
    file << "Rho,RhoU,RhoV,RhoW,RhoE\n";

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
    std::cout << std::endl;
    std::cout << "Log residuals written to " << filename << std::endl;
    std::cout << std::endl;
}

void SolverEuler::writeTurboPerformanceToCsvFile() const {

    std::string filename = "turbo.csv";
    std::ofstream file(filename); // open in truncate (default) mode

    if (!file.is_open()) {
        std::cerr << "Error: Could not open turbo performance file: " << filename << std::endl;
        return;
    }

    // Write header
    file << "Time[μs],Massflow[kg/s],PRtt,TRtt,ETAtt\n";

    size_t size = _turboPerformance.at(TurboPerformance::MASS_FLOW).size();
    for (size_t i = 0; i < size; i++) {
        file << _time.at(i)*1E6 << ",";
        file    << _turboPerformance.at(TurboPerformance::MASS_FLOW)[i] << "," 
                << _turboPerformance.at(TurboPerformance::TOTAL_PRESSURE_RATIO)[i] << "," 
                << _turboPerformance.at(TurboPerformance::TOTAL_TEMPERATURE_RATIO)[i] << "," 
                << _turboPerformance.at(TurboPerformance::TOTAL_EFFICIENCY)[i] << std::endl; 
    }

    file.close();
    std::cout << std::endl;
    std::cout << "Turbo performance written to " << filename << std::endl;
    std::cout << std::endl;
}


void SolverEuler::writeGreitzerDynamicsToCsvFile() const {

    std::string filename = "greitzer_dynamics.csv";
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open turbo performance file: " << filename << std::endl;
        return;
    }

    // Write header
    file << "Time[s],PlenumPressure[Pa],PlenumInletMassflow[kg/s],PlenumOutletMassflow[kg/s]\n";
    size_t size = _greitzerModel->getSize();
    for (size_t i = 0; i < size; i++) {
        file << _greitzerModel->getTime(i) << ",";
        file << _greitzerModel->getPlenumPressure(i) << ",";
        file << _greitzerModel->getPlenumInletMassflow(i) << ",";
        file << _greitzerModel->getPlenumOutletMassflow(i) << std::endl;
    }

    file.close();
    std::cout << std::endl;
    std::cout << "Greitzer dynamics written to " << filename << std::endl;
    std::cout << std::endl;
}



void SolverEuler::writeMonitorPointsToCsvFile() const {

    std::string folder = "Monitor_Points";
    std::filesystem::create_directories(folder); // Ensure the folder exists

    for (size_t iPoint = 0; iPoint < _monitorPoints.size(); iPoint++) {
        std::string filename = folder + "/Monitor_Point_" + std::to_string(iPoint) + ".csv";
        std::ofstream file(filename); 

        if (!file.is_open()) {
            std::cerr << "Error: Could not open turbo performance file: " << filename << std::endl;
            return;
        }

        // Write header
        file << "Time[s],Pressure[Pa],Velocity_X[m/s],Velocity_Y[m/s],Velocity_Z[m/s]\n";
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
    }

    std::cout << std::endl;
    std::cout << "Written monitor points to " << folder  << std::endl;
    std::cout << std::endl;

}


void SolverEuler::updateRadialProfiles(FlowSolution &solution){
    StateVector conservative, primitive;
    Vector3D velocityCart, velocityCyl;
    std::vector<FloatType> densityProfile(_nPointsJ);
    std::vector<FloatType> velTangProfile(_nPointsJ);
    FloatType theta;

    for (size_t j = 0; j < _nPointsJ; j++) {
        conservative = solution.at(_nPointsI-1, j, 0);
        primitive = getPrimitiveVariablesFromConservative(conservative);    
        velocityCart(0) = primitive[1];
        velocityCart(1) = primitive[2];
        velocityCart(2) = primitive[3];
        densityProfile[j] = primitive[0];
        theta = _mesh.getTheta(_nPointsI-1,j,0);
        velocityCyl = computeCylindricalComponentsFromCartesian(velocityCart, theta);
        velTangProfile[j] = std::abs(velocityCyl.z());
    }

    if (_boundaryTypes[BoundaryIndex::I_END] == BoundaryType::THROTTLE && _isGreitzerModelingActive==false){
        FloatType mflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
        FloatType totPressureInlet = _config.getInletBCValues().at(0);
        FloatType throttleCoeff = _config.getOutletBCValues()[0];
        _hubStaticPressure = totPressureInlet + throttleCoeff * mflow*mflow;
    }
    else if (_boundaryTypes[BoundaryIndex::I_END] == BoundaryType::THROTTLE && _isGreitzerModelingActive==true){
        FloatType mflow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
        _hubStaticPressure = _greitzerModel->computePlenumPressure(mflow);
    }
    else {
        // nothing needed in other cases
        }

    integrateRadialEquilibrium(
        densityProfile, 
        velTangProfile, 
        _radialProfileRadialCoords, 
        _hubStaticPressure, 
        _radialProfilePressure);

}


void SolverEuler::computeSourceResiduals(
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

                    sourceCyl.x() = 0.0; 
                    sourceCyl.y() = (+ primitive[0] * velocityCyl.z() * velocityCyl.z() + pressure) / radius; 
                    sourceCyl.z() = - primitive[0] * velocityCyl.y() * velocityCyl.z() / radius; 

                    sourceCart = computeCartesianComponentsFromCylindrical(sourceCyl, theta);
                    sourceGeometrical[0] = 0.0;
                    sourceGeometrical[1] = sourceCart.x();
                    sourceGeometrical[2] = sourceCart.y();
                    sourceGeometrical[3] = sourceCart.z();
                    sourceGeometrical[4] = 0.0;
                    
                    residuals.subtract(i, j, k, sourceGeometrical*volume);
                }

                if (gongSourceFlag){
                    omega = _mesh.getInputFields(InputField::RPM, i, j, k) * 2 * M_PI / 60;
                    FloatType scalingFactor = _config.getRotationalSpeedScalingFactor(timePhysical);
                    omega *= scalingFactor;
                    gongSource = computeGongSource(radius, theta, omega, i, j, k, volume);
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

void SolverEuler::understandWhatSourcesAreNeeded(
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

void SolverEuler::initializeMonitorPoints(){

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

void SolverEuler::updateMonitorPoints(const FlowSolution &solution){
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
        _monitorPoints[i][MonitorOutputField::TIME].push_back(_time.back());
    }
}

void SolverEuler::checkConvergence(bool &exitLoop, bool &isSteady) const {
    if (!isSteady) return;

    StateVector current = _logResiduals.back();
    StateVector initial = _logResiduals.front();

    if (current[0] < initial[0] - _residualsDropConvergence &&
        current[1] < initial[1] - _residualsDropConvergence &&
        current[2] < initial[2] - _residualsDropConvergence &&
        current[4] < initial[4] - _residualsDropConvergence) {
        std::cout << "\nConvergence reached at iteration " << _logResiduals.size() << std::endl;
        std::cout << std::endl;
        exitLoop = true;
    } 
}


void SolverEuler::computeGradientOfField(const Matrix3D<FloatType> &var, Matrix3D<Vector3D> &grad) const {
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


void SolverEuler::computeSolutionGradient(FlowSolution &sol, std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad){
    if (!_config.isViscosityActive()){
        return;
    }

    Matrix3D<FloatType> rho = sol.getDensity();
    Matrix3D<FloatType> ux = sol.getVelocityX();
    Matrix3D<FloatType> uy = sol.getVelocityY();
    Matrix3D<FloatType> uz = sol.getVelocityZ();
    Matrix3D<FloatType> et = sol.getTotalEnergy();

    computeGradientOfField(ux, solutionGrad[SolutionName::VELOCITY_X]);
    computeGradientOfField(uy, solutionGrad[SolutionName::VELOCITY_Y]);
    computeGradientOfField(uz, solutionGrad[SolutionName::VELOCITY_Z]);
    Matrix3D<FloatType> temperature = _fluid->computeTemperature_conservative(rho, ux, uy, uz, et);
    computeGradientOfField(temperature, solutionGrad[SolutionName::TEMPERATURE]);

    
}


StateVector SolverEuler::computeGongSource(
    const FloatType& radius, 
    const FloatType& theta, 
    const FloatType& omega, 
    const size_t i, 
    const size_t j, 
    const size_t k, 
    const FloatType& volume) const{

    if (std::abs(omega)<1E-3){
        // stator blades -> no source term 
        return StateVector({0,0,0,0,0});
    }

    if (_config.getPeriodicityAngleDeg() != 0){
        std::cerr << "Error: Periodicity angle different from full annulus is not supported for Gong BF formulation\n";
    }

    StateVector U0, U1, U2;
    if (omega > 0) {
        if (k==0){ // first point
            U0 = _conservativeSolution.at(i, j, k);
            U1 = _conservativeSolution.at(i, j, _nPointsK-2);
            U2 = _conservativeSolution.at(i, j, _nPointsK-3);
        }
        else if (k==1){ // second point
            U0 = _conservativeSolution.at(i, j, k);
            U1 = _conservativeSolution.at(i, j, 0);
            U2 = _conservativeSolution.at(i, j, _nPointsK-2);
        }
        else { // internals
            U0 = _conservativeSolution.at(i, j, k);
            U1 = _conservativeSolution.at(i, j, k-1);
            U2 = _conservativeSolution.at(i, j, k-2);
        }
    } else {
        if (k==_nPointsK-1){ // last point
            U0 = _conservativeSolution.at(i, j, 0);
            U1 = _conservativeSolution.at(i, j, 1);
            U2 = _conservativeSolution.at(i, j, 2);
        }
        else if (k==_nPointsK-2){ // second to last
            U0 = _conservativeSolution.at(i, j, k);
            U1 = _conservativeSolution.at(i, j, 0);
            U2 = _conservativeSolution.at(i, j, 1);
        }
        else {
            U0 = _conservativeSolution.at(i, j, k);
            U1 = _conservativeSolution.at(i, j, k+1);
            U2 = _conservativeSolution.at(i, j, k+2);
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



void SolverEuler::writeSolution(size_t iterationCounter, bool alsoGradients){
    _output->writeSolution(iterationCounter, alsoGradients);
}




void SolverEuler::setMomentumSolutionOnViscousWalls(
    FlowSolution &sol, 
    const BoundaryIndex &boundaryIndex, 
    const Vector3D& wallVelocity) const{
    
    size_t iStart, iLast, jStart, jLast, kStart, kLast;
    getBoundarySliceIndices(boundaryIndex, iStart, iLast, jStart, jLast, kStart, kLast);

    FloatType density{0.0};
    for (size_t i = iStart; i < iLast; i++) {
        for (size_t j = jStart; j < jLast; j++) {
            for (size_t k = kStart; k < kLast; k++) {
                density = sol.at(i, j, k)[0];
                for (int eq = 1; eq <= 3; ++eq){
                    sol.set(i, j, k, eq, density * wallVelocity(eq-1));
                }
            }
        }
    }
}


void SolverEuler::preprocessSolution(FlowSolution &sol, bool updateRadialProf) {
    
    if (updateRadialProf) {
        updateRadialProfiles(sol);
    };

    // if there are no-slip walls impose zero velocity
    Vector3D wallVel(0.0, 0.0, 0.0);
    if (_config.isViscosityActive()){
        for (auto& bc : _boundaryTypes){
            if (bc.second == BoundaryType::NO_SLIP_WALL){
                wallVel = _config.getNoSlipWallVelocity(bc.first);
                setMomentumSolutionOnViscousWalls(sol, bc.first, wallVel);
            }
        }
    }
    
}