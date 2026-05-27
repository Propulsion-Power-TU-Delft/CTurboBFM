#include "output.hpp"

namespace {
    constexpr const char* DENSITY                   = "Density";
    constexpr const char* VELOCITY_X                = "Velocity X";
    constexpr const char* VELOCITY_Y                = "Velocity Y";
    constexpr const char* VELOCITY_Z                = "Velocity Z";
    constexpr const char* TOTAL_ENERGY              = "Total Energy";
    constexpr const char* PRESSURE                  = "Pressure";
    constexpr const char* TEMPERATURE               = "Temperature";
    constexpr const char* MACH                      = "Mach";
    constexpr const char* TOTAL_PRESSURE            = "Total Pressure";
    constexpr const char* TOTAL_TEMPERATURE         = "Total Temperature";
    constexpr const char* ENTROPY                   = "Entropy";
    constexpr const char* RELATIVE_VELOCITY_X       = "Relative Velocity X";
    constexpr const char* RELATIVE_VELOCITY_Y       = "Relative Velocity Y";
    constexpr const char* RELATIVE_VELOCITY_Z       = "Relative Velocity Z";
    constexpr const char* GRID_VELOCITY_X           = "Grid Velocity X";
    constexpr const char* GRID_VELOCITY_Y           = "Grid Velocity Y";
    constexpr const char* GRID_VELOCITY_Z           = "Grid Velocity Z";
    constexpr const char* RELATIVE_MACH             = "Relative Mach";
    constexpr const char* VISCOUS_BODY_FORCE_X      = "Viscous Body Force X";
    constexpr const char* VISCOUS_BODY_FORCE_Y      = "Viscous Body Force Y";
    constexpr const char* VISCOUS_BODY_FORCE_Z      = "Viscous Body Force Z";
    constexpr const char* INVISCID_BODY_FORCE_X     = "Inviscid Body Force X";
    constexpr const char* INVISCID_BODY_FORCE_Y     = "Inviscid Body Force Y";
    constexpr const char* INVISCID_BODY_FORCE_Z     = "Inviscid Body Force Z";
    constexpr const char* BLOCKAGE                  = "Blockage";
    constexpr const char* DEVIATION_ANGLE           = "Deviation Angle";
}

Output::Output(
    const Config &config, 
    const Mesh &mesh, 
    const FlowSolution &solution, 
    const FluidBase &fluid, 
    const Matrix3D<Vector3D> &inviscidForce, 
    const Matrix3D<Vector3D> &viscousForce,
    const Matrix3D<FloatType> &deviationAngle)
    : _config(config), 
    _mesh(mesh), 
    _solution(solution), 
    _fluid(fluid), 
    _inviscidForce(inviscidForce), 
    _viscousForce(viscousForce), 
    _deviationAngle(deviationAngle) {    

    std::filesystem::create_directory(_outputVolumeDirectory);
    _outputFieldsType = _config.getOutputFieldsType();
    _isUnsteadyOutput = _config.saveUnsteadySolution();
    _ni = _mesh.getNumberPointsI();
    _nj = _mesh.getNumberPointsJ();
    _nk = _mesh.getNumberPointsK();
    _isBfmActive = _config.isBFMActive();
    allocateOutputArrays();
    }

void Output::allocateOutputArrays() {
    
    // Primary solution variables
    _outputFields.emplace(DENSITY,                      Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_X,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_Y,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_Z,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(TOTAL_ENERGY,                 Matrix3D<FloatType>(_ni, _nj, _nk));
    
    // Secondary solution variables
    if (_outputFieldsType == OutputFieldsType::SECONDARY || _outputFieldsType == OutputFieldsType::TURBO_BFM){
        _outputFields.emplace(PRESSURE,                 Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(TEMPERATURE,              Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(MACH,                     Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(TOTAL_PRESSURE,           Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(TOTAL_TEMPERATURE,        Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(ENTROPY,                  Matrix3D<FloatType>(_ni, _nj, _nk));
    }

    // Turbo BFM variables
    const bool isBFMActive = _config.isBFMActive();
    if (isBFMActive && _outputFieldsType == OutputFieldsType::TURBO_BFM){
        _outputFields.emplace(BLOCKAGE,                 Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(RELATIVE_MACH,            Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(GRID_VELOCITY_X,          Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(GRID_VELOCITY_Y,          Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(GRID_VELOCITY_Z,          Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(RELATIVE_VELOCITY_X,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(RELATIVE_VELOCITY_Y,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(RELATIVE_VELOCITY_Z,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(VISCOUS_BODY_FORCE_X,     Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(VISCOUS_BODY_FORCE_Y,     Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(VISCOUS_BODY_FORCE_Z,     Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(INVISCID_BODY_FORCE_X,    Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(INVISCID_BODY_FORCE_Y,    Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(INVISCID_BODY_FORCE_Z,    Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(DEVIATION_ANGLE,          Matrix3D<FloatType>(_ni, _nj, _nk));
    }
}
    

void Output::updateOutputArrays() {
    
    // Primary
    _outputFields[DENSITY] = _solution.getDensity();
    _outputFields[VELOCITY_X] = _solution.getVelocityX();
    _outputFields[VELOCITY_Y] = _solution.getVelocityY();
    _outputFields[VELOCITY_Z] = _solution.getVelocityZ();
    _outputFields[TOTAL_ENERGY] = _solution.getTotalEnergy();

    // Others
    if (_outputFieldsType != OutputFieldsType::PRIMARY){
        Vector3D vel, gridVelCyl, gridVelCart, relVel;
        FloatType rho, et, omega, radius, theta;
        
        for (size_t i = 0; i < _ni; ++i) {
            for (size_t j = 0; j < _nj; ++j) {
                for (size_t k = 0; k < _nk; ++k) {
                
                    rho = _outputFields[DENSITY](i, j, k);
                    vel.x() = _outputFields[VELOCITY_X](i, j, k);
                    vel.y() = _outputFields[VELOCITY_Y](i, j, k);
                    vel.z() = _outputFields[VELOCITY_Z](i, j, k);
                    et  = _outputFields[TOTAL_ENERGY](i, j, k);

                    _outputFields[PRESSURE](i, j, k) = _fluid.computePressure_rho_u_et(rho, vel, et);
                    _outputFields[TEMPERATURE](i, j, k) = _fluid.computeTemperature_rho_u_et(rho, vel, et);
                    _outputFields[MACH](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, vel, et);
                    _outputFields[TOTAL_PRESSURE](i, j, k) = _fluid.computeTotalPressure_rho_u_et(rho, vel, et);
                    _outputFields[TOTAL_TEMPERATURE](i, j, k) = _fluid.computeTotalTemperature_rho_u_et(rho, vel, et);
                    _outputFields[ENTROPY](i, j, k) = _fluid.computeEntropy_rho_u_et(rho, vel, et);
                    
                    // only for BFM
                    if (_isBfmActive && _outputFieldsType == OutputFieldsType::TURBO_BFM){

                        _outputFields[BLOCKAGE](i, j, k) = _mesh.getInputFields(InputField::BLOCKAGE)(i, j, k);
                        omega = _mesh.getInputFields(InputField::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                        FloatType scalingFactor = _config.getRotationalSpeedScalingFactor();
                        omega *= scalingFactor;
                        radius = std::sqrt(_mesh.getVertex(i, j, k).z() * _mesh.getVertex(i, j, k).z() +
                                                    _mesh.getVertex(i, j, k).y() * _mesh.getVertex(i, j, k).y());
                        theta = std::atan2(_mesh.getVertex(i, j, k).z(), _mesh.getVertex(i, j, k).y());
                        gridVelCyl = {0.0, 0.0, omega * radius};
                        gridVelCart = computeCartesianComponentsFromCylindrical(gridVelCyl, theta);

                        _outputFields[GRID_VELOCITY_X](i, j, k) = gridVelCart.x();
                        _outputFields[GRID_VELOCITY_Y](i, j, k) = gridVelCart.y();
                        _outputFields[GRID_VELOCITY_Z](i, j, k) = gridVelCart.z();
                        
                        relVel = vel - gridVelCart;
                        _outputFields[RELATIVE_VELOCITY_X](i, j, k) = relVel.x();
                        _outputFields[RELATIVE_VELOCITY_Y](i, j, k) = relVel.y();
                        _outputFields[RELATIVE_VELOCITY_Z](i, j, k) = relVel.z();

                        _outputFields[RELATIVE_MACH](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, relVel, et);

                        _outputFields[VISCOUS_BODY_FORCE_X](i, j, k) = _viscousForce(i, j, k).x();
                        _outputFields[VISCOUS_BODY_FORCE_Y](i, j, k) = _viscousForce(i, j, k).y();
                        _outputFields[VISCOUS_BODY_FORCE_Z](i, j, k) = _viscousForce(i, j, k).z();

                        _outputFields[INVISCID_BODY_FORCE_X](i, j, k) = _inviscidForce(i, j, k).x();
                        _outputFields[INVISCID_BODY_FORCE_Y](i, j, k) = _inviscidForce(i, j, k).y();
                        _outputFields[INVISCID_BODY_FORCE_Z](i, j, k) = _inviscidForce(i, j, k).z();

                        _outputFields[DEVIATION_ANGLE](i, j, k) = _deviationAngle(i, j, k);
                    }
                }
            }
        }
    }
}




std::string Output::getOutputFilename(size_t iterationCounter) {
    std::string filename;
    if (_isUnsteadyOutput){
        std::ostringstream oss;
        oss << _config.getSolutionName() << "_" << std::setw(6) << std::setfill('0') << iterationCounter;
        filename = oss.str();
    }
    else {
        filename = _config.getSolutionName();
    }

    return filename;
}

void Output::writeSolution(size_t iterationCounter){
    updateOutputArrays();

    std::string filename = getOutputFilename(iterationCounter);
    std::ofstream file(_outputVolumeDirectory + "/" + filename + ".csv");
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    };
    writeHeader(file);
    writeData(file);
    file.close();
    std::cout << std::endl;
    std::cout << "Solution written to file: " << filename << std::endl;
    std::cout << std::endl;}


void Output::writeHeader(std::ofstream& file) const {
    file << "NI=" << _ni << "\n";
    file << "NJ=" << _nj << "\n";
    file << "NK=" << _nk << "\n";

    // write coordinates header
    file << "x,y,z";

    // write scalar fields header
    for (auto& field : _outputFields) {
        file << "," << field.first;
    }
    file << "\n";
}

void Output::writeData(std::ofstream& file) const {
    for (size_t i=0; i<_ni; ++i){
        for (size_t j=0; j<_nj; ++j){
            for (size_t k=0; k<_nk; ++k){
                file << _mesh.getVertex(i,j,k).x() 
                     << "," 
                     << _mesh.getVertex(i,j,k).y() 
                     << "," 
                     << _mesh.getVertex(i,j,k).z() ;
                for (auto& field : _outputFields) {
                    file << "," << field.second(i,j,k);
                }
                file << "\n";
            }
        }
    }
}