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
    constexpr const char* EDDY_VISCOSITY            = "Eddy Viscosity";
    constexpr const char* MOLECULAR_VISCOSITY       = "Molecular Viscosity";
    constexpr const char* WALL_SHEAR_STRESS_X       = "Wall Shear Stress X";
    constexpr const char* WALL_SHEAR_STRESS_Y       = "Wall Shear Stress Y";
    constexpr const char* WALL_SHEAR_STRESS_Z       = "Wall Shear Stress Z";
    constexpr const char* Y_PLUS                    = "Y Plus";
}

Output::Output(
    const Config &config, 
    const Mesh &mesh, 
    const FlowSolution &solution, 
    const std::map<SolutionName, Matrix3D<Vector3D>>& solutionGrad,
    const FluidBase &fluid, 
    const std::vector<Boundary> &boundaries,
    const TurbulenceModelBase &turbulenceModel,
    const Matrix3D<Vector3D> &inviscidForce, 
    const Matrix3D<Vector3D> &viscousForce,
    const Matrix3D<FloatType> &deviationAngle,
    const Matrix3D<FloatType> &wallDistance)
    : _config(config), 
    _mesh(mesh), 
    _solution(solution), 
    _solutionGrad(solutionGrad),
    _fluid(fluid), 
    _boundaries(boundaries),
    _turbulenceModel(turbulenceModel),
    _surfacesI(mesh.getSurfacesI()),
    _surfacesJ(mesh.getSurfacesJ()),
    _surfacesK(mesh.getSurfacesK()),
    _inviscidForce(inviscidForce), 
    _viscousForce(viscousForce), 
    _deviationAngle(deviationAngle),
    _wallDistance(wallDistance) {    

    std::filesystem::create_directory(_outputVolumeDirectory);
    _outputFieldsType = _config.getOutputFieldsType();
    _isUnsteadyOutput = _config.saveUnsteadySolution();
    
    _ni = _mesh.getNumberPointsI();
    _nj = _mesh.getNumberPointsJ();
    _nk = _mesh.getNumberPointsK();
    
    _isBfmActive = _config.isBFMActive();
    _isTurbulenceActive = _config.isTurbulenceActive();
    _isViscosityActive = _config.isViscosityActive();
    
    allocateOutputArrays();
    }

void Output::allocateOutputArrays() {
    
    // Primary solution variables
    _outputFields.emplace(DENSITY,                      Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_X,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_Y,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(VELOCITY_Z,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(TOTAL_ENERGY,                 Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(PRESSURE,                     Matrix3D<FloatType>(_ni, _nj, _nk));
    _outputFields.emplace(TEMPERATURE,                  Matrix3D<FloatType>(_ni, _nj, _nk));
    
    // Secondary solution variables
    if (_outputFieldsType == OutputFieldsType::SECONDARY || _outputFieldsType == OutputFieldsType::TURBO_BFM){
        _outputFields.emplace(MACH,                     Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(TOTAL_PRESSURE,           Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(TOTAL_TEMPERATURE,        Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(ENTROPY,                  Matrix3D<FloatType>(_ni, _nj, _nk));
    }

    // Turbo BFM variables
    if (_isBfmActive && _outputFieldsType == OutputFieldsType::TURBO_BFM){
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

    // RANS variables
    if (_isTurbulenceActive || _isViscosityActive){
        _outputFields.emplace(EDDY_VISCOSITY,           Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(MOLECULAR_VISCOSITY,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(WALL_SHEAR_STRESS_X,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(WALL_SHEAR_STRESS_Y,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(WALL_SHEAR_STRESS_Z,      Matrix3D<FloatType>(_ni, _nj, _nk));
        _outputFields.emplace(Y_PLUS,                   Matrix3D<FloatType>(_ni, _nj, _nk));
    }
}
    

void Output::updateOutputArrays() {
    
    // Primary
    updatePrimaryFields();

    // Secondary
    if (_outputFieldsType != OutputFieldsType::PRIMARY || _isTurbulenceActive){
        updateSecondaryFields();
    }

    // Turbulence
    if (_isTurbulenceActive || _isViscosityActive){
        updateViscousFields();
    }

}

void Output::updatePrimaryFields() {
    _outputFields[DENSITY] = _solution.getDensity();
    _outputFields[VELOCITY_X] = _solution.getVelocityX();
    _outputFields[VELOCITY_Y] = _solution.getVelocityY();
    _outputFields[VELOCITY_Z] = _solution.getVelocityZ();
    _outputFields[TOTAL_ENERGY] = _solution.getTotalEnergy();
    
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                _outputFields[PRESSURE](i, j, k) = _fluid.computePressure_rho_u_et(
                    _outputFields[DENSITY](i, j, k),
                    {_outputFields[VELOCITY_X](i, j, k),
                    _outputFields[VELOCITY_Y](i, j, k),
                    _outputFields[VELOCITY_Z](i, j, k)},
                    _outputFields[TOTAL_ENERGY](i, j, k)
                );
                
                _outputFields[TEMPERATURE](i, j, k) = _fluid.computeTemperature_rho_u_et(
                    _outputFields[DENSITY](i, j, k),
                    {_outputFields[VELOCITY_X](i, j, k),
                    _outputFields[VELOCITY_Y](i, j, k),
                    _outputFields[VELOCITY_Z](i, j, k)},
                    _outputFields[TOTAL_ENERGY](i, j, k)
                );
            }
        }
    }
}

void Output::updateSecondaryFields() {
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

                _outputFields[MACH](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, vel, et);
                _outputFields[TOTAL_PRESSURE](i, j, k) = _fluid.computeTotalPressure_rho_u_et(rho, vel, et);
                _outputFields[TOTAL_TEMPERATURE](i, j, k) = _fluid.computeTotalTemperature_rho_u_et(rho, vel, et);
                _outputFields[ENTROPY](i, j, k) = _fluid.computeEntropy_rho_u_et(rho, vel, et);
                
                // only for BFM
                if (_isBfmActive && _outputFieldsType == OutputFieldsType::TURBO_BFM){

                    _outputFields[BLOCKAGE](i, j, k) = _mesh.getInputFields(InputField::BLOCKAGE)(i, j, k);
                    omega = _mesh.getInputFields(InputField::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                    omega *= _config.getRotationalSpeedScalingFactor();
                    radius = _mesh.getRadius(i, j, k);
                    theta = _mesh.getTheta(i, j, k);
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

void Output::updateViscousFields() {

    // volumetric fields
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                _outputFields[MOLECULAR_VISCOSITY](i, j, k) = _fluid.computeMolecularDynamicViscosity(
                    _outputFields[TEMPERATURE](i, j, k)
                );

                _outputFields[EDDY_VISCOSITY](i, j, k) = _turbulenceModel.getEddyViscosity(_outputFields[DENSITY](i, j, k), i, j, k);
            }
        }
    }

    // no-slip wall fields
    for (auto& bound: _boundaries) {
        if (bound.type == BoundaryType::NO_SLIP_WALL) {
            BoundaryNodesIndexRange range = fetchBoundaryNodesIndexRange(bound, _ni, _nj, _nk);
            BoundaryOrientation orientation = bound.orientation;
            std::array<int, 3> step = {0, 0, 0};
            if (orientation == BoundaryOrientation::I_START)
            {
                step[0] = 1;
            }
            else if (orientation == BoundaryOrientation::I_END)
            {
                step[0] = -1;
            }
            else if (orientation == BoundaryOrientation::J_START)
            {
                step[1] = 1;
            }
            else if (orientation == BoundaryOrientation::J_END)
            {
                step[1] = -1;
            }
            else if (orientation == BoundaryOrientation::K_START)
            {
                step[2] = 1;
            }
            else if (orientation == BoundaryOrientation::K_END)
            {
                step[2] = -1;
            }

            FloatType d = 0.0;
            for (size_t i = range.iStart; i < range.iLast; ++i) {
                for (size_t j = range.jStart; j < range.jLast; ++j) {
                    for (size_t k = range.kStart; k < range.kLast; ++k) {

                        // the node is already on the boundary, no need for interpolation
                        Vector3D tauWall = computeWallShearStress(bound, i, j, k);
                        _outputFields[WALL_SHEAR_STRESS_X](i, j, k) = tauWall.x();
                        _outputFields[WALL_SHEAR_STRESS_Y](i, j, k) = tauWall.y();
                        _outputFields[WALL_SHEAR_STRESS_Z](i, j, k) = tauWall.z();

                        d = _wallDistance(i+step[0], j+step[1], k+step[2]) - _wallDistance(i, j, k);
                        _outputFields[Y_PLUS](i, j, k) = std::sqrt(
                            tauWall.magnitude() * _outputFields[DENSITY](i, j, k)
                            ) * d / _outputFields[MOLECULAR_VISCOSITY](i, j, k);
                    }
                }
            }
        }
    }
}

Vector3D Output::computeWallShearStress(Boundary boundary, size_t i, size_t j, size_t k) {
    
    // velocity gradient
    Vector3D gradU, gradV, gradW;
    gradU = _solutionGrad.at(SolutionName::VELOCITY_X)(i, j, k);
    gradV = _solutionGrad.at(SolutionName::VELOCITY_Y)(i, j, k);
    gradW = _solutionGrad.at(SolutionName::VELOCITY_Z)(i, j, k);

    // viscous stress with Stokes hypothesis
    FloatType mu = _outputFields[MOLECULAR_VISCOSITY](i, j, k);
    FloatType muEddy = _outputFields[EDDY_VISCOSITY](i, j, k);
    FloatType muEffective = mu + muEddy;
    FloatType lambda = -2.0 * mu / 3.0;
    ViscousStressTensor tau = computeViscousStressTensor(
        muEffective,
        lambda,
        gradU,
        gradV,
        gradW
    );

    // Extract the surface vector, oriented towards the fluid
    Vector3D surface;
    switch (boundary.orientation) {
        case BoundaryOrientation::I_START : {
            surface = _surfacesI(0, j, k);
            break;
        }
        case BoundaryOrientation::I_END : {
            surface = -_surfacesI(boundary.i_max, j, k);
            break;
        }
        case BoundaryOrientation::J_START : {
            surface = _surfacesJ(i, 0, k);
            break;
        }
        case BoundaryOrientation::J_END : {
            surface = -_surfacesJ(i, boundary.j_max, k);
            break;
        }
        case BoundaryOrientation::K_START : {
            surface = _surfacesK(i, j, 0);
            break;
        }
        case BoundaryOrientation::K_END : {
            surface = -_surfacesK(i, j, boundary.k_max);
            break;
        }
    }
    surface = surface.normalized();
    Vector3D tauWall = computeViscousStressVector(tau, surface);

    return tauWall;
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