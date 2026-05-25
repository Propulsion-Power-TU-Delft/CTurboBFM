#include "output_base.hpp"
#include "math_utils.hpp"


OutputBase::OutputBase(
    const Config &config, 
    const Mesh &mesh, 
    const FlowSolution &solution, 
    const FluidBase &fluid, 
    const Matrix3D<Vector3D> &inviscidForce, 
    const Matrix3D<Vector3D> &viscousForce,
    const Matrix3D<FloatType> &deviationAngle,
    const Matrix3D<FloatType> &wallDistance)
    : _config(config), 
    _mesh(mesh), 
    _solution(solution), 
    _fluid(fluid), 
    _inviscidForce(inviscidForce), 
    _viscousForce(viscousForce), 
    _deviationAngle(deviationAngle),
    _wallDistance(wallDistance),
    _isUnsteadyOutput(_config.saveUnsteadySolution()) {    
    std::filesystem::create_directory(_outputDirectory);
    _outputFields = _config.getOutputFields();
    }



void OutputBase::getOutputFieldsMap(
    std::map<std::string, Matrix3D<FloatType>>& fieldsMap, 
    bool alsoGradients) const {
    
    allocateSpaceForOutput(fieldsMap, alsoGradients);
    storeFields(fieldsMap, alsoGradients);
}

void OutputBase::allocateSpaceForOutput(
    std::map<std::string, Matrix3D<FloatType>>& fieldsMap, 
    bool alsoGradients) const {

    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();
    
    // Primary solution variables
    fieldsMap["Density"] = _solution.getDensity();
    fieldsMap["Velocity X"] = _solution.getVelocityX();
    fieldsMap["Velocity Y"] = _solution.getVelocityY();
    fieldsMap["Velocity Z"] = _solution.getVelocityZ();
    fieldsMap["Total Energy"] = _solution.getTotalEnergy();
    
    // Secondary solution variables
    if (_outputFields == OutputFields::SECONDARY || _outputFields == OutputFields::TURBO_BFM){
        fieldsMap.emplace("Pressure",          Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Temperature",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Mach",           Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Total Pressure",    Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Total Temperature", Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Entropy",           Matrix3D<FloatType>(ni, nj, nk));
    }
    
    // allocate space also for gradients if needed
    if (alsoGradients){
        fieldsMap.emplace("Density Gradient X",          Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Density Gradient Y",          Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Density Gradient Z",          Matrix3D<FloatType>(ni, nj, nk));

        fieldsMap.emplace("Velocity X Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity X Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity X Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        fieldsMap.emplace("Velocity Y Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity Y Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity Y Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        fieldsMap.emplace("Velocity Z Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity Z Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Velocity Z Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        fieldsMap.emplace("Pressure Gradient X",     Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Pressure Gradient Y",     Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Pressure Gradient Z",     Matrix3D<FloatType>(ni, nj, nk));
    }

    const bool isBFMActive = _config.isBFMActive();
    if (isBFMActive && _outputFields == OutputFields::TURBO_BFM){
        fieldsMap["Blockage"] = _mesh.getInputFields(InputField::BLOCKAGE);
        fieldsMap.emplace("Relative Mach",          Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Grid Velocity X",        Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Grid Velocity Y",        Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Grid Velocity Z",        Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Relative Velocity X",    Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Relative Velocity Y",    Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Relative Velocity Z",    Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Viscous Body Force X",   Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Viscous Body Force Y",   Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Viscous Body Force Z",   Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Inviscid Body Force X",  Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Inviscid Body Force Y",  Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Inviscid Body Force Z",  Matrix3D<FloatType>(ni, nj, nk));
        fieldsMap.emplace("Deviation Angle",        Matrix3D<FloatType>(ni, nj, nk));
    }

    if (_config.isTurbulenceActive()){
        fieldsMap.emplace("Wall Distance",        Matrix3D<FloatType>(ni, nj, nk));
    }
}
    

void OutputBase::storeFields(
    std::map<std::string, Matrix3D<FloatType>>& fieldsMap, 
    bool alsoGradients) const {
    
    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();
    const bool isBFMActive = _config.isBFMActive();

    if (_outputFields == OutputFields::SECONDARY || _outputFields == OutputFields::TURBO_BFM){
        Vector3D vel, gridVelCyl, gridVelCart, relVel;
        FloatType rho, et, omega, radius, theta;
        
        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {
                for (size_t k = 0; k < nk; ++k) {
                
                    rho = fieldsMap["Density"](i, j, k);
                    vel.x() = fieldsMap["Velocity X"](i, j, k);
                    vel.y() = fieldsMap["Velocity Y"](i, j, k);
                    vel.z() = fieldsMap["Velocity Z"](i, j, k);
                    et  = fieldsMap["Total Energy"](i, j, k);

                    fieldsMap["Pressure"](i, j, k) = _fluid.computePressure_rho_u_et(rho, vel, et);
                    fieldsMap["Temperature"](i, j, k) = _fluid.computeTemperature_rho_u_et(rho, vel, et);
                    fieldsMap["Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, vel, et);
                    fieldsMap["Total Pressure"](i, j, k) = _fluid.computeTotalPressure_rho_u_et(rho, vel, et);
                    fieldsMap["Total Temperature"](i, j, k) = _fluid.computeTotalTemperature_rho_u_et(rho, vel, et);
                    fieldsMap["Entropy"](i, j, k) = _fluid.computeEntropy_rho_u_et(rho, vel, et);

                    if (isBFMActive && _outputFields == OutputFields::TURBO_BFM){
                        omega = _mesh.getInputFields(InputField::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                        FloatType scalingFactor = _config.getRotationalSpeedScalingFactor();
                        omega *= scalingFactor;
                        radius = std::sqrt(_mesh.getVertex(i, j, k).z() * _mesh.getVertex(i, j, k).z() +
                                                    _mesh.getVertex(i, j, k).y() * _mesh.getVertex(i, j, k).y());
                        theta = std::atan2(_mesh.getVertex(i, j, k).z(), _mesh.getVertex(i, j, k).y());
                        gridVelCyl = {0.0, 0.0, omega * radius};
                        gridVelCart = computeCartesianComponentsFromCylindrical(gridVelCyl, theta);

                        fieldsMap["Grid Velocity X"](i, j, k) = gridVelCart.x();
                        fieldsMap["Grid Velocity Y"](i, j, k) = gridVelCart.y();
                        fieldsMap["Grid Velocity Z"](i, j, k) = gridVelCart.z();
                        
                        relVel = vel - gridVelCart;
                        fieldsMap["Relative Velocity X"](i, j, k) = relVel.x();
                        fieldsMap["Relative Velocity Y"](i, j, k) = relVel.y();
                        fieldsMap["Relative Velocity Z"](i, j, k) = relVel.z();

                        fieldsMap["Relative Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, relVel, et);

                        fieldsMap["Viscous Body Force X"](i, j, k) = _viscousForce(i, j, k).x();
                        fieldsMap["Viscous Body Force Y"](i, j, k) = _viscousForce(i, j, k).y();
                        fieldsMap["Viscous Body Force Z"](i, j, k) = _viscousForce(i, j, k).z();

                        fieldsMap["Inviscid Body Force X"](i, j, k) = _inviscidForce(i, j, k).x();
                        fieldsMap["Inviscid Body Force Y"](i, j, k) = _inviscidForce(i, j, k).y();
                        fieldsMap["Inviscid Body Force Z"](i, j, k) = _inviscidForce(i, j, k).z();

                        fieldsMap["Deviation Angle"](i, j, k) = _deviationAngle(i, j, k);
                    }

                    if (_config.isTurbulenceActive()){
                        fieldsMap["Wall Distance"](i, j, k) = _wallDistance(i, j, k);
                    }

                }
            }
        }
    }
    
    if (alsoGradients){
        // compute gradients
        Matrix3D<Vector3D>  rhoGrad(ni, nj, nk), velXGrad(ni, nj, nk), velYGrad(ni, nj, nk), 
                            velZGrad(ni, nj, nk), pressGrad(ni, nj, nk);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), fieldsMap["Density"], rhoGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), fieldsMap["Velocity X"], velXGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), fieldsMap["Velocity Y"], velYGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), fieldsMap["Velocity Z"], velZGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), fieldsMap["Pressure"], pressGrad);
    


        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {
                for (size_t k = 0; k < nk; ++k) {
                    fieldsMap["Density Gradient X"](i, j, k) = rhoGrad(i, j, k).x();
                    fieldsMap["Density Gradient Y"](i, j, k) = rhoGrad(i, j, k).y();
                    fieldsMap["Density Gradient Z"](i, j, k) = rhoGrad(i, j, k).z();

                    fieldsMap["Velocity X Gradient X"](i, j, k) = velXGrad(i, j, k).x();
                    fieldsMap["Velocity X Gradient Y"](i, j, k) = velXGrad(i, j, k).y();
                    fieldsMap["Velocity X Gradient Z"](i, j, k) = velXGrad(i, j, k).z();

                    fieldsMap["Velocity Y Gradient X"](i, j, k) = velYGrad(i, j, k).x();
                    fieldsMap["Velocity Y Gradient Y"](i, j, k) = velYGrad(i, j, k).y();
                    fieldsMap["Velocity Y Gradient Z"](i, j, k) = velYGrad(i, j, k).z();

                    fieldsMap["Velocity Z Gradient X"](i, j, k) = velZGrad(i, j, k).x();
                    fieldsMap["Velocity Z Gradient Y"](i, j, k) = velZGrad(i, j, k).y();
                    fieldsMap["Velocity Z Gradient Z"](i, j, k) = velZGrad(i, j, k).z();

                    fieldsMap["Pressure Gradient X"](i, j, k) = pressGrad(i, j, k).x();
                    fieldsMap["Pressure Gradient Y"](i, j, k) = pressGrad(i, j, k).y();
                    fieldsMap["Pressure Gradient Z"](i, j, k) = pressGrad(i, j, k).z();
                }
            }
        }
    }

}


std::string OutputBase::getOutputFilename(size_t iterationCounter) {
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