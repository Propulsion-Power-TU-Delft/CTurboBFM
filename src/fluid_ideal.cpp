#include "fluid_ideal.hpp"
#include <cmath>

FluidIdeal::FluidIdeal(FloatType gamma, FloatType R)
    : _gamma(gamma), _R(R) {
    _cp = (_gamma * _R) / (_gamma - 1);
    _cv = _cp / _gamma;
}

FloatType FluidIdeal::computeStaticEnergy_p_rho(FloatType p, FloatType rho) const {
    return p / ((_gamma - 1) * rho);
}

FloatType FluidIdeal::computePressure_rho_e(FloatType rho, FloatType e) const {
    return (_gamma - 1) * rho * e;
}

FloatType FluidIdeal::computePressure_rho_T(FloatType rho, FloatType Temp) const {
    return rho * _R * Temp;
}

FloatType FluidIdeal::computeSoundSpeed_p_rho(FloatType p, FloatType rho) const {
    return std::sqrt(_gamma * p / rho);
}

FloatType FluidIdeal::computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const {
    FloatType velMag2 = std::pow(vel.magnitude(), 2);
    return et - 0.5 * velMag2;
}

FloatType FluidIdeal::computeStaticEnergy_u_et(FloatType velMag, FloatType et) const {
    return et - 0.5 * velMag * velMag;
}

FloatType FluidIdeal::computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType p = computePressure_rho_e(rho, e);
    return computeSoundSpeed_p_rho(p, rho);
}

FloatType FluidIdeal::computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    return computePressure_rho_e(rho, e);
}

FloatType FluidIdeal::computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    FloatType M = computeMachNumber_rho_u_et(rho, u, et);
    return p * std::pow(1.0 + (_gamma - 1.0) / 2.0 * M * M, _gamma / (_gamma - 1));
}

FloatType FluidIdeal::computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType T = computeTemperature_rho_u_et(rho, u, et);
    FloatType M = computeMachNumber_rho_u_et(rho, u, et);
    return T * (1.0 + (_gamma - 1.0) / 2.0 * M * M);
}

FloatType FluidIdeal::computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    return p / (rho * _R);
}

FloatType FluidIdeal::computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType p = computePressure_rho_e(rho, e);
    return et + p / rho;
}

FloatType FluidIdeal::computeStaticPressure_pt_M(FloatType pt, FloatType M) const {
    return pt * std::pow(1.0 + (_gamma - 1.0) / 2.0 * M * M, -_gamma / (_gamma - 1.0));
}

FloatType FluidIdeal::computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const {
    return Tt / (1 + (_gamma - 1) / 2 * M * M);
}

FloatType FluidIdeal::computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType p = computePressure_rho_u_et(rho, u, et);
    // FloatType s = computeEntropy_p_rho(p, rho);

    FloatType T = computeTemperature_rho_u_et(rho, u, et);
    FloatType s = computeEntropy_p_T(p, T);
    return s;
}

FloatType FluidIdeal::computeDensity_p_T(FloatType p, FloatType T) const {
    return p / (_R * T);
}

FloatType FluidIdeal::computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType a = computeSoundSpeed_rho_u_et(rho, u, et);
    FloatType umag = u.magnitude();
    return umag / a;
}


FloatType FluidIdeal::computeTotalPressure_p_M(FloatType pressure, FloatType mach) const {
    return pressure * std::pow(1.0 + (_gamma - 1.0) / 2.0 * mach * mach, _gamma / (_gamma - 1));
}

FloatType FluidIdeal::computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const {
    return temperature * (1.0 + (_gamma - 1.0) / 2.0 * mach * mach);
}

FloatType FluidIdeal::computeEntropy_p_rho(FloatType pressure, FloatType density) const {
    return pressure / std::pow(density, _gamma);
}

void FluidIdeal::computeInitFields(
    FloatType initMach, 
    FloatType initTemperature, 
    FloatType initPressure, 
    Vector3D flowDirection, 
    FloatType &density, 
    Vector3D &velocity, 
    FloatType &totEnergy){

    FloatType soundSpeed = std::sqrt(_gamma*_R*initTemperature);
    velocity = (flowDirection / flowDirection.magnitude()) * soundSpeed * initMach;
    density = computeDensity_p_T(initPressure, initTemperature);
    FloatType energy = computeStaticEnergy_p_rho(initPressure, density);
    totEnergy = energy + 0.5 * pow(velocity.magnitude(), 2);
}

FloatType FluidIdeal::computePressure_primitive(StateVector primitive) const {
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto pressure = computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    return pressure;
}

FloatType FluidIdeal::computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const{
    FloatType eta = (std::pow(pressureRatio, (_gamma - 1.0) / _gamma) - 1.0) / (temperatureRatio - 1.0);
    return eta;
}

FloatType FluidIdeal::computeEntropy_p_T(FloatType pressure, FloatType temperature) const{
    FloatType entropy = _cp*log(temperature/288.15) - _R*log(pressure/101325);
    return entropy;
}

Matrix3D<FloatType> FluidIdeal::computeTemperature_conservative(
    Matrix3D<FloatType>& rho, 
    Matrix3D<FloatType>& ux, 
    Matrix3D<FloatType>& uy, 
    Matrix3D<FloatType>& uz, 
    Matrix3D<FloatType>& et) const {
    
    Matrix3D<FloatType> temperature = (et - (ux * ux + uy * uy + uz * uz) * 0.5) * ((_gamma-1.0)/_R);
    return temperature;
}

void FluidIdeal::setTransportProperties(const Config &config){
    _sutherlandMuRef = config.getSutherlandMuRef();
    _sutherlandTemperatureRef = config.getSutherlandTemperatureRef();
    _sutherlandSconstant = config.getSutherlandSconstant();

    _cp = config.getFluidHeatCapacity(); 
    _Pr = config.getFluidPrandtlNumber();      
}

FloatType FluidIdeal::computeMolecularDynamicViscosity(FloatType temperature) const {
    return _sutherlandMuRef * std::pow(temperature/_sutherlandTemperatureRef, 1.5) * (
        (_sutherlandTemperatureRef + _sutherlandSconstant) / (temperature + _sutherlandSconstant)
    );   
}

FloatType FluidIdeal::computeThermalConductivity(FloatType dynamicViscosity) const {
    return _cp * dynamicViscosity / _Pr; 
}