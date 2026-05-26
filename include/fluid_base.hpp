#pragma once
#include "types.hpp"
#include "config.hpp"

// Class for ideal gases
class FluidBase {
public:

    FluidBase() = default;
    
    virtual ~FluidBase() = default;

    virtual FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const = 0;
    
    virtual FloatType computePressure_rho_e(FloatType rho, FloatType e) const = 0;

    virtual FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const = 0;
    
    virtual FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const = 0;

    virtual FloatType computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const = 0;
    
    virtual FloatType computeStaticEnergy_u_et(FloatType velMag, FloatType et) const = 0;

    virtual FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    virtual FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const = 0;

    virtual FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const = 0;

    virtual FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeDensity_p_T(FloatType p, FloatType T) const = 0;

    virtual FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const = 0;

    virtual FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const = 0;

    virtual FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const = 0;

    virtual FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const = 0;

    virtual void computeInitFields(
        FloatType initMach, 
        FloatType initTemperature, 
        FloatType initPressure, 
        Vector3D flowDirection, 
        FloatType &density, 
        Vector3D &velocity, 
        FloatType &totEnergy) = 0;

    virtual FloatType computePressure_primitive(StateVector primitive) const = 0;

    virtual FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const = 0;

    virtual FloatType getGamma() const = 0;

    virtual FloatType getRconstant() const = 0;

    virtual Matrix3D<FloatType> computeTemperature_conservative(
        Matrix3D<FloatType>& rho, 
        Matrix3D<FloatType>& ux, 
        Matrix3D<FloatType>& uy, 
        Matrix3D<FloatType>& uz, 
        Matrix3D<FloatType>& et) const = 0;
    
    virtual void setTransportProperties(const Config &config) = 0;

    virtual FloatType computeMolecularDynamicViscosity(FloatType temperature) = 0;

    virtual FloatType computeThermalConductivity(FloatType dynamicViscosity) = 0;
    
};
