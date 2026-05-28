#pragma once
#include "types.hpp"
#include "fluid_base.hpp"
#include "config.hpp"

// Class for ideal gases
class FluidIdeal : public FluidBase {
public:

    FluidIdeal(FloatType gamma, FloatType R);

    FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const override;

    FloatType computePressure_rho_e(FloatType rho, FloatType e) const override;

    FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const override;
    
    FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const override;

    FloatType computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const override;
    
    FloatType computeStaticEnergy_u_et(FloatType velMag, FloatType et) const override;

    FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const override;
        
    FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const override;
    
    FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    FloatType computeDensity_p_T(FloatType p, FloatType T) const override;

    FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const override;
    
    FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const override;

    FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const override;

    FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const override;

    void computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D flowDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy) override;

    FloatType computePressure_primitive(StateVector primitive) const override;

    FloatType getGamma() const override { return _gamma; }

    FloatType getRconstant() const override { return _R; }

    /** Compute compressor total-to-total efficiency
     * @param pressureRatio Total pressure ratio P0_exit / P0_inlet [-]
     * @param temperatureRatio Total temperature ratio T0_exit / T0_inlet [-]
     * @return Total efficiency [-]
    */
    FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const override;

    Matrix3D<FloatType> computeTemperature_conservative(
        Matrix3D<FloatType>& rho, 
        Matrix3D<FloatType>& ux, 
        Matrix3D<FloatType>& uy, 
        Matrix3D<FloatType>& uz, 
        Matrix3D<FloatType>& et) const override;
    
    void setTransportProperties(const Config &config) override;

    FloatType computeMolecularDynamicViscosity(FloatType temperature) const override;

    FloatType computeThermalConductivity(FloatType dynamicViscosity) const override;

private:
    FloatType _gamma;  
    FloatType _R;      
    FloatType _cp;     
    FloatType _cv;    

    FloatType _sutherlandMuRef;
    FloatType _sutherlandTemperatureRef;
    FloatType _sutherlandSconstant;

    FloatType _Pr;
};
