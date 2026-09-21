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

    FloatType computeTemperature_rho_e(FloatType rho, FloatType e) const override;

    FloatType computeSoundSpeed_rho_e(FloatType rho, FloatType e) const override;

    FloatType computeEntropy_rho_e(FloatType rho, FloatType e) const override;

    FloatType computeDensity_p_s(FloatType p, FloatType s) const override;

    FloatType computeTemperature_p_s(FloatType p, FloatType s) const override;

    FloatType computeInternalEnergy_p_s(FloatType p, FloatType s) const override;

    FloatType computeFundamentalDerivative_rho_e(FloatType rho, FloatType e) const override;

    FloatType computeDpDrho_e(FloatType rho, FloatType e) const override;

    FloatType computeDpDe_rho(FloatType rho, FloatType e) const override;

    FloatType computeInternalEnergy_p_T(FloatType p, FloatType T) const override;

    FloatType computeSoundSpeed_p_T(FloatType p, FloatType T) const override;

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

    ViscosityModel getViscosityModel() const override { return _viscosityModel; }

    void setConstantViscosity(FloatType mu, FloatType cp, FloatType Pr) {
        _viscosityModel = ViscosityModel::CONSTANT;
        _muConstant = mu;
        _cp = cp;
        _Pr = Pr;
    }

    void setSutherlandViscosity(FloatType muRef, FloatType TRef, FloatType S, FloatType cp, FloatType Pr) {
        _viscosityModel = ViscosityModel::SUTHERLAND;
        _sutherlandMuRef = muRef;
        _sutherlandTemperatureRef = TRef;
        _sutherlandSconstant = S;
        _cp = cp;
        _Pr = Pr;
    }

    FloatType computeMolecularDynamicViscosity(FloatType temperature) const override;

    FloatType computeThermalConductivity(FloatType dynamicViscosity) const override;

private:
    FloatType _gamma;  
    FloatType _R;      
    FloatType _cp;     
    FloatType _cv;    

    ViscosityModel _viscosityModel {ViscosityModel::SUTHERLAND};
    FloatType _muConstant {0.0};

    FloatType _sutherlandMuRef {0.0};
    FloatType _sutherlandTemperatureRef {0.0};
    FloatType _sutherlandSconstant {0.0};

    FloatType _Pr {0.72};
};
