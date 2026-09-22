#pragma once
#include "types.hpp"
#include "config.hpp"
#include <cmath>

// Base class for fluid thermodynamic models
class FluidBase {
public:

    FluidBase() = default;
    
    virtual ~FluidBase() = default;

    // Thermodynamic evaluations from (rho, e)
    virtual FloatType computePressure_rho_e(FloatType rho, FloatType e) const = 0;

    virtual FloatType computeTemperature_rho_e(FloatType rho, FloatType e) const = 0;

    virtual FloatType computeSoundSpeed_rho_e(FloatType rho, FloatType e) const = 0;

    virtual FloatType computeEntropy_rho_e(FloatType rho, FloatType e) const = 0;

    virtual FloatType computeEnthalpy_rho_e(FloatType rho, FloatType e) const {
        return e + computePressure_rho_e(rho, e) / rho;
    }

    virtual FloatType computeIsentropicExponent_rho_e(FloatType rho, FloatType e) const {
        FloatType p = computePressure_rho_e(rho, e);
        FloatType a = computeSoundSpeed_rho_e(rho, e);
        return rho * a * a / p;
    }

    virtual FloatType computeFundamentalDerivative_rho_e(FloatType rho, FloatType e) const {
        return 0.5 * (getGamma() + 1.0);
    }

    // Thermodynamic derivatives chi = (dp/drho)_e and kappa = (dp/de)_rho
    virtual FloatType computeDpDrho_e(FloatType rho, FloatType e) const = 0;
    virtual FloatType computeDpDe_rho(FloatType rho, FloatType e) const = 0;

    // Thermodynamic evaluations from (p, rho)
    virtual FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const = 0;
    
    virtual FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const = 0;

    virtual FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const = 0;

    // Thermodynamic evaluations from (p, T) and (rho, T)
    virtual FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const = 0;

    virtual FloatType computeDensity_p_T(FloatType p, FloatType T) const = 0;

    virtual FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const = 0;

    virtual FloatType computeInternalEnergy_p_T(FloatType p, FloatType T) const {
        FloatType rho = computeDensity_p_T(p, T);
        return computeStaticEnergy_p_rho(p, rho);
    }

    virtual FloatType computeSoundSpeed_p_T(FloatType p, FloatType T) const {
        FloatType rho = computeDensity_p_T(p, T);
        FloatType e = computeInternalEnergy_p_T(p, T);
        return computeSoundSpeed_rho_e(rho, e);
    }

    // Thermodynamic evaluations from (p, s) - useful for isentropic expansions, stagnation states, outlet BCs
    virtual FloatType computeDensity_p_s(FloatType p, FloatType s) const = 0;

    virtual FloatType computeTemperature_p_s(FloatType p, FloatType s) const = 0;

    virtual FloatType computeInternalEnergy_p_s(FloatType p, FloatType s) const = 0;

    virtual FloatType computeEnthalpy_p_s(FloatType p, FloatType s) const {
        return computeInternalEnergy_p_s(p, s) + p / computeDensity_p_s(p, s);
    }

    // Kinematics and state vector conversions
    virtual FloatType computeStaticEnergy_u_et(const Vector3D& vel, FloatType et) const {
        return et - 0.5 * vel.dot(vel);
    }
    
    virtual FloatType computeStaticEnergy_u_et(FloatType velMag, FloatType et) const {
        return et - 0.5 * velMag * velMag;
    }

    virtual FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType e = computeStaticEnergy_u_et(u, et);
        return computeSoundSpeed_rho_e(rho, e);
    }
    
    virtual FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType e = computeStaticEnergy_u_et(u, et);
        return computePressure_rho_e(rho, e);
    }

    virtual FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType e = computeStaticEnergy_u_et(u, et);
        return computeTemperature_rho_e(rho, e);
    }

    virtual FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType e = computeStaticEnergy_u_et(u, et);
        return computeEntropy_rho_e(rho, e);
    }

    virtual FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType p = computePressure_rho_u_et(rho, u, et);
        return et + p / rho;
    }

    virtual FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
        FloatType a = computeSoundSpeed_rho_u_et(rho, u, et);
        return u.magnitude() / a;
    }

    virtual FloatType computePressure_primitive(StateVector primitive) const {
        Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
        return computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    }

    // Stagnation and isentropic relations
    virtual FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    virtual FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const = 0;

    virtual FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const = 0;

    virtual FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const = 0;

    virtual FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const = 0;

    virtual void computeInitFields(
        FloatType initMach, 
        FloatType initTemperature, 
        FloatType initPressure, 
        Vector3D flowDirection, 
        FloatType &density, 
        Vector3D &velocity, 
        FloatType &totEnergy) = 0;

    virtual FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const = 0;

    virtual FloatType computeTotalEfficiency_h(FloatType ht_in, FloatType ht_out, FloatType ht_out_s) const {
        FloatType delta_h = ht_out - ht_in;
        FloatType delta_h_s = ht_out_s - ht_in;
        if (std::abs(delta_h) < 1e-12) return 1.0;
        if (delta_h >= 0.0) {
            // Compressor / pump: isentropic work / actual work
            return delta_h_s / delta_h;
        } else {
            // Turbine / expander: actual work / isentropic work
            if (std::abs(delta_h_s) < 1e-12) return 1.0;
            return delta_h / delta_h_s;
        }
    }

    virtual FloatType getGamma() const = 0;

    virtual FloatType getRconstant() const = 0;

    virtual Matrix3D<FloatType> computeTemperature_conservative(
        Matrix3D<FloatType>& rho, 
        Matrix3D<FloatType>& ux, 
        Matrix3D<FloatType>& uy, 
        Matrix3D<FloatType>& uz, 
        Matrix3D<FloatType>& et) const {
        Matrix3D<FloatType> temperature(rho.sizeI(), rho.sizeJ(), rho.sizeK());
        for (size_t i = 0; i < rho.sizeI(); ++i) {
            for (size_t j = 0; j < rho.sizeJ(); ++j) {
                for (size_t k = 0; k < rho.sizeK(); ++k) {
                    FloatType e = et(i,j,k) - 0.5 * (ux(i,j,k)*ux(i,j,k) + uy(i,j,k)*uy(i,j,k) + uz(i,j,k)*uz(i,j,k));
                    temperature(i,j,k) = computeTemperature_rho_e(rho(i,j,k), e);
                }
            }
        }
        return temperature;
    }
    
    virtual void setTransportProperties(const Config &config) = 0;

    virtual ViscosityModel getViscosityModel() const { return ViscosityModel::SUTHERLAND; }

    virtual FloatType computeMolecularDynamicViscosity(FloatType temperature) const = 0;

    virtual FloatType computeMolecularDynamicViscosity_rho_e(FloatType rho, FloatType e) const {
        FloatType T = computeTemperature_rho_e(rho, e);
        return computeMolecularDynamicViscosity(T);
    }

    virtual FloatType computeThermalConductivity(FloatType dynamicViscosity) const = 0;
    
    // Bounds checking
    virtual bool isStateInBounds_rho_e(FloatType rho, FloatType e) const { return true; }
    virtual FloatType getRhoMin() const { return 0.0; }
    virtual FloatType getRhoMax() const { return 1e9; }
    virtual FloatType getEMin() const { return 0.0; }
    virtual FloatType getEMax() const { return 1e9; }
};

