#include "fluid_real.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>

FluidReal::FluidReal(const std::string& tablePath) {
    loadTable(tablePath);
}

void FluidReal::loadTable(const std::string& tablePath) {
    std::ifstream file(tablePath);
    if (!file.is_open()) {
        throw std::runtime_error("FluidReal: Could not open table file at: " + tablePath);
    }

    std::string line;
    std::string currentSection = "";

    // Variables for RHO_E section
    size_t nRho_re = 0, nE_re = 0;
    FloatType rhoMin_re = 0, rhoMax_re = 0, eMin_re = 0, eMax_re = 0;

    // Variables for P_T section
    size_t nP_pt = 0, nT_pt = 0;
    FloatType pMin_pt = 0, pMax_pt = 0, TMin_pt = 0, TMax_pt = 0;

    // Variables for P_S section
    size_t nP_ps = 0, nS_ps = 0;
    FloatType pMin_ps = 0, pMax_ps = 0, sMin_ps = 0, sMax_ps = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '\r') continue;

        // Check for section markers
        if (line.find("# SECTION: RHO_E") != std::string::npos) {
            currentSection = "RHO_E";
            continue;
        } else if (line.find("# SECTION: P_T") != std::string::npos) {
            currentSection = "P_T";
            continue;
        } else if (line.find("# SECTION: P_S") != std::string::npos) {
            currentSection = "P_S";
            continue;
        }

        // Skip comment lines
        if (line[0] == '#') continue;

        // Header metadata
        if (currentSection.empty()) {
            auto colonPos = line.find(':');
            if (colonPos != std::string::npos) {
                std::string key = line.substr(0, colonPos);
                std::string val = line.substr(colonPos + 1);
                // trim whitespace
                key.erase(key.find_last_not_of(" \t\r\n") + 1);
                val.erase(0, val.find_first_not_of(" \t\r\n"));

                if (key == "FLUID_NAME") _fluidName = val;
                else if (key == "MOLAR_MASS") _molarMass = std::stod(val);
                else if (key == "GAS_CONSTANT") _gasConstant = std::stod(val);
                else if (key == "CRITICAL_PRESSURE") _pCrit = std::stod(val);
                else if (key == "CRITICAL_TEMPERATURE") _TCrit = std::stod(val);
                else if (key == "CRITICAL_DENSITY") _rhoCrit = std::stod(val);
                else if (key == "NOMINAL_GAMMA") _nominalGamma = std::stod(val);
            }
            continue;
        }

        // Section RHO_E metadata and data
        if (currentSection == "RHO_E") {
            auto colonPos = line.find(':');
            if (colonPos != std::string::npos) {
                std::string key = line.substr(0, colonPos);
                std::string val = line.substr(colonPos + 1);
                key.erase(key.find_last_not_of(" \t\r\n") + 1);
                val.erase(0, val.find_first_not_of(" \t\r\n"));

                if (key == "N_RHO") nRho_re = std::stoul(val);
                else if (key == "N_E") nE_re = std::stoul(val);
                else if (key == "RHO_MIN") rhoMin_re = std::stod(val);
                else if (key == "RHO_MAX") rhoMax_re = std::stod(val);
                else if (key == "E_MIN") eMin_re = std::stod(val);
                else if (key == "E_MAX") {
                    eMax_re = std::stod(val);

                    // Initialize all RHO_E tables
                    _table_p.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_T.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_a.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_s.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_mu.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_kappa.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_dp_drho.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);
                    _table_dp_de.init(nRho_re, nE_re, rhoMin_re, rhoMax_re, eMin_re, eMax_re);

                    // Read RHO_E data
                    for (size_t i = 0; i < nRho_re; ++i) {
                        for (size_t j = 0; j < nE_re; ++j) {
                            if (!std::getline(file, line)) break;
                            while (line.empty() || line[0] == '#' || line[0] == '\r') {
                                if (!std::getline(file, line)) break;
                            }
                            std::istringstream iss(line);
                            FloatType p, T, a, s, mu, kappa, dp_drho, dp_de;
                            iss >> p >> T >> a >> s >> mu >> kappa >> dp_drho >> dp_de;
                            size_t idx = i * nE_re + j;
                            _table_p.data[idx] = p;
                            _table_T.data[idx] = T;
                            _table_a.data[idx] = a;
                            _table_s.data[idx] = s;
                            _table_mu.data[idx] = mu;
                            _table_kappa.data[idx] = kappa;
                            _table_dp_drho.data[idx] = dp_drho;
                            _table_dp_de.data[idx] = dp_de;
                        }
                    }
                }
            }
            continue;
        }

        // Section P_T metadata and data
        if (currentSection == "P_T") {
            auto colonPos = line.find(':');
            if (colonPos != std::string::npos) {
                std::string key = line.substr(0, colonPos);
                std::string val = line.substr(colonPos + 1);
                key.erase(key.find_last_not_of(" \t\r\n") + 1);
                val.erase(0, val.find_first_not_of(" \t\r\n"));

                if (key == "N_P") nP_pt = std::stoul(val);
                else if (key == "N_T") nT_pt = std::stoul(val);
                else if (key == "P_MIN") pMin_pt = std::stod(val);
                else if (key == "P_MAX") pMax_pt = std::stod(val);
                else if (key == "T_MIN") TMin_pt = std::stod(val);
                else if (key == "T_MAX") {
                    TMax_pt = std::stod(val);

                    _table_rho_pt.init(nP_pt, nT_pt, pMin_pt, pMax_pt, TMin_pt, TMax_pt);
                    _table_e_pt.init(nP_pt, nT_pt, pMin_pt, pMax_pt, TMin_pt, TMax_pt);
                    _table_a_pt.init(nP_pt, nT_pt, pMin_pt, pMax_pt, TMin_pt, TMax_pt);
                    _table_s_pt.init(nP_pt, nT_pt, pMin_pt, pMax_pt, TMin_pt, TMax_pt);

                    for (size_t i = 0; i < nP_pt; ++i) {
                        for (size_t j = 0; j < nT_pt; ++j) {
                            if (!std::getline(file, line)) break;
                            while (line.empty() || line[0] == '#' || line[0] == '\r') {
                                if (!std::getline(file, line)) break;
                            }
                            std::istringstream iss(line);
                            FloatType rho, e, a, s;
                            iss >> rho >> e >> a >> s;
                            size_t idx = i * nT_pt + j;
                            _table_rho_pt.data[idx] = rho;
                            _table_e_pt.data[idx] = e;
                            _table_a_pt.data[idx] = a;
                            _table_s_pt.data[idx] = s;
                        }
                    }
                }
            }
            continue;
        }

        // Section P_S metadata and data
        if (currentSection == "P_S") {
            auto colonPos = line.find(':');
            if (colonPos != std::string::npos) {
                std::string key = line.substr(0, colonPos);
                std::string val = line.substr(colonPos + 1);
                key.erase(key.find_last_not_of(" \t\r\n") + 1);
                val.erase(0, val.find_first_not_of(" \t\r\n"));

                if (key == "N_P") nP_ps = std::stoul(val);
                else if (key == "N_S") nS_ps = std::stoul(val);
                else if (key == "P_MIN") pMin_ps = std::stod(val);
                else if (key == "P_MAX") pMax_ps = std::stod(val);
                else if (key == "S_MIN") sMin_ps = std::stod(val);
                else if (key == "S_MAX") {
                    sMax_ps = std::stod(val);

                    _table_rho_ps.init(nP_ps, nS_ps, pMin_ps, pMax_ps, sMin_ps, sMax_ps);
                    _table_T_ps.init(nP_ps, nS_ps, pMin_ps, pMax_ps, sMin_ps, sMax_ps);
                    _table_e_ps.init(nP_ps, nS_ps, pMin_ps, pMax_ps, sMin_ps, sMax_ps);

                    for (size_t i = 0; i < nP_ps; ++i) {
                        for (size_t j = 0; j < nS_ps; ++j) {
                            if (!std::getline(file, line)) break;
                            while (line.empty() || line[0] == '#' || line[0] == '\r') {
                                if (!std::getline(file, line)) break;
                            }
                            std::istringstream iss(line);
                            FloatType rho, T, e;
                            iss >> rho >> T >> e;
                            size_t idx = i * nS_ps + j;
                            _table_rho_ps.data[idx] = rho;
                            _table_T_ps.data[idx] = T;
                            _table_e_ps.data[idx] = e;
                        }
                    }
                }
            }
            continue;
        }
    }
}

FloatType FluidReal::computePressure_rho_e(FloatType rho, FloatType e) const {
    return _table_p.interpolate(rho, e);
}

FloatType FluidReal::computeTemperature_rho_e(FloatType rho, FloatType e) const {
    return _table_T.interpolate(rho, e);
}

FloatType FluidReal::computeSoundSpeed_rho_e(FloatType rho, FloatType e) const {
    return _table_a.interpolate(rho, e);
}

FloatType FluidReal::computeEntropy_rho_e(FloatType rho, FloatType e) const {
    return _table_s.interpolate(rho, e);
}

FloatType FluidReal::computeStaticEnergy_p_rho(FloatType p, FloatType rho) const {
    // 1D bisection on e such that computePressure_rho_e(rho, e) == p
    FloatType e_low = _table_p.yMin;
    FloatType e_high = _table_p.yMax;
    FloatType p_low = computePressure_rho_e(rho, e_low);
    FloatType p_high = computePressure_rho_e(rho, e_high);

    if (p <= p_low) return e_low;
    if (p >= p_high) return e_high;

    for (int iter = 0; iter < 30; ++iter) {
        FloatType e_mid = 0.5 * (e_low + e_high);
        FloatType p_mid = computePressure_rho_e(rho, e_mid);
        if (std::abs(p_mid - p) <= 1e-4 * p) {
            return e_mid;
        }
        if (p_mid < p) {
            e_low = e_mid;
        } else {
            e_high = e_mid;
        }
    }
    return 0.5 * (e_low + e_high);
}

FloatType FluidReal::computeSoundSpeed_p_rho(FloatType p, FloatType rho) const {
    FloatType e = computeStaticEnergy_p_rho(p, rho);
    return computeSoundSpeed_rho_e(rho, e);
}

FloatType FluidReal::computeEntropy_p_rho(FloatType pressure, FloatType density) const {
    FloatType e = computeStaticEnergy_p_rho(pressure, density);
    return computeEntropy_rho_e(density, e);
}

FloatType FluidReal::computePressure_rho_T(FloatType rho, FloatType Temp) const {
    // 1D bisection on e such that computeTemperature_rho_e(rho, e) == Temp
    FloatType e_low = _table_T.yMin;
    FloatType e_high = _table_T.yMax;
    for (int iter = 0; iter < 25; ++iter) {
        FloatType e_mid = 0.5 * (e_low + e_high);
        FloatType T_mid = computeTemperature_rho_e(rho, e_mid);
        if (T_mid < Temp) {
            e_low = e_mid;
        } else {
            e_high = e_mid;
        }
    }
    FloatType e = 0.5 * (e_low + e_high);
    return computePressure_rho_e(rho, e);
}

FloatType FluidReal::computeDensity_p_T(FloatType p, FloatType T) const {
    return _table_rho_pt.interpolate(p, T);
}

FloatType FluidReal::computeEntropy_p_T(FloatType pressure, FloatType temperature) const {
    return _table_s_pt.interpolate(pressure, temperature);
}

FloatType FluidReal::computeDensity_p_s(FloatType p, FloatType s) const {
    return _table_rho_ps.interpolate(p, s);
}

FloatType FluidReal::computeTemperature_p_s(FloatType p, FloatType s) const {
    return _table_T_ps.interpolate(p, s);
}

FloatType FluidReal::computeInternalEnergy_p_s(FloatType p, FloatType s) const {
    return _table_e_ps.interpolate(p, s);
}

FloatType FluidReal::computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType p = computePressure_rho_e(rho, e);
    FloatType s = computeEntropy_rho_e(rho, e);
    FloatType ht = et + p / rho;

    // Stagnation state: find p_t >= p such that h(p_t, s) == ht
    FloatType p_low = p;
    FloatType p_high = _table_rho_ps.xMax;

    for (int iter = 0; iter < 25; ++iter) {
        FloatType p_mid = 0.5 * (p_low + p_high);
        FloatType rho_mid = computeDensity_p_s(p_mid, s);
        FloatType e_mid = computeInternalEnergy_p_s(p_mid, s);
        FloatType h_mid = e_mid + p_mid / rho_mid;
        if (h_mid < ht) {
            p_low = p_mid;
        } else {
            p_high = p_mid;
        }
    }
    return 0.5 * (p_low + p_high);
}

FloatType FluidReal::computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const {
    FloatType pt = computeTotalPressure_rho_u_et(rho, u, et);
    FloatType e = computeStaticEnergy_u_et(u, et);
    FloatType s = computeEntropy_rho_e(rho, e);
    return computeTemperature_p_s(pt, s);
}

FloatType FluidReal::computeStaticPressure_pt_M(FloatType pt, FloatType M) const {
    return pt * std::pow(1.0 + (_nominalGamma - 1.0) / 2.0 * M * M, -_nominalGamma / (_nominalGamma - 1.0));
}

FloatType FluidReal::computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const {
    return Tt / (1.0 + (_nominalGamma - 1.0) / 2.0 * M * M);
}

FloatType FluidReal::computeTotalPressure_p_M(FloatType pressure, FloatType mach) const {
    return pressure * std::pow(1.0 + (_nominalGamma - 1.0) / 2.0 * mach * mach, _nominalGamma / (_nominalGamma - 1.0));
}

FloatType FluidReal::computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const {
    return temperature * (1.0 + (_nominalGamma - 1.0) / 2.0 * mach * mach);
}

void FluidReal::computeInitFields(
    FloatType initMach, 
    FloatType initTemperature, 
    FloatType initPressure, 
    Vector3D flowDirection, 
    FloatType &density, 
    Vector3D &velocity, 
    FloatType &totEnergy) {

    density = computeDensity_p_T(initPressure, initTemperature);
    FloatType energy = _table_e_pt.interpolate(initPressure, initTemperature);
    FloatType soundSpeed = computeSoundSpeed_rho_e(density, energy);
    FloatType dirMag = flowDirection.magnitude();
    Vector3D dirNorm = (dirMag > 1e-12) ? flowDirection / dirMag : Vector3D(1.0, 0.0, 0.0);
    velocity = dirNorm * soundSpeed * initMach;
    totEnergy = energy + 0.5 * velocity.dot(velocity);
}

FloatType FluidReal::computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const {
    FloatType eta = (std::pow(pressureRatio, (_nominalGamma - 1.0) / _nominalGamma) - 1.0) / (temperatureRatio - 1.0);
    return eta;
}

void FluidReal::setTransportProperties(const Config &config) {
    _viscosityModel = config.getViscosityModel();
    if (_viscosityModel == ViscosityModel::CONSTANT) {
        _muConstant = config.getFluidMuConstant();
    }
    _cpConstant = config.getFluidHeatCapacity();
    _PrConstant = config.getFluidPrandtlNumber();
}

FloatType FluidReal::computeMolecularDynamicViscosity(FloatType temperature) const {
    if (_viscosityModel == ViscosityModel::CONSTANT) {
        return _muConstant;
    }
    // Return nominal or critical viscosity
    return _muConstant;
}

FloatType FluidReal::computeMolecularDynamicViscosity_rho_e(FloatType rho, FloatType e) const {
    if (_viscosityModel == ViscosityModel::CONSTANT) {
        return _muConstant;
    }
    return _table_mu.interpolate(rho, e);
}

FloatType FluidReal::computeThermalConductivity(FloatType dynamicViscosity) const {
    return _cpConstant * dynamicViscosity / _PrConstant;
}

FloatType FluidReal::computeThermalConductivity_rho_e(FloatType rho, FloatType e) const {
    return _table_kappa.interpolate(rho, e);
}

FloatType FluidReal::computeDpDrho_e(FloatType rho, FloatType e) const {
    return _table_dp_drho.interpolate(rho, e);
}

FloatType FluidReal::computeDpDe_rho(FloatType rho, FloatType e) const {
    return _table_dp_de.interpolate(rho, e);
}

FloatType FluidReal::computeInternalEnergy_p_T(FloatType p, FloatType T) const {
    return _table_e_pt.interpolate(p, T);
}

FloatType FluidReal::computeSoundSpeed_p_T(FloatType p, FloatType T) const {
    return _table_a_pt.interpolate(p, T);
}

bool FluidReal::isStateInBounds_rho_e(FloatType rho, FloatType e) const {
    return (rho >= _table_p.xMin && rho <= _table_p.xMax &&
            e   >= _table_p.yMin && e   <= _table_p.yMax);
}
