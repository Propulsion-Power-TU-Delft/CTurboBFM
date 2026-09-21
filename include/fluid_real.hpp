#pragma once
#include "types.hpp"
#include "fluid_base.hpp"
#include "config.hpp"
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

struct Table2D {
    size_t nX {0};
    size_t nY {0};
    FloatType xMin {0.0};
    FloatType xMax {0.0};
    FloatType yMin {0.0};
    FloatType yMax {0.0};
    FloatType dx {1.0};
    FloatType dy {1.0};
    std::vector<FloatType> data;

    void init(size_t nx, size_t ny, FloatType xmin, FloatType xmax, FloatType ymin, FloatType ymax) {
        nX = nx;
        nY = ny;
        xMin = xmin;
        xMax = xmax;
        yMin = ymin;
        yMax = ymax;
        dx = (nX > 1) ? (xMax - xMin) / static_cast<FloatType>(nX - 1) : 1.0;
        dy = (nY > 1) ? (yMax - yMin) / static_cast<FloatType>(nY - 1) : 1.0;
        data.assign(nX * nY, 0.0);
    }

    inline FloatType interpolate(FloatType x, FloatType y) const {
        FloatType xc = std::clamp(x, xMin, xMax);
        FloatType yc = std::clamp(y, yMin, yMax);

        FloatType rx = (xc - xMin) / dx;
        FloatType ry = (yc - yMin) / dy;

        size_t i = static_cast<size_t>(rx);
        size_t j = static_cast<size_t>(ry);

        if (i >= nX - 1) i = (nX >= 2) ? nX - 2 : 0;
        if (j >= nY - 1) j = (nY >= 2) ? nY - 2 : 0;

        FloatType wx = rx - static_cast<FloatType>(i);
        FloatType wy = ry - static_cast<FloatType>(j);

        FloatType f00 = data[i * nY + j];
        FloatType f10 = data[(i + 1) * nY + j];
        FloatType f01 = data[i * nY + (j + 1)];
        FloatType f11 = data[(i + 1) * nY + (j + 1)];

        return (1.0 - wx) * (1.0 - wy) * f00 +
               wx * (1.0 - wy) * f10 +
               (1.0 - wx) * wy * f01 +
               wx * wy * f11;
    }
};

class FluidReal : public FluidBase {
public:
    explicit FluidReal(const std::string& tablePath);

    // Thermodynamic evaluations from (rho, e)
    FloatType computePressure_rho_e(FloatType rho, FloatType e) const override;
    FloatType computeTemperature_rho_e(FloatType rho, FloatType e) const override;
    FloatType computeSoundSpeed_rho_e(FloatType rho, FloatType e) const override;
    FloatType computeEntropy_rho_e(FloatType rho, FloatType e) const override;

    // Thermodynamic evaluations from (p, rho)
    FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const override;
    FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const override;
    FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const override;

    // Thermodynamic evaluations from (p, T) and (rho, T)
    FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const override;
    FloatType computeDensity_p_T(FloatType p, FloatType T) const override;
    FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const override;

    // Thermodynamic evaluations from (p, s)
    FloatType computeDensity_p_s(FloatType p, FloatType s) const override;
    FloatType computeTemperature_p_s(FloatType p, FloatType s) const override;
    FloatType computeInternalEnergy_p_s(FloatType p, FloatType s) const override;

    // Stagnation and isentropic relations
    FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const override;
    FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const override;
    FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const override;
    FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const override;

    void computeInitFields(
        FloatType initMach, 
        FloatType initTemperature, 
        FloatType initPressure, 
        Vector3D flowDirection, 
        FloatType &density, 
        Vector3D &velocity, 
        FloatType &totEnergy) override;

    FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const override;

    FloatType getGamma() const override { return _nominalGamma; }
    FloatType getRconstant() const override { return _gasConstant; }

    void setTransportProperties(const Config &config) override;
    ViscosityModel getViscosityModel() const override { return _viscosityModel; }
    FloatType computeMolecularDynamicViscosity(FloatType temperature) const override;
    FloatType computeMolecularDynamicViscosity_rho_e(FloatType rho, FloatType e) const;
    FloatType computeThermalConductivity(FloatType dynamicViscosity) const override;
    FloatType computeThermalConductivity_rho_e(FloatType rho, FloatType e) const;

    // Thermodynamic derivatives
    FloatType computeDpDrho_e(FloatType rho, FloatType e) const;
    FloatType computeDpDe_rho(FloatType rho, FloatType e) const;

    const std::string& getFluidName() const { return _fluidName; }
    FloatType getMolarMass() const { return _molarMass; }
    FloatType getCriticalPressure() const { return _pCrit; }
    FloatType getCriticalTemperature() const { return _TCrit; }
    FloatType getCriticalDensity() const { return _rhoCrit; }

private:
    void loadTable(const std::string& tablePath);

    std::string _fluidName {"Unknown"};
    FloatType _molarMass {0.0};
    FloatType _gasConstant {287.0};
    FloatType _pCrit {0.0};
    FloatType _TCrit {0.0};
    FloatType _rhoCrit {0.0};
    FloatType _nominalGamma {1.3};

    // Primary (rho, e) tables
    Table2D _table_p;       // pressure [Pa]
    Table2D _table_T;       // temperature [K]
    Table2D _table_a;       // sound speed [m/s]
    Table2D _table_s;       // entropy [J/(kg K)]
    Table2D _table_mu;      // dynamic viscosity [Pa s]
    Table2D _table_kappa;   // thermal conductivity [W/(m K)]
    Table2D _table_dp_drho; // dp/drho|e [Pa m^3 / kg]
    Table2D _table_dp_de;   // dp/de|rho [Pa kg / J]

    // Auxiliary (p, T) tables
    Table2D _table_rho_pt;  // density from (p, T)
    Table2D _table_e_pt;    // energy from (p, T)
    Table2D _table_a_pt;    // sound speed from (p, T)
    Table2D _table_s_pt;    // entropy from (p, T)

    // Auxiliary (p, s) tables
    Table2D _table_rho_ps;  // density from (p, s)
    Table2D _table_T_ps;    // temperature from (p, s)
    Table2D _table_e_ps;    // energy from (p, s)

    ViscosityModel _viscosityModel {ViscosityModel::CONSTANT};
    FloatType _muConstant {2.0e-5};
    FloatType _cpConstant {1000.0};
    FloatType _PrConstant {0.72};
};
