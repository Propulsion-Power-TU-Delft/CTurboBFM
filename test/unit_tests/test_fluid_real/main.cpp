#include <iostream>
#include <cmath>
#include "../../../include/fluid_real.hpp"
#include "gtest/gtest.h"

class FluidRealTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Path relative to execution in test/unit_tests/test_fluid_real
        fluid = std::make_unique<FluidReal>("fluid_table_CO2.lut");
    }

    std::unique_ptr<FluidReal> fluid;
};

TEST_F(FluidRealTest, TestMetadata) {
    EXPECT_EQ(fluid->getFluidName(), "CO2");
    EXPECT_NEAR(fluid->getMolarMass(), 0.04401, 1e-4);
    EXPECT_NEAR(fluid->getCriticalPressure(), 7.3773e6, 1e4);
    EXPECT_NEAR(fluid->getCriticalTemperature(), 304.13, 0.5);
}

TEST_F(FluidRealTest, TestForwardEvaluation_rho_e) {
    // Reference supercritical CO2 state: p=8.0 MPa, T=330 K
    FloatType rho_ref = 198.8735;
    FloatType e_ref = 411627.60;
    FloatType p_expected = 8.0e6;
    FloatType T_expected = 330.0;
    FloatType a_expected = 233.83;

    FloatType p_calc = fluid->computePressure_rho_e(rho_ref, e_ref);
    FloatType T_calc = fluid->computeTemperature_rho_e(rho_ref, e_ref);
    FloatType a_calc = fluid->computeSoundSpeed_rho_e(rho_ref, e_ref);
    FloatType s_calc = fluid->computeEntropy_rho_e(rho_ref, e_ref);

    // Bilinear interpolation tolerance (< 1.5% due to discretization grid)
    EXPECT_NEAR(p_calc, p_expected, p_expected * 0.015);
    EXPECT_NEAR(T_calc, T_expected, T_expected * 0.01);
    EXPECT_NEAR(a_calc, a_expected, a_expected * 0.02);
    EXPECT_GT(s_calc, 1000.0);
}

TEST_F(FluidRealTest, TestInversion_p_rho) {
    FloatType rho_ref = 198.8735;
    FloatType e_ref = 411627.60;
    FloatType p_ref = fluid->computePressure_rho_e(rho_ref, e_ref);

    // Invert (p_ref, rho_ref) -> e
    FloatType e_recovered = fluid->computeStaticEnergy_p_rho(p_ref, rho_ref);
    EXPECT_NEAR(e_recovered, e_ref, e_ref * 0.005);
}

TEST_F(FluidRealTest, TestAuxiliaryQueries_p_T) {
    FloatType p_ref = 8.0e6;
    FloatType T_ref = 330.0;
    FloatType rho_expected = 198.87;

    FloatType rho_calc = fluid->computeDensity_p_T(p_ref, T_ref);
    EXPECT_NEAR(rho_calc, rho_expected, rho_expected * 0.02);
}

TEST_F(FluidRealTest, TestIsentropic_p_s) {
    FloatType p_ref = 8.0e6;
    FloatType T_ref = 330.0;
    FloatType s_ref = fluid->computeEntropy_p_T(p_ref, T_ref);

    FloatType rho_rec = fluid->computeDensity_p_s(p_ref, s_ref);
    FloatType T_rec = fluid->computeTemperature_p_s(p_ref, s_ref);

    EXPECT_NEAR(rho_rec, fluid->computeDensity_p_T(p_ref, T_ref), 5.0);
    EXPECT_NEAR(T_rec, T_ref, 3.0);
}

TEST_F(FluidRealTest, TestStagnationState) {
    FloatType rho = 200.0;
    FloatType e = 410000.0;
    Vector3D vel = {50.0, 0.0, 0.0};
    FloatType et = e + 0.5 * vel.dot(vel);

    FloatType p_static = fluid->computePressure_rho_e(rho, e);
    FloatType p_total = fluid->computeTotalPressure_rho_u_et(rho, vel, et);
    FloatType T_static = fluid->computeTemperature_rho_e(rho, e);
    FloatType T_total = fluid->computeTotalTemperature_rho_u_et(rho, vel, et);

    // Total pressure and temperature must be strictly greater than static
    EXPECT_GT(p_total, p_static);
    EXPECT_GT(T_total, T_static);
}

TEST_F(FluidRealTest, TestInletCharacteristicPhysics) {
    // Reservoir / Inflow total state: pt = 9.0 MPa, Tt = 340 K
    FloatType pt_in = 9.0e6;
    FloatType Tt_in = 340.0;
    FloatType s_in = fluid->computeEntropy_p_T(pt_in, Tt_in);
    FloatType ht_in = fluid->computeEnthalpy_p_s(pt_in, s_in);

    // For any subsonic inlet static pressure pb < pt_in:
    FloatType pb = 8.5e6;
    FloatType hb = fluid->computeEnthalpy_p_s(pb, s_in);
    EXPECT_LT(hb, ht_in);

    // Inflow velocity from enthalpy conservation
    FloatType Vb = std::sqrt(2.0 * (ht_in - hb));
    EXPECT_GT(Vb, 0.0);

    // Density and energy at inlet boundary
    FloatType rhob = fluid->computeDensity_p_s(pb, s_in);
    FloatType eb = fluid->computeInternalEnergy_p_s(pb, s_in);
    EXPECT_GT(rhob, 0.0);
    EXPECT_GT(eb, 0.0);

    // Speed of sound at inlet boundary
    FloatType ab = fluid->computeSoundSpeed_rho_e(rhob, eb);
    FloatType Mb = Vb / ab;
    EXPECT_LT(Mb, 1.0); // Subsonic inlet
}

TEST_F(FluidRealTest, TestOutletIsentropicExpansion) {
    // Interior state
    FloatType p_int = 8.5e6;
    FloatType T_int = 330.0;
    FloatType rho_int = fluid->computeDensity_p_T(p_int, T_int);
    FloatType e_int = fluid->computeInternalEnergy_p_s(p_int, fluid->computeEntropy_p_T(p_int, T_int));
    FloatType s_int = fluid->computeEntropy_rho_e(rho_int, e_int);

    // Imposed lower outlet pressure (expansion)
    FloatType p_exit = 8.0e6;
    FloatType rho_exit = fluid->computeDensity_p_s(p_exit, s_int);
    FloatType e_exit = fluid->computeInternalEnergy_p_s(p_exit, s_int);

    // Density must decrease in isentropic expansion
    EXPECT_LT(rho_exit, rho_int);
    EXPECT_LT(e_exit, e_int);

    // Entropy at exit must match interior entropy
    FloatType s_exit = fluid->computeEntropy_rho_e(rho_exit, e_exit);
    EXPECT_NEAR(s_exit, s_int, 2.0); // within 0.1% of entropy value (~1800 J/kg K)
}

TEST_F(FluidRealTest, TestThermodynamicDerivatives_and_pT) {
    FloatType p_ref = 8.0e6;
    FloatType T_ref = 330.0;
    FloatType rho_ref = fluid->computeDensity_p_T(p_ref, T_ref);
    FloatType e_ref = fluid->computeInternalEnergy_p_T(p_ref, T_ref);

    FloatType dp_drho = fluid->computeDpDrho_e(rho_ref, e_ref);
    FloatType dp_de = fluid->computeDpDe_rho(rho_ref, e_ref);

    // Both derivatives must be strictly positive for single-phase stable gas
    EXPECT_GT(dp_drho, 0.0);
    EXPECT_GT(dp_de, 0.0);

    // Speed of sound consistency: a^2 approx dp/drho + (p/rho^2) * dp/de
    FloatType a2_deriv = dp_drho + (p_ref / (rho_ref * rho_ref)) * dp_de;
    FloatType a_table = fluid->computeSoundSpeed_rho_e(rho_ref, e_ref);
    EXPECT_GT(a2_deriv, 0.0);
    EXPECT_NEAR(std::sqrt(a2_deriv), a_table, a_table * 0.02);

    // (p, T) sound speed vs (rho, e) sound speed
    FloatType a_pt = fluid->computeSoundSpeed_p_T(p_ref, T_ref);
    EXPECT_NEAR(a_pt, a_table, 1.0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
