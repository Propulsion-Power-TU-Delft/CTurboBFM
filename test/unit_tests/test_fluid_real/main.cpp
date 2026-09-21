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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
