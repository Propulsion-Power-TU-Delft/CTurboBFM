#include <iostream>
#include "../../../include/fluid_ideal.hpp"
#include "gtest/gtest.h"


TEST(FluidIdealTest, TestStaticEnergy_p_rho) {
    FloatType gamma {2.0}, R {1.0};
    FluidIdeal fluid(gamma, R);

    std::vector<FloatType> testDensity {1.0};
    std::vector<FloatType> testPressure {1.0};
    std::vector<FloatType> expectedStaticEnergy {1.0};

    for (int i=0; i<testDensity.size(); i++){
        ASSERT_DOUBLE_EQ(fluid.computeStaticEnergy_p_rho(testPressure[i], testDensity[i]), expectedStaticEnergy[i]);
    }
    
}

TEST(FluidIdealTest, TestSoundSpeed_p_rho) {
    FloatType gamma {2.0}, R {1.0};
    FluidIdeal fluid(gamma, R);
    
    FloatType rho {1.0};
    FloatType p {1.0};

    FloatType expectedValue {std::sqrt(2.0)};

    std::vector<FloatType> computedStaticEnergy;
    ASSERT_DOUBLE_EQ(fluid.computeSoundSpeed_p_rho(p, rho), expectedValue);
}

TEST(FluidIdealTest, TestThermodynamics_rho_e) {
    FloatType gamma {1.4}, R {287.05};
    FluidIdeal fluid(gamma, R);

    FloatType T_ref = 300.0;
    FloatType p_ref = 101325.0;
    FloatType rho_ref = fluid.computeDensity_p_T(p_ref, T_ref);
    FloatType cv = R / (gamma - 1.0);
    FloatType e_ref = cv * T_ref;

    // Temperature from (rho, e)
    EXPECT_NEAR(fluid.computeTemperature_rho_e(rho_ref, e_ref), T_ref, 1e-10);

    // Pressure from (rho, e)
    EXPECT_NEAR(fluid.computePressure_rho_e(rho_ref, e_ref), p_ref, 1e-6);

    // Sound speed from (rho, e) vs (p, rho)
    FloatType a1 = fluid.computeSoundSpeed_rho_e(rho_ref, e_ref);
    FloatType a2 = fluid.computeSoundSpeed_p_rho(p_ref, rho_ref);
    EXPECT_DOUBLE_EQ(a1, a2);
    EXPECT_NEAR(a1, std::sqrt(gamma * R * T_ref), 1e-10);

    // Isentropic exponent and fundamental derivative
    EXPECT_NEAR(fluid.computeIsentropicExponent_rho_e(rho_ref, e_ref), gamma, 1e-10);
    EXPECT_NEAR(fluid.computeFundamentalDerivative_rho_e(rho_ref, e_ref), 0.5 * (gamma + 1.0), 1e-10);
}

TEST(FluidIdealTest, TestIsentropicState_p_s) {
    FloatType gamma {1.4}, R {287.05};
    FluidIdeal fluid(gamma, R);

    FloatType T_orig = 350.0;
    FloatType p_orig = 200000.0;
    FloatType s_orig = fluid.computeEntropy_p_T(p_orig, T_orig);

    // Round-trip (p_orig, s_orig) -> T, rho, e
    FloatType T_recovered = fluid.computeTemperature_p_s(p_orig, s_orig);
    FloatType rho_recovered = fluid.computeDensity_p_s(p_orig, s_orig);
    FloatType e_recovered = fluid.computeInternalEnergy_p_s(p_orig, s_orig);

    EXPECT_NEAR(T_recovered, T_orig, 1e-10);
    EXPECT_NEAR(rho_recovered, fluid.computeDensity_p_T(p_orig, T_orig), 1e-10);
    EXPECT_NEAR(e_recovered, fluid.computeStaticEnergy_p_rho(p_orig, rho_recovered), 1e-10);

    // Isentropic expansion to p_expanded: s remains constant
    FloatType p_expanded = 100000.0;
    FloatType T_expanded = fluid.computeTemperature_p_s(p_expanded, s_orig);
    FloatType expectedT_expanded = T_orig * std::pow(p_expanded / p_orig, (gamma - 1.0) / gamma);
    EXPECT_NEAR(T_expanded, expectedT_expanded, 1e-8);
}



TEST(FluidIdealTest, TestConstantViscosity) {
    FloatType gamma {1.4}, R {287.05};
    FluidIdeal fluid(gamma, R);

    FloatType muConst {1.5e-5};
    FloatType cp {1006.0};
    FloatType Pr {0.71};
    fluid.setConstantViscosity(muConst, cp, Pr);

    EXPECT_EQ(fluid.getViscosityModel(), ViscosityModel::CONSTANT);
    EXPECT_DOUBLE_EQ(fluid.computeMolecularDynamicViscosity(300.0), muConst);
    EXPECT_DOUBLE_EQ(fluid.computeMolecularDynamicViscosity(500.0), muConst);
    EXPECT_DOUBLE_EQ(fluid.computeMolecularDynamicViscosity(100.0), muConst);
    EXPECT_DOUBLE_EQ(fluid.computeThermalConductivity(muConst), cp * muConst / Pr);
}

TEST(FluidIdealTest, TestSutherlandViscosity) {
    FloatType gamma {1.4}, R {287.05};
    FluidIdeal fluid(gamma, R);

    FloatType muRef {1.716e-5};
    FloatType TRef {273.15};
    FloatType S {110.4};
    FloatType cp {1006.0};
    FloatType Pr {0.71};
    fluid.setSutherlandViscosity(muRef, TRef, S, cp, Pr);

    EXPECT_EQ(fluid.getViscosityModel(), ViscosityModel::SUTHERLAND);
    // At T = TRef, mu should equal muRef
    EXPECT_NEAR(fluid.computeMolecularDynamicViscosity(TRef), muRef, 1e-12);
    // At T = 300 K
    FloatType expectedMu = muRef * std::pow(300.0 / TRef, 1.5) * ((TRef + S) / (300.0 + S));
    EXPECT_NEAR(fluid.computeMolecularDynamicViscosity(300.0), expectedMu, 1e-12);
}

TEST(FluidIdealTest, TestThermodynamicDerivatives_and_pT) {
    FloatType gamma {1.4}, R {287.05};
    FluidIdeal fluid(gamma, R);

    FloatType T = 300.0;
    FloatType p = 101325.0;
    FloatType rho = fluid.computeDensity_p_T(p, T);
    FloatType cv = R / (gamma - 1.0);
    FloatType e = cv * T;

    // dp/drho|e = (gamma - 1) * e
    EXPECT_NEAR(fluid.computeDpDrho_e(rho, e), (gamma - 1.0) * e, 1e-10);
    // dp/de|rho = (gamma - 1) * rho
    EXPECT_NEAR(fluid.computeDpDe_rho(rho, e), (gamma - 1.0) * rho, 1e-10);

    // Consistency: a^2 = dp/drho|e + (p / rho^2) * dp/de|rho
    FloatType a2_deriv = fluid.computeDpDrho_e(rho, e) + (p / (rho * rho)) * fluid.computeDpDe_rho(rho, e);
    FloatType a_expected = std::sqrt(gamma * R * T);
    EXPECT_NEAR(std::sqrt(a2_deriv), a_expected, 1e-8);

    // (p, T) queries
    EXPECT_NEAR(fluid.computeInternalEnergy_p_T(p, T), e, 1e-10);
    EXPECT_NEAR(fluid.computeSoundSpeed_p_T(p, T), a_expected, 1e-10);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}