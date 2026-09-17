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

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}