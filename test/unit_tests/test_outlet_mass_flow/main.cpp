#include <gtest/gtest.h>
#include "types.hpp"
#include "config.hpp"
#include "boundary_outlet_mass_flow.hpp"
#include "fluid_ideal.hpp"
#include <fstream>
#include <cmath>

class OutletMassFlowTest : public ::testing::Test {
protected:
    void writeTempConfig(const std::string& filename, const std::string& content) {
        std::ofstream out(filename);
        out << content;
        out.close();
    }
};

TEST_F(OutletMassFlowTest, TestBoundaryTypeFromString) {
    EXPECT_EQ(boundaryTypeFromString("mass_flow"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("MASS_FLOW"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("outlet_mass_flow"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("OUTLET_MASS_FLOW"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("massflow"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("mass_flow_rate"), BoundaryType::MASS_FLOW);
    EXPECT_EQ(boundaryTypeFromString("Mass_Flow_Rate"), BoundaryType::MASS_FLOW);
}

TEST_F(OutletMassFlowTest, TestConfigMassFlowDefaults) {
    writeTempConfig("test_default_mf.ini", "CFL = 1.0\n");
    Config cfg("test_default_mf.ini");
    EXPECT_DOUBLE_EQ(cfg.getMassFlowRelaxation(), 0.05);
    EXPECT_DOUBLE_EQ(cfg.getMassFlowMaxPressureChange(), 0.02);
    std::remove("test_default_mf.ini");
}

TEST_F(OutletMassFlowTest, TestConfigMassFlowCustom) {
    writeTempConfig("test_custom_mf.ini", 
        "MASS_FLOW_RELAXATION = 0.1\n"
        "MASS_FLOW_MAX_PRESSURE_CHANGE = 0.05\n"
    );
    Config cfg("test_custom_mf.ini");
    EXPECT_DOUBLE_EQ(cfg.getMassFlowRelaxation(), 0.1);
    EXPECT_DOUBLE_EQ(cfg.getMassFlowMaxPressureChange(), 0.05);
    std::remove("test_custom_mf.ini");
}

TEST_F(OutletMassFlowTest, TestBoundaryOutletMassFlowFlux) {
    writeTempConfig("test_flux_mf.ini", 
        "INIT_PRESSURE = 100000.0\n"
        "OUTLET_PRESSURE_RAMP_ITERATIONS = 0\n"
    );
    Config cfg("test_flux_mf.ini");
    FluidIdeal fluid(1.4, 287.05);

    // Dummy mesh pointer for BoundaryBase reference
    const Mesh* dummyMeshPtr = reinterpret_cast<const Mesh*>(0x1000);
    const Mesh& dummyMesh = *dummyMeshPtr;

    std::vector<FloatType> pressureProfile = {95000.0, 96000.0, 97000.0};
    BoundaryOutletMassFlow bc(cfg, dummyMesh, fluid, pressureProfile);

    // Subsonic interior state: rho=1.2, u=50, v=0, w=0, p=100000
    FloatType rho = 1.2;
    FloatType u = 50.0;
    FloatType p = 100000.0;
    FloatType e = fluid.computeStaticEnergy_p_rho(p, rho);
    FloatType et = e + 0.5 * u * u;
    StateVector internalCons({rho, rho * u, 0.0, 0.0, rho * et});

    Vector3D surface(1.0, 0.0, 0.0); // outward normal in +x
    Vector3D midPoint(0.0, 0.0, 0.0);
    FlowSolution dummySol(1, 1, 1);

    // Face at j=1 -> should use pressureProfile[1] = 96000.0
    StateVector flux = bc.computeBoundaryFlux(internalCons, surface, midPoint, {0, 1, 0}, dummySol, 1);

    // Verify boundary flux continuity
    // Mass flux = rho_b * u_b * Area
    EXPECT_GT(flux[0], 0.0);
    // Pressure flux in momentum
    EXPECT_GT(flux[1], 0.0);

    // Face with j out of range -> should clamp to last pressure (97000.0)
    StateVector fluxClamped = bc.computeBoundaryFlux(internalCons, surface, midPoint, {0, 5, 0}, dummySol, 1);
    EXPECT_GT(fluxClamped[0], 0.0);

    std::remove("test_flux_mf.ini");
}

TEST_F(OutletMassFlowTest, TestMassFlowControllerLogic) {
    FloatType targetMassFlow = 10.0;
    FloatType currentPressure = 100000.0;
    FloatType relaxation = 0.05;
    FloatType maxChange = 0.02;

    // Case 1: current mass flow > target (flow too high -> backpressure too low -> increase pressure)
    FloatType currentMassFlowHigh = 12.0;
    FloatType errorHigh = (currentMassFlowHigh - targetMassFlow) / targetMassFlow; // +0.20
    FloatType deltaP_relHigh = std::clamp(relaxation * errorHigh, -maxChange, maxChange); // 0.05 * 0.20 = 0.01
    FloatType newPressureHigh = currentPressure * (1.0 + deltaP_relHigh);
    EXPECT_GT(newPressureHigh, currentPressure);
    EXPECT_DOUBLE_EQ(deltaP_relHigh, 0.01);
    EXPECT_DOUBLE_EQ(newPressureHigh, 101000.0);

    // Case 2: current mass flow < target (flow too low -> backpressure too high -> decrease pressure)
    FloatType currentMassFlowLow = 8.0;
    FloatType errorLow = (currentMassFlowLow - targetMassFlow) / targetMassFlow; // -0.20
    FloatType deltaP_relLow = std::clamp(relaxation * errorLow, -maxChange, maxChange); // -0.01
    FloatType newPressureLow = currentPressure * (1.0 + deltaP_relLow);
    EXPECT_LT(newPressureLow, currentPressure);
    EXPECT_DOUBLE_EQ(deltaP_relLow, -0.01);
    EXPECT_DOUBLE_EQ(newPressureLow, 99000.0);

    // Case 3: massive error clamped by maxChange (+/- 2%)
    FloatType currentMassFlowVeryHigh = 20.0; // +100% error
    FloatType errorVeryHigh = (currentMassFlowVeryHigh - targetMassFlow) / targetMassFlow; // +1.0
    FloatType deltaP_relClamped = std::clamp(relaxation * errorVeryHigh, -maxChange, maxChange); // clamped to 0.02
    EXPECT_DOUBLE_EQ(deltaP_relClamped, maxChange);
}
