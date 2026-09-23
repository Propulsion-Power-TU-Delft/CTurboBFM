#include <gtest/gtest.h>
#include "config.hpp"
#include <fstream>

class ImplicitConfigTest : public ::testing::Test {
protected:
    void writeTempConfig(const std::string& filename, const std::string& content) {
        std::ofstream out(filename);
        out << content;
        out.close();
    }
};

TEST_F(ImplicitConfigTest, TestLegacyNumericValues) {
    writeTempConfig("test_legacy_0.ini", "TIME_INTEGRATION_TYPE = 0\n");
    Config c0("test_legacy_0.ini");
    EXPECT_EQ(c0.getTimeIntegration(), TimeIntegration::RUNGE_KUTTA_4);
    EXPECT_FALSE(c0.isTimeIntegrationImplicit());

    writeTempConfig("test_legacy_1.ini", "TIME_INTEGRATION_TYPE = 1\n");
    Config c1("test_legacy_1.ini");
    EXPECT_EQ(c1.getTimeIntegration(), TimeIntegration::RUNGE_KUTTA_3);
    EXPECT_FALSE(c1.isTimeIntegrationImplicit());

    writeTempConfig("test_legacy_2.ini", "TIME_INTEGRATION_TYPE = 2\n");
    Config c2("test_legacy_2.ini");
    EXPECT_EQ(c2.getTimeIntegration(), TimeIntegration::IMPLICIT_LU_SGS);
    EXPECT_TRUE(c2.isTimeIntegrationImplicit());

    writeTempConfig("test_legacy_3.ini", "TIME_INTEGRATION_TYPE = 3\n");
    Config c3("test_legacy_3.ini");
    EXPECT_EQ(c3.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_TRUE(c3.isTimeIntegrationImplicit());

    std::remove("test_legacy_0.ini");
    std::remove("test_legacy_1.ini");
    std::remove("test_legacy_2.ini");
    std::remove("test_legacy_3.ini");
}

TEST_F(ImplicitConfigTest, TestStringValues) {
    writeTempConfig("test_str_rk4.ini", "TIME_INTEGRATION_TYPE = explicit_rk4\n");
    Config c_rk4("test_str_rk4.ini");
    EXPECT_EQ(c_rk4.getTimeIntegration(), TimeIntegration::RUNGE_KUTTA_4);

    writeTempConfig("test_str_rk3.ini", "TIME_INTEGRATION_TYPE = explicit_rk3\n");
    Config c_rk3("test_str_rk3.ini");
    EXPECT_EQ(c_rk3.getTimeIntegration(), TimeIntegration::RUNGE_KUTTA_3);

    writeTempConfig("test_str_lusgs.ini", "TIME_INTEGRATION_TYPE = implicit_lu_sgs\n");
    Config c_lusgs("test_str_lusgs.ini");
    EXPECT_EQ(c_lusgs.getTimeIntegration(), TimeIntegration::IMPLICIT_LU_SGS);
    EXPECT_TRUE(c_lusgs.isTimeIntegrationImplicit());

    writeTempConfig("test_str_krylov.ini", "TIME_INTEGRATION_TYPE = implicit_krylov\n");
    Config c_krylov("test_str_krylov.ini");
    EXPECT_EQ(c_krylov.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_TRUE(c_krylov.isTimeIntegrationImplicit());

    // Test whitespace tolerance (e.g. "explicit_ rk4")
    writeTempConfig("test_str_space.ini", "TIME_INTEGRATION_TYPE = explicit_ rk4\n");
    Config c_space("test_str_space.ini");
    EXPECT_EQ(c_space.getTimeIntegration(), TimeIntegration::RUNGE_KUTTA_4);

    // Test case-insensitivity
    writeTempConfig("test_str_caps.ini", "TIME_INTEGRATION = IMPLICIT_LU_SGS\n");
    Config c_caps("test_str_caps.ini");
    EXPECT_EQ(c_caps.getTimeIntegration(), TimeIntegration::IMPLICIT_LU_SGS);

    std::remove("test_str_rk4.ini");
    std::remove("test_str_rk3.ini");
    std::remove("test_str_lusgs.ini");
    std::remove("test_str_krylov.ini");
    std::remove("test_str_space.ini");
    std::remove("test_str_caps.ini");
}

TEST_F(ImplicitConfigTest, TestImplicitParameters) {
    std::string configContent = 
        "TIME_INTEGRATION_TYPE = implicit_lu_sgs\n"
        "CFL = 5.0\n"
        "CFL_START = 1.5\n"
        "CFL_MAX = 50.0\n"
        "CFL_RAMP_ITERATIONS = 100\n"
        "IMPLICIT_UNDER_RELAXATION = 0.8\n"
        "LINEAR_SOLVER_TOL = 1e-4\n"
        "LINEAR_SOLVER_MAX_ITER = 30\n";

    writeTempConfig("test_params.ini", configContent);
    Config c("test_params.ini");

    EXPECT_DOUBLE_EQ(c.getCFL(), 5.0);
    EXPECT_DOUBLE_EQ(c.getCFLStart(), 1.5);
    EXPECT_DOUBLE_EQ(c.getCFLMax(), 50.0);
    EXPECT_EQ(c.getCFLRampIterations(), 100);
    EXPECT_DOUBLE_EQ(c.getImplicitUnderRelaxation(), 0.8);
    EXPECT_DOUBLE_EQ(c.getLinearSolverTol(), 1e-4);
    EXPECT_EQ(c.getLinearSolverMaxIter(), 30);

    std::remove("test_params.ini");
}

TEST_F(ImplicitConfigTest, TestImplicitCFLFallbackWithoutCFLKey) {
    std::string configContent = 
        "CFL_START = 5.0\n"
        "CFL_MAX = 50.0\n"
        "CFL_RAMP_ITERATIONS = 200\n"
        "IMPLICIT_UNDER_RELAXATION = 1.0\n"
        "LINEAR_SOLVER_TOL = 1e-2\n"
        "LINEAR_SOLVER_MAX_ITER = 20\n"
        "TIME_INTEGRATION_TYPE = implicit_lu_sgs\n";

    writeTempConfig("test_no_cfl.ini", configContent);
    Config c("test_no_cfl.ini");

    // Must not throw an exception when querying getCFL()
    EXPECT_NO_THROW({
        EXPECT_DOUBLE_EQ(c.getCFL(), 5.0);
    });
    EXPECT_DOUBLE_EQ(c.getCFLStart(), 5.0);
    EXPECT_DOUBLE_EQ(c.getCFLMax(), 50.0);
    EXPECT_EQ(c.getCFLRampIterations(), 200);
    EXPECT_DOUBLE_EQ(c.getImplicitUnderRelaxation(), 1.0);
    EXPECT_DOUBLE_EQ(c.getLinearSolverTol(), 1e-2);
    EXPECT_EQ(c.getLinearSolverMaxIter(), 20);
    EXPECT_EQ(c.getTimeIntegration(), TimeIntegration::IMPLICIT_LU_SGS);

    std::remove("test_no_cfl.ini");
}

TEST_F(ImplicitConfigTest, TestGMRESConfig) {
    std::string configContent = 
        "TIME_INTEGRATION_TYPE = implicit_krylov\n"
        "LINEAR_SOLVER_TYPE = gmres\n"
        "KRYLOV_RESTART = 15\n"
        "LINEAR_SOLVER_PRECONDITIONER = diagonal\n";

    writeTempConfig("test_gmres.ini", configContent);
    Config c("test_gmres.ini");

    EXPECT_EQ(c.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_EQ(c.getLinearSolverType(), LinearSolverType::GMRES);
    EXPECT_EQ(c.getLinearPreconditionerType(), LinearPreconditionerType::DIAGONAL);
    EXPECT_EQ(c.getKrylovRestart(), 15);
    std::remove("test_gmres.ini");

    // Test alias TIME_INTEGRATION_TYPE = implicit_gmres
    writeTempConfig("test_gmres_alias.ini", "TIME_INTEGRATION_TYPE = implicit_gmres\n");
    Config c_alias("test_gmres_alias.ini");
    EXPECT_EQ(c_alias.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_EQ(c_alias.getLinearSolverType(), LinearSolverType::GMRES);
    std::remove("test_gmres_alias.ini");
}

TEST_F(ImplicitConfigTest, TestFGMRESConfig) {
    std::string configContent = 
        "TIME_INTEGRATION_TYPE = implicit_krylov\n"
        "LINEAR_SOLVER_TYPE = fgmres\n"
        "KRYLOV_RESTART = 25\n"
        "LINEAR_SOLVER_PRECONDITIONER = lu_sgs\n";

    writeTempConfig("test_fgmres.ini", configContent);
    Config c("test_fgmres.ini");

    EXPECT_EQ(c.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_EQ(c.getLinearSolverType(), LinearSolverType::FGMRES);
    EXPECT_EQ(c.getLinearPreconditionerType(), LinearPreconditionerType::LU_SGS);
    EXPECT_EQ(c.getKrylovRestart(), 25);
    std::remove("test_fgmres.ini");

    // Test alias TIME_INTEGRATION_TYPE = implicit_fgmres with default preconditioner
    writeTempConfig("test_fgmres_alias.ini", "TIME_INTEGRATION_TYPE = implicit_fgmres\n");
    Config c_alias("test_fgmres_alias.ini");
    EXPECT_EQ(c_alias.getTimeIntegration(), TimeIntegration::IMPLICIT_KRYLOV);
    EXPECT_EQ(c_alias.getLinearSolverType(), LinearSolverType::FGMRES);
    // For FGMRES, default preconditioner should be LU_SGS
    EXPECT_EQ(c_alias.getLinearPreconditionerType(), LinearPreconditionerType::LU_SGS);
    std::remove("test_fgmres_alias.ini");
}
