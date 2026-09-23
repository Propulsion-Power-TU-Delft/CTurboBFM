#pragma once

#include "implicit_solver_base.hpp"
#include "math_utils.hpp"
#include <vector>

/**
 * @brief Lower-Upper Symmetric Gauss-Seidel (LU-SGS) implicit solver.
 *
 * Implements the Yoon & Jameson (1988) / Blazek LU-SGS algorithm for structured grids.
 * Inverts the linearized pseudo-time system (D + L) D^-1 (D + U) deltaU = - R
 * in two sweeps:
 *   1) Forward sweep along diagonal hyperplanes (i + j + k = const)
 *   2) Backward sweep in reverse order
 *
 * Uses scalar diagonal inversion: D_i = (V_i / dt_i + sigma_c + sigma_v) * I
 * which avoids storing large sparse matrices and has O(N) memory complexity.
 */
class ImplicitSolverLUSGS : public ImplicitSolverBase {
public:
    ImplicitSolverLUSGS(
        const Config& config, 
        const Mesh& mesh, 
        const FluidBase& fluid,
        const TurbulenceModelBase* turbModel = nullptr
    );

    ~ImplicitSolverLUSGS() override = default;

    void solveCorrection(
        const FlowSolution& solution,
        const FlowSolution& residuals,
        const Matrix3D<FloatType>& dt,
        FlowSolution& deltaU
    ) override;

private:
    FlowSolution _deltaUStar;
};
