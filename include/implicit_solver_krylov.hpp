#pragma once

#include "implicit_solver_base.hpp"
#include "math_utils.hpp"
#include <Eigen/Sparse>
#include <vector>

/**
 * @brief Krylov subspace implicit linear solver using Eigen.
 *
 * Supports:
 * - BiCGSTAB (with diagonal preconditioner)
 * - GMRES(m) (via Eigen unsupported GMRES)
 * - FGMRES(m) (Flexible GMRES with matrix-free LU-SGS right preconditioning)
 */
class ImplicitSolverKrylov : public ImplicitSolverBase {
public:
    ImplicitSolverKrylov(
        const Config& config, 
        const Mesh& mesh, 
        const FluidBase& fluid,
        const TurbulenceModelBase* turbModel = nullptr
    );

    ~ImplicitSolverKrylov() override = default;

    void solveCorrection(
        const FlowSolution& solution,
        const FlowSolution& residuals,
        const Matrix3D<FloatType>& dt,
        FlowSolution& deltaU
    ) override;

private:
    LinearSolverType _solverType{LinearSolverType::BICGSTAB};
    LinearPreconditionerType _precondType{LinearPreconditionerType::DIAGONAL};
    FloatType _tolerance{1e-2};
    size_t _maxIterations{20};
    size_t _restart{20};

    void applyLUSGSPreconditioner(
        const Eigen::VectorXd& v,
        Eigen::VectorXd& z,
        const Matrix3D<StateVector>& prim,
        const Matrix3D<FloatType>& soundSpeed,
        const Matrix3D<FloatType>& diagD,
        FloatType gamma
    ) const;

    Eigen::VectorXd solveFGMRES(
        const Eigen::SparseMatrix<FloatType>& A,
        const Eigen::VectorXd& b,
        const Matrix3D<StateVector>& prim,
        const Matrix3D<FloatType>& soundSpeed,
        const Matrix3D<FloatType>& diagD,
        FloatType gamma
    );

    mutable FlowSolution _lusgsDeltaUStar;
    mutable FlowSolution _lusgsDeltaU;
};
