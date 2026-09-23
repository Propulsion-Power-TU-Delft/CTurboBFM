#include "implicit_solver_krylov.hpp"
#include "turbulence_model_base.hpp"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <cmath>

ImplicitSolverKrylov::ImplicitSolverKrylov(
    const Config& config, 
    const Mesh& mesh, 
    const FluidBase& fluid,
    const TurbulenceModelBase* turbModel
) : ImplicitSolverBase(config, mesh, fluid, turbModel),
    _solverType(config.getLinearSolverType()),
    _precondType(config.getLinearPreconditionerType()),
    _tolerance(config.getLinearSolverTol()),
    _maxIterations(config.getLinearSolverMaxIter()),
    _restart(config.getKrylovRestart()),
    _lusgsDeltaUStar(mesh.getNumberPointsI(), mesh.getNumberPointsJ(), mesh.getNumberPointsK()),
    _lusgsDeltaU(mesh.getNumberPointsI(), mesh.getNumberPointsJ(), mesh.getNumberPointsK()) {
}

void ImplicitSolverKrylov::solveCorrection(
    const FlowSolution& solution,
    const FlowSolution& residuals,
    const Matrix3D<FloatType>& dt,
    FlowSolution& deltaU
) {
    const FloatType gamma = _fluid.getGamma();
    const bool isViscous = _config.isViscosityActive();
    const bool is3D = (_topology == Topology::THREE_DIMENSIONAL);
    const bool is1D = (_topology == Topology::ONE_DIMENSIONAL);

    const auto& surfacesI = _mesh.getSurfacesI();
    const auto& surfacesJ = _mesh.getSurfacesJ();
    const auto& surfacesK = _mesh.getSurfacesK();
    const auto& volumes = _mesh.getVolumes();

    const size_t totalCells = _ni * _nj * _nk;
    const size_t numDofs = totalCells * 5;

    // Helper lambda for cell flat index
    auto cellIndex = [&](size_t i, size_t j, size_t k) -> size_t {
        return (i * _nj + j) * _nk + k;
    };

    // Precompute primitive variables and sound speeds
    Matrix3D<StateVector> prim(_ni, _nj, _nk);
    Matrix3D<FloatType> soundSpeed(_ni, _nj, _nk);
    Matrix3D<FloatType> diagD(_ni, _nj, _nk);

    #pragma omp parallel for collapse(2) schedule(static) if(_ni * _nj >= 64)
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                StateVector cons = solution.at(i, j, k);
                StateVector p = getPrimitiveVariablesFromConservative(cons);
                prim(i, j, k) = p;

                Vector3D vel(p[1], p[2], p[3]);
                FloatType a = _fluid.computeSoundSpeed_rho_u_et(p[0], vel, p[4]);
                soundSpeed(i, j, k) = a;

                FloatType lamW = computeFaceSpectralRadius(p, surfacesI(i, j, k), a);
                FloatType lamE = computeFaceSpectralRadius(p, surfacesI(i + 1, j, k), a);

                FloatType lamS = 0.0, lamN = 0.0;
                if (!is1D) {
                    lamS = computeFaceSpectralRadius(p, surfacesJ(i, j, k), a);
                    lamN = computeFaceSpectralRadius(p, surfacesJ(i, j + 1, k), a);
                }

                FloatType lamB = 0.0, lamT = 0.0;
                if (is3D) {
                    lamB = computeFaceSpectralRadius(p, surfacesK(i, j, k), a);
                    lamT = computeFaceSpectralRadius(p, surfacesK(i, j, k + 1), a);
                }

                FloatType sigmaC = 0.5 * (lamW + lamE + lamS + lamN + lamB + lamT);

                FloatType sigmaV = 0.0;
                if (isViscous) {
                    FloatType T = _fluid.computeTemperature_rho_u_et(p[0], vel, p[4]);
                    FloatType muLam = _fluid.computeMolecularDynamicViscosity(T);
                    FloatType muEddy = 0.0;
                    if (_turbModel != nullptr) {
                        muEddy = _turbModel->getEddyViscosity(p[0], i, j, k);
                    }
                    FloatType muTot = muLam + muEddy;
                    FloatType vol = volumes(i, j, k);
                    FloatType sTotSq = surfacesI(i, j, k).magnitudeSquared() + surfacesI(i + 1, j, k).magnitudeSquared();
                    if (!is1D) {
                        sTotSq += surfacesJ(i, j, k).magnitudeSquared() + surfacesJ(i, j + 1, k).magnitudeSquared();
                    }
                    if (is3D) {
                        sTotSq += surfacesK(i, j, k).magnitudeSquared() + surfacesK(i, j, k + 1).magnitudeSquared();
                    }
                    sigmaV = (muTot / (p[0] * vol)) * sTotSq;
                }

                FloatType vol = volumes(i, j, k);
                FloatType dtLocal = dt(i, j, k);
                diagD(i, j, k) = (vol / dtLocal) + sigmaC + sigmaV;
            }
        }
    }

    // Assemble sparse matrix triplets and RHS vector
    std::vector<Eigen::Triplet<FloatType>> triplets;
    // Reserve estimation: ~7 blocks of 5x5 = 175 non-zeros per cell
    triplets.reserve(totalCells * 175);

    Eigen::VectorXd rhs(numDofs);

    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                const size_t cIdx = cellIndex(i, j, k);
                const size_t baseRow = cIdx * 5;

                // RHS = - residuals
                StateVector res = residuals.at(i, j, k);
                for (size_t v = 0; v < 5; ++v) {
                    rhs(baseRow + v) = -res[v];
                }

                // Center diagonal block: diagD * I_5
                const FloatType dVal = diagD(i, j, k);
                for (size_t v = 0; v < 5; ++v) {
                    triplets.emplace_back(baseRow + v, baseRow + v, dVal);
                }

                // Helper to add a 5x5 block
                auto addBlock = [&](size_t neighborCellIdx, const Matrix2D<FloatType>& J, FloatType lam, bool isLower) {
                    const size_t baseCol = neighborCellIdx * 5;
                    for (size_t r = 0; r < 5; ++r) {
                        for (size_t c = 0; c < 5; ++c) {
                            FloatType val = 0.0;
                            if (isLower) {
                                // -0.5 * (A + lambda * I)
                                val = -0.5 * (J(r, c) + (r == c ? lam : 0.0));
                            } else {
                                // 0.5 * (A - lambda * I)
                                val = 0.5 * (J(r, c) - (r == c ? lam : 0.0));
                            }
                            triplets.emplace_back(baseRow + r, baseCol + c, val);
                        }
                    }
                };

                // Lower neighbor in I: (i - 1, j, k)
                if (i > 0) {
                    const Vector3D& S = surfacesI(i, j, k);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i - 1, j, k), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i - 1, j, k), S, soundSpeed(i - 1, j, k));
                    addBlock(cellIndex(i - 1, j, k), J, lam, true);
                }

                // Upper neighbor in I: (i + 1, j, k)
                if (i + 1 < _ni) {
                    const Vector3D& S = surfacesI(i + 1, j, k);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i + 1, j, k), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i + 1, j, k), S, soundSpeed(i + 1, j, k));
                    addBlock(cellIndex(i + 1, j, k), J, lam, false);
                }

                // Lower neighbor in J: (i, j - 1, k)
                if (!is1D && j > 0) {
                    const Vector3D& S = surfacesJ(i, j, k);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i, j - 1, k), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j - 1, k), S, soundSpeed(i, j - 1, k));
                    addBlock(cellIndex(i, j - 1, k), J, lam, true);
                }

                // Upper neighbor in J: (i, j + 1, k)
                if (!is1D && j + 1 < _nj) {
                    const Vector3D& S = surfacesJ(i, j + 1, k);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i, j + 1, k), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j + 1, k), S, soundSpeed(i, j + 1, k));
                    addBlock(cellIndex(i, j + 1, k), J, lam, false);
                }

                // Lower neighbor in K: (i, j, k - 1)
                if (is3D && k > 0) {
                    const Vector3D& S = surfacesK(i, j, k);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i, j, k - 1), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k - 1), S, soundSpeed(i, j, k - 1));
                    addBlock(cellIndex(i, j, k - 1), J, lam, true);
                }

                // Upper neighbor in K: (i, j, k + 1)
                if (is3D && k + 1 < _nk) {
                    const Vector3D& S = surfacesK(i, j, k + 1);
                    Matrix2D<FloatType> J = computeAdvectionJacobian(prim(i, j, k + 1), S, _fluid);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k + 1), S, soundSpeed(i, j, k + 1));
                    addBlock(cellIndex(i, j, k + 1), J, lam, false);
                }
            }
        }
    }

    Eigen::SparseMatrix<FloatType> A(numDofs, numDofs);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd sol;
    if (_solverType == LinearSolverType::GMRES) {
        Eigen::GMRES<Eigen::SparseMatrix<FloatType>, Eigen::DiagonalPreconditioner<FloatType>> solver;
        solver.set_restart(_restart);
        solver.setTolerance(_tolerance);
        solver.setMaxIterations(_maxIterations);
        solver.compute(A);
        sol = solver.solve(rhs);
    } else if (_solverType == LinearSolverType::FGMRES) {
        sol = solveFGMRES(A, rhs, prim, soundSpeed, diagD, gamma);
    } else {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<FloatType>, Eigen::DiagonalPreconditioner<FloatType>> solver;
        solver.setTolerance(_tolerance);
        solver.setMaxIterations(_maxIterations);
        solver.compute(A);
        sol = solver.solve(rhs);
    }

    // Unpack solution into deltaU
    deltaU.setToZero();
    #pragma omp parallel for collapse(2) schedule(static) if(_ni * _nj >= 64)
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                const size_t baseRow = cellIndex(i, j, k) * 5;
                StateVector dU;
                for (size_t v = 0; v < 5; ++v) {
                    dU[v] = sol(baseRow + v);
                }
                deltaU.set(i, j, k, dU);
            }
        }
    }
}

void ImplicitSolverKrylov::applyLUSGSPreconditioner(
    const Eigen::VectorXd& v,
    Eigen::VectorXd& z,
    const Matrix3D<StateVector>& prim,
    const Matrix3D<FloatType>& soundSpeed,
    const Matrix3D<FloatType>& diagD,
    FloatType gamma
) const {
    const bool is3D = (_topology == Topology::THREE_DIMENSIONAL);
    const bool is1D = (_topology == Topology::ONE_DIMENSIONAL);

    const auto& surfacesI = _mesh.getSurfacesI();
    const auto& surfacesJ = _mesh.getSurfacesJ();
    const auto& surfacesK = _mesh.getSurfacesK();

    auto cellIndex = [&](size_t i, size_t j, size_t k) -> size_t {
        return (i * _nj + j) * _nk + k;
    };

    _lusgsDeltaUStar.setToZero();
    _lusgsDeltaU.setToZero();

    // Forward Sweep: L * deltaU^* = v - LowerNeighbors
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                const size_t baseRow = cellIndex(i, j, k) * 5;
                StateVector rhs;
                for (size_t c = 0; c < 5; ++c) {
                    rhs[c] = v(baseRow + c);
                }

                if (i > 0) {
                    const Vector3D& S = surfacesI(i, j, k);
                    StateVector dU_m = _lusgsDeltaUStar.at(i - 1, j, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i - 1, j, k), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i - 1, j, k), S, soundSpeed(i - 1, j, k));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                if (!is1D && j > 0) {
                    const Vector3D& S = surfacesJ(i, j, k);
                    StateVector dU_m = _lusgsDeltaUStar.at(i, j - 1, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j - 1, k), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j - 1, k), S, soundSpeed(i, j - 1, k));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                if (is3D && k > 0) {
                    const Vector3D& S = surfacesK(i, j, k);
                    StateVector dU_m = _lusgsDeltaUStar.at(i, j, k - 1);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j, k - 1), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k - 1), S, soundSpeed(i, j, k - 1));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                FloatType d = diagD(i, j, k);
                _lusgsDeltaUStar.set(i, j, k, rhs / d);
            }
        }
    }

    // Backward Sweep: U * deltaU = D * deltaU^* - UpperNeighbors
    for (int i = static_cast<int>(_ni) - 1; i >= 0; --i) {
        for (int j = static_cast<int>(_nj) - 1; j >= 0; --j) {
            for (int k = static_cast<int>(_nk) - 1; k >= 0; --k) {
                FloatType d = diagD(i, j, k);
                StateVector rhs = _lusgsDeltaUStar.at(i, j, k) * d;

                if (static_cast<size_t>(i + 1) < _ni) {
                    const Vector3D& S = surfacesI(i + 1, j, k);
                    StateVector dU_p = _lusgsDeltaU.at(i + 1, j, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i + 1, j, k), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i + 1, j, k), S, soundSpeed(i + 1, j, k));
                    rhs += (AdU - dU_p * lam) * 0.5;
                }

                if (!is1D && static_cast<size_t>(j + 1) < _nj) {
                    const Vector3D& S = surfacesJ(i, j + 1, k);
                    StateVector dU_p = _lusgsDeltaU.at(i, j + 1, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j + 1, k), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j + 1, k), S, soundSpeed(i, j + 1, k));
                    rhs += (AdU - dU_p * lam) * 0.5;
                }

                if (is3D && static_cast<size_t>(k + 1) < _nk) {
                    const Vector3D& S = surfacesK(i, j, k + 1);
                    StateVector dU_p = _lusgsDeltaU.at(i, j, k + 1);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j, k + 1), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k + 1), S, soundSpeed(i, j, k + 1));
                    rhs += (AdU - dU_p * lam) * 0.5;
                }

                _lusgsDeltaU.set(i, j, k, rhs / d);

                const size_t baseRow = cellIndex(i, j, k) * 5;
                StateVector resVal = _lusgsDeltaU.at(i, j, k);
                for (size_t c = 0; c < 5; ++c) {
                    z(baseRow + c) = resVal[c];
                }
            }
        }
    }
}

Eigen::VectorXd ImplicitSolverKrylov::solveFGMRES(
    const Eigen::SparseMatrix<FloatType>& A,
    const Eigen::VectorXd& b,
    const Matrix3D<StateVector>& prim,
    const Matrix3D<FloatType>& soundSpeed,
    const Matrix3D<FloatType>& diagD,
    FloatType gamma
) {
    const size_t n = b.size();
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    FloatType bnorm = b.norm();
    if (bnorm < 1e-15) {
        return x;
    }

    size_t totalIters = 0;
    const size_t restartLen = (_restart > 0) ? _restart : 20;
    const size_t maxOuterIters = (_maxIterations + restartLen - 1) / restartLen;

    for (size_t outer = 0; outer < maxOuterIters && totalIters < _maxIterations; ++outer) {
        Eigen::VectorXd r = b - A * x;
        FloatType r0_norm = r.norm();
        if (r0_norm / bnorm < _tolerance) {
            return x;
        }

        size_t m = std::min(restartLen, _maxIterations - totalIters);
        std::vector<Eigen::VectorXd> V(m + 1);
        std::vector<Eigen::VectorXd> Z(m);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m + 1, m);
        std::vector<FloatType> cs(m, 1.0), sn(m, 0.0);
        Eigen::VectorXd g = Eigen::VectorXd::Zero(m + 1);

        V[0] = r / r0_norm;
        g(0) = r0_norm;

        size_t k = 0;
        for (size_t j = 0; j < m; ++j) {
            totalIters++;
            k = j;

            // Apply right preconditioner: z_j = M^{-1} v_j
            Z[j].resize(n);
            if (_precondType == LinearPreconditionerType::LU_SGS) {
                applyLUSGSPreconditioner(V[j], Z[j], prim, soundSpeed, diagD, gamma);
            } else if (_precondType == LinearPreconditionerType::DIAGONAL) {
                for (size_t i = 0; i < _ni; ++i) {
                    for (size_t j_idx = 0; j_idx < _nj; ++j_idx) {
                        for (size_t k_idx = 0; k_idx < _nk; ++k_idx) {
                            size_t baseRow = ((i * _nj + j_idx) * _nk + k_idx) * 5;
                            FloatType dVal = diagD(i, j_idx, k_idx);
                            for (size_t v = 0; v < 5; ++v) {
                                Z[j](baseRow + v) = V[j](baseRow + v) / dVal;
                            }
                        }
                    }
                }
            } else {
                Z[j] = V[j];
            }

            // Matrix-vector product: w = A * z_j
            Eigen::VectorXd w = A * Z[j];

            // Modified Gram-Schmidt orthogonalization
            for (size_t i = 0; i <= j; ++i) {
                H(i, j) = V[i].dot(w);
                w -= H(i, j) * V[i];
            }
            H(j + 1, j) = w.norm();
            if (H(j + 1, j) > 1e-15) {
                V[j + 1] = w / H(j + 1, j);
            }

            // Apply prior Givens rotations to column j
            for (size_t i = 0; i < j; ++i) {
                FloatType temp = cs[i] * H(i, j) + sn[i] * H(i + 1, j);
                H(i + 1, j) = -sn[i] * H(i, j) + cs[i] * H(i + 1, j);
                H(i, j) = temp;
            }

            // Compute new Givens rotation
            FloatType h1 = H(j, j);
            FloatType h2 = H(j + 1, j);
            FloatType denom = std::hypot(h1, h2);
            if (denom < 1e-15) {
                cs[j] = 1.0;
                sn[j] = 0.0;
            } else {
                cs[j] = h1 / denom;
                sn[j] = h2 / denom;
            }

            // Eliminate subdiagonal
            H(j, j) = cs[j] * h1 + sn[j] * h2;
            H(j + 1, j) = 0.0;

            // Apply rotation to g
            FloatType g1 = g(j);
            g(j) = cs[j] * g1;
            g(j + 1) = -sn[j] * g1;

            FloatType res_norm = std::abs(g(j + 1));
            if (res_norm / bnorm < _tolerance) {
                break;
            }
        }

        // Solve upper triangular system H(0..k, 0..k) * y = g(0..k)
        Eigen::VectorXd y(k + 1);
        for (int i = static_cast<int>(k); i >= 0; --i) {
            FloatType sum = g(i);
            for (size_t j_col = i + 1; j_col <= k; ++j_col) {
                sum -= H(i, j_col) * y(j_col);
            }
            y(i) = sum / H(i, i);
        }

        // Update x = x + sum_{i=0}^k y_i * Z_i
        for (size_t i = 0; i <= k; ++i) {
            x += y(i) * Z[i];
        }

        FloatType cur_err = (b - A * x).norm() / bnorm;
        if (cur_err < _tolerance) {
            return x;
        }
    }

    return x;
}
