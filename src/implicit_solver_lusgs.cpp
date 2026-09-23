#include "implicit_solver_lusgs.hpp"
#include "turbulence_model_base.hpp"
#include <cmath>

ImplicitSolverLUSGS::ImplicitSolverLUSGS(
    const Config& config, 
    const Mesh& mesh, 
    const FluidBase& fluid,
    const TurbulenceModelBase* turbModel
) : ImplicitSolverBase(config, mesh, fluid, turbModel),
    _deltaUStar(mesh.getNumberPointsI(), mesh.getNumberPointsJ(), mesh.getNumberPointsK()) {
}

void ImplicitSolverLUSGS::solveCorrection(
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

    // Precompute primitive variables, sound speed, and diagonal matrix D
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

                // Face spectral radii
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

    _deltaUStar.setToZero();
    deltaU.setToZero();

    // -------------------------------------------------------------
    // Forward Sweep: L * deltaU^* = -R - LowerNeighbors
    // Order: i increasing, j increasing, k increasing
    // -------------------------------------------------------------
    for (size_t i = 0; i < _ni; ++i) {
        for (size_t j = 0; j < _nj; ++j) {
            for (size_t k = 0; k < _nk; ++k) {
                StateVector rhs = -residuals.at(i, j, k);

                // Lower neighbor in I: cell (i-1, j, k) across face (i, j, k)
                if (i > 0) {
                    const Vector3D& S = surfacesI(i, j, k);
                    StateVector dU_m = _deltaUStar.at(i - 1, j, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i - 1, j, k), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i - 1, j, k), S, soundSpeed(i - 1, j, k));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                // Lower neighbor in J: cell (i, j-1, k) across face (i, j, k)
                if (!is1D && j > 0) {
                    const Vector3D& S = surfacesJ(i, j, k);
                    StateVector dU_m = _deltaUStar.at(i, j - 1, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j - 1, k), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j - 1, k), S, soundSpeed(i, j - 1, k));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                // Lower neighbor in K: cell (i, j, k-1) across face (i, j, k)
                if (is3D && k > 0) {
                    const Vector3D& S = surfacesK(i, j, k);
                    StateVector dU_m = _deltaUStar.at(i, j, k - 1);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j, k - 1), S, dU_m, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k - 1), S, soundSpeed(i, j, k - 1));
                    rhs += (AdU + dU_m * lam) * 0.5;
                }

                FloatType invD = 1.0 / diagD(i, j, k);
                _deltaUStar.set(i, j, k, rhs * invD);
            }
        }
    }

    // -------------------------------------------------------------
    // Backward Sweep: (D + U) * deltaU = D * deltaU^*
    // deltaU = deltaU^* + D^-1 * UpperNeighbors
    // Order: i decreasing, j decreasing, k decreasing
    // -------------------------------------------------------------
    for (size_t i = _ni; i-- > 0; ) {
        for (size_t j = _nj; j-- > 0; ) {
            for (size_t k = _nk; k-- > 0; ) {
                StateVector sum{};

                // Upper neighbor in I: cell (i+1, j, k) across face (i+1, j, k)
                if (i + 1 < _ni) {
                    const Vector3D& S = surfacesI(i + 1, j, k);
                    StateVector dU_p = deltaU.at(i + 1, j, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i + 1, j, k), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i + 1, j, k), S, soundSpeed(i + 1, j, k));
                    sum += (dU_p * lam - AdU) * 0.5;
                }

                // Upper neighbor in J: cell (i, j+1, k) across face (i, j+1, k)
                if (!is1D && j + 1 < _nj) {
                    const Vector3D& S = surfacesJ(i, j + 1, k);
                    StateVector dU_p = deltaU.at(i, j + 1, k);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j + 1, k), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j + 1, k), S, soundSpeed(i, j + 1, k));
                    sum += (dU_p * lam - AdU) * 0.5;
                }

                // Upper neighbor in K: cell (i, j, k+1) across face (i, j, k+1)
                if (is3D && k + 1 < _nk) {
                    const Vector3D& S = surfacesK(i, j, k + 1);
                    StateVector dU_p = deltaU.at(i, j, k + 1);
                    StateVector AdU = evaluateJacobianVectorProduct(prim(i, j, k + 1), S, dU_p, gamma);
                    FloatType lam = computeFaceSpectralRadius(prim(i, j, k + 1), S, soundSpeed(i, j, k + 1));
                    sum += (dU_p * lam - AdU) * 0.5;
                }

                FloatType invD = 1.0 / diagD(i, j, k);
                deltaU.set(i, j, k, _deltaUStar.at(i, j, k) + sum * invD);
            }
        }
    }
}
