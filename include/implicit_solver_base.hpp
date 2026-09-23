#pragma once

#include "types.hpp"
#include "mesh.hpp"
#include "config.hpp"
#include "fluid_base.hpp"
#include <memory>

class TurbulenceModelBase;

/**
 * @brief Base class for implicit linear system solvers / time-stepping operators.
 * Solves the correction linear system: A * deltaU = - residuals
 */
class ImplicitSolverBase {
public:
    ImplicitSolverBase(
        const Config& config, 
        const Mesh& mesh, 
        const FluidBase& fluid,
        const TurbulenceModelBase* turbModel = nullptr
    ) : _config(config), _mesh(mesh), _fluid(fluid), _turbModel(turbModel),
        _ni(mesh.getNumberPointsI()),
        _nj(mesh.getNumberPointsJ()),
        _nk(mesh.getNumberPointsK()),
        _topology(config.getTopology()) {}

    virtual ~ImplicitSolverBase() = default;

    void setTurbulenceModel(const TurbulenceModelBase* turbModel) {
        _turbModel = turbModel;
    }

    /**
     * @brief Solve for the conservative variable correction deltaU.
     * @param solution Current conservative flow solution U^n
     * @param residuals Residuals R(U^n) computed by computeResiduals()
     * @param dt Local pseudo-time step array (dt)
     * @param deltaU Output correction array deltaU = U^{n+1} - U^n
     */
    virtual void solveCorrection(
        const FlowSolution& solution,
        const FlowSolution& residuals,
        const Matrix3D<FloatType>& dt,
        FlowSolution& deltaU
    ) = 0;

protected:
    const Config& _config;
    const Mesh& _mesh;
    const FluidBase& _fluid;
    const TurbulenceModelBase* _turbModel{nullptr};
    size_t _ni{0}, _nj{0}, _nk{0};
    Topology _topology;

    static FloatType computeFaceSpectralRadius(
        const StateVector& primitive,
        const Vector3D& S,
        FloatType soundSpeed
    ) {
        const FloatType u = primitive[1];
        const FloatType v = primitive[2];
        const FloatType w = primitive[3];
        const FloatType sMag = S.magnitude();
        if (sMag < 1e-15) return 0.0;

        const FloatType vn = std::abs(u * S.x() + v * S.y() + w * S.z()) / sMag;
        return (vn + soundSpeed) * sMag;
    }

    static StateVector evaluateJacobianVectorProduct(
        const StateVector& primitive,
        const Vector3D& S,
        const StateVector& dU,
        FloatType gamma
    ) {
        const FloatType u = primitive[1];
        const FloatType v = primitive[2];
        const FloatType w = primitive[3];
        const FloatType et = primitive[4];

        const FloatType nx = S.x();
        const FloatType ny = S.y();
        const FloatType nz = S.z();

        const FloatType umag2 = u*u + v*v + w*w;
        const FloatType V = u*nx + v*ny + w*nz;
        const FloatType phi = 0.5 * (gamma - 1.0) * umag2;
        const FloatType gm1 = gamma - 1.0;

        StateVector res;

        // Row 0: continuity
        res[0] = nx * dU[1] + ny * dU[2] + nz * dU[3];

        // Row 1: x-momentum
        res[1] = ((phi - u*u)*nx - u*v*ny - u*w*nz) * dU[0]
               + ((3.0 - gamma)*u*nx + V) * dU[1]
               + (-gm1*v*nx + u*ny) * dU[2]
               + (-gm1*w*nx + u*nz) * dU[3]
               + (gm1*nx) * dU[4];

        // Row 2: y-momentum
        res[2] = (-u*v*nx + (phi - v*v)*ny - v*w*nz) * dU[0]
               + (v*nx - gm1*u*ny) * dU[1]
               + (u*nx + (3.0 - gamma)*v*ny + w*nz) * dU[2]
               + (-gm1*w*ny + v*nz) * dU[3]
               + (gm1*ny) * dU[4];

        // Row 3: z-momentum
        res[3] = (-u*w*nx - v*w*ny + (phi - w*w)*nz) * dU[0]
               + (w*nx - gm1*u*nz) * dU[1]
               + (w*ny - gm1*v*nz) * dU[2]
               + (u*nx + v*ny + (3.0 - gamma)*w*nz) * dU[3]
               + (gm1*nz) * dU[4];

        // Row 4: energy
        const FloatType h0_term = gamma * et - 2.0 * phi;
        const FloatType et_phi = gamma * et - phi;

        res[4] = (-V * h0_term) * dU[0]
               + (et_phi * nx - gm1 * u * V) * dU[1]
               + (et_phi * ny - gm1 * v * V) * dU[2]
               + (et_phi * nz - gm1 * w * V) * dU[3]
               + (gamma * V) * dU[4];

        return res;
    }
};
