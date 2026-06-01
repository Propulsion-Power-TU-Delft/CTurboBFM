#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "boundary_base.hpp"
#include "math_utils.hpp"

class TurbulenceModelBase {
public:
    explicit TurbulenceModelBase(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        const std::vector<Boundary> &boundaries,
        const FlowSolution &initialSolution
    );
    
    virtual ~TurbulenceModelBase() = default;

    virtual void setupModelEquations() = 0;

    virtual void setupInitialValues() = 0;
    
    virtual void updateBoundaryValues(const FlowSolution &sol) = 0;
    
    virtual void enforceNoSlipWallCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) = 0;
    
    virtual void enforceOutletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) = 0;
    
    virtual void enforceInletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) = 0;

    virtual void solve(
        const FlowSolution &sol,
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad, 
        const Matrix3D<FloatType> &dt) = 0;
    
    virtual void updateSolutionGradient() = 0;

    const std::array<int, 3> getStepMask(FluxDirection direction) const;

    virtual FloatType getEddyViscosity(const FloatType &density, size_t i, size_t j, size_t k) const = 0;

    virtual FloatType getEddyThermalConductivity(const FloatType &mu, const FloatType &cp, const FloatType &Pr) const = 0;

protected:
    const Config& _config;
    const FluidBase& _fluid;
    const Mesh& _mesh;
    const std::vector<Boundary>& _boundaries;
    const FlowSolution& _initialSolution;
    size_t _ni, _nj, _nk;
    Topology _topology;
};
