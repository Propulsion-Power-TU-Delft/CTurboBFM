#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "boundary_base.hpp"

class TurbulenceModelBase {
public:
    explicit TurbulenceModelBase(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        const std::vector<Boundary> &boundaries
    );
    
    virtual ~TurbulenceModelBase() = default;

    virtual void setupModelEquations() = 0;

    virtual void setupInitialValues() = 0;
    
    virtual void setupBoundaryValues(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) = 0;
    
    virtual void enforceNoSlipWallCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) = 0;
    
    virtual void enforceInletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) = 0;
    
    virtual void enforceOutletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) = 0;

    virtual void solve(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) = 0;

protected:
    const Config& _config;
    const FluidBase& _fluid;
    const Mesh& _mesh;
    const std::vector<Boundary>& _boundaries;
    size_t _ni, _nj, _nk;
};
