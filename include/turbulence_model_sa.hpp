#pragma once

#include "turbulence_model_base.hpp"

class TurbulenceModelSA : public TurbulenceModelBase {
public:
    explicit TurbulenceModelSA(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        const std::vector<Boundary> &boundaries,
        const Matrix3D<FloatType> &wallDistance
    );
    
    ~TurbulenceModelSA() override = default;

protected:
    void setupModelEquations() override;
    
    void setupBoundaryValues(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) override;
    
    void setupInitialValues() override;
    
    void enforceNoSlipWallCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) override;
    
    void enforceInletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) override;
    
    void enforceOutletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) override;
    
    virtual void solve(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad) override;

private:
    const FloatType _cb1 = 0.1355;
    const FloatType _cb2 = 0.622;
    const FloatType _sigma = 2/3;
    const FloatType _kappa = 0.41;
    const FloatType _cw2 = 0.3;
    const FloatType _cw3 = 2.0;
    const FloatType _cv1 = 7.1;
    const FloatType _ct1 = 1.0;
    const FloatType _ct2 = 2.0;
    const FloatType _ct3 = 1.25; // nasa 1.2, blazek 1.3
    const FloatType _ct4 = 0.5;
    const FloatType _cw1 = _cb1/(_kappa * _kappa) + (1.0 + _cb2) / _sigma;
    const Matrix3D<FloatType> &_wallDistance;
    Matrix3D<FloatType> _nuHat;
};