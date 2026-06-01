#pragma once

#include "turbulence_model_base.hpp"

class TurbulenceModelNone : public TurbulenceModelBase {
public:
    explicit TurbulenceModelNone(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        const std::vector<Boundary> &boundaries,
        const Matrix3D<FloatType> &wallDistance,
        const FlowSolution &initialSolution
    ) : TurbulenceModelBase(config, fluid, mesh, boundaries, initialSolution) {};
    
    ~TurbulenceModelNone() override = default;

protected:
    void setupModelEquations() override {};
    
    void updateBoundaryValues(
        const FlowSolution &sol) override {};
    
    void setupInitialValues() override {};
    
    void enforceNoSlipWallCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override {};
    
    void enforceOutletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override {};
    
    void enforceInletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override {};
    
    void solve(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad,
        const Matrix3D<FloatType> &dt) override {};
    
    void updateSolutionGradient() override {};

    void updateMeanFlowTerms(const FlowSolution &sol) {};

    void computeFluxTerms(Matrix3D<FloatType> &residual, const FlowSolution &sol) {};

    FloatType getEddyViscosity(const FloatType &density, size_t i, size_t j, size_t k) const override {return 0.0;}

    FloatType getEddyThermalConductivity(const FloatType &mu, const FloatType &cp, const FloatType &Pr) const override {return 0.0;}

private:
    
};