#pragma once

#include "turbulence_model_base.hpp"

class TurbulenceModelSA : public TurbulenceModelBase {
public:
    explicit TurbulenceModelSA(
        const Config &config, 
        const FluidBase &fluid, 
        const Mesh &mesh, 
        const std::vector<Boundary> &boundaries,
        const Matrix3D<FloatType> &wallDistance,
        const FlowSolution &initialSolution
    );
    
    ~TurbulenceModelSA() override = default;

protected:
    void setupModelEquations() override;
    
    void updateBoundaryValues(
        const FlowSolution &sol) override;
    
    void setupInitialValues() override;
    
    void enforceNoSlipWallCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override;
    
    void enforceOutletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override;
    
    void enforceInletCondition(
        const Boundary& boundary, 
        const FlowSolution &sol) override;
    
    void solve(
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad,
        const Matrix3D<FloatType> &dt) override;
    
    void updateSolutionGradient() override;

    void updateMeanFlowTerms(const FlowSolution &sol);

    void computeFluxTerms(Matrix3D<FloatType> &residual, const FlowSolution &sol);

    FloatType getEddyViscosity(const FloatType &density, size_t i, size_t j, size_t k) const override;

    void computeSourceTerms(
        Matrix3D<FloatType> &residual, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad);
    
    FloatType computeSource(
        const std::array<size_t, 3> &idx, 
        const FlowSolution &sol, 
        const std::map<SolutionName, Matrix3D<Vector3D>> &solGrad);

    void computeFluxContribution(FluxDirection direction, const FlowSolution &sol, Matrix3D<FloatType> &residual);

    FloatType computeAdvectionFlux(
        const std::array<size_t, 3> &idxLeft, 
        const std::array<size_t, 3> &idxRight,
        const Vector3D &surface,
        const FlowSolution &sol);
    
    FloatType computeViscousFlux(
        const std::array<size_t, 3> &idxLeft, 
        const std::array<size_t, 3> &idxRight,
        const Vector3D &surface);
    
    void initializeFromZero();

    void initializeFromRestartFile();
    
    void updateSolution(const Matrix3D<FloatType> &residual, const Matrix3D<FloatType> &dt);

    FloatType getEddyThermalConductivity(const FloatType &mu, const FloatType &cp, const FloatType &Pr) const override {return mu*cp/Pr;}

private:
    const FloatType _cb1 = 0.1355;
    const FloatType _cb2 = 0.622;
    const FloatType _sigma = 2.0/3.0;
    const FloatType _kappa = 0.41;
    const FloatType _cw2 = 0.3;
    const FloatType _cw3 = 2.0;
    const FloatType _cv1 = 7.1;
    const FloatType _cv2 = 5.0;
    const FloatType _ct3 = 1.25; // nasa 1.2, blazek 1.3
    const FloatType _ct4 = 0.5;
    const FloatType _cw1 = _cb1/(_kappa * _kappa) + (1.0 + _cb2) / _sigma;
    const Matrix3D<FloatType> &_wallDistance;
    
    Matrix3D<FloatType> _nuHat;
    Matrix3D<FloatType> _nuLaminar;
    Matrix3D<Vector3D> _nuHatGrad;
    Matrix3D<FloatType> _fv1;
    Matrix3D<FloatType> _initNuHat;
};