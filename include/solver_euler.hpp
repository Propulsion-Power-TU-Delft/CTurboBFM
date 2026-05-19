#pragma once
#include "types.hpp"
#include "solver_base.hpp"
#include "output_base.hpp"
#include "output_csv.hpp"
#include "source_bfm_base.hpp"
#include "source_bfm_hall.hpp"
#include "source_bfm_thollet.hpp"
#include "source_bfm_chima.hpp"
#include "source_bfm_correlations.hpp"
#include "source_bfm_lift_drag.hpp"
#include "greitzer_model.hpp"   
#include "turbulence_model_base.hpp"
#include "turbulence_model_sa.hpp"
#include <unordered_map>

class SolverEuler : public SolverBase {

public:

    SolverEuler(Config& config, Mesh& mesh);

    ~SolverEuler() = default;

    void solve() override;

    virtual void checkConvergence(bool &exitLoop, bool &isSteady) const override;

    void writeSolution(size_t iterationCounter, bool alsoGradients=false) override;


private:

    void initializeSolutionArrays() override;

    void computeTimestepArray(const FlowSolution &solution, Matrix3D<FloatType> &timestep);

    void updateMassFlows(const FlowSolution &solution);

    void updateTurboPerformance(const FlowSolution &solution);

    void updateMonitorPoints(const FlowSolution &solution);

    void initializeMonitorPoints();


    /** @brief compute the global residual, defined as V*dU/dt + R = 0 -> R = fluxes - source terms */
    void computeResiduals(
        FlowSolution& solution, 
        const std::map<SolutionName, Matrix3D<Vector3D>>& solutionGrad, 
        const size_t itCounter, 
        const FloatType timePhysical, 
        Matrix3D<FloatType> &timestep,
        FlowSolution &residuals);
    
    /** @brief compute the advection flux contribution to the residuals */
    void computeAdvectionFluxResiduals(
        FluxDirection direction, 
        const FlowSolution& solution, 
        size_t itCounter, 
        FlowSolution &residuals) const;
    
    /** @brief compute the advection flux contribution to the residuals */
    void computeViscousFluxResiduals(
        FluxDirection direction, 
        const FlowSolution& solution, 
        const std::map<SolutionName, Matrix3D<Vector3D>>& solutionGrad, 
        size_t itCounter, 
        FlowSolution &residuals) const;

    void updateSolution(
        const FlowSolution &solOld, 
        FlowSolution &solNew, 
        const FlowSolution &residuals, 
        const FloatType &integrationCoeff, 
        const Matrix3D<FloatType> &dt);

    void enforcePeriodicityOnSolution(FlowSolution &sol);

    void computeSolutionGradient(FlowSolution &sol, std::map<SolutionName, Matrix3D<Vector3D>> &solutionGrad);

    void printInfoResiduals(FlowSolution &residuals, size_t it);

    void initializeSolutionFromScratch();

    void initializeSolutionFromRestart();

    void printCheckOfMassFlowConservation(size_t it) const;

    void printTurboPerformance(size_t it) const;

    void printLogResidualsHeader() const;

    StateVector computeLogResidualNorm(const FlowSolution &residuals) const;

    void printLogResiduals(const StateVector &logRes, unsigned long int it) const;

    void writeLogResidualsToCsvFile() const;

    void writeTurboPerformanceToCsvFile() const;

    void writeGreitzerDynamicsToCsvFile() const;

    void writeMonitorPointsToCsvFile() const;

    void updateRadialProfiles(FlowSolution &solution);
    
    /** @brief compute all the source terms and apply them to the residuals 
     * Contains the BFM source terms, axisymmetric geometric terms, and Gong BF formulation terms
    */
    void computeSourceResiduals(
        FlowSolution& solution, 
        const std::map<SolutionName, Matrix3D<Vector3D>>& solutionGrad, 
        size_t itCounter, 
        FlowSolution &residuals, 
        Matrix3D<Vector3D> &inviscidForce, 
        Matrix3D<Vector3D> &viscousForce, 
        Matrix3D<FloatType> &deviationAngle, 
        FloatType timePhysical,
        Matrix3D<FloatType> &timestep);

    void readRestartFile(
        const std::string &restartFileName, 
        size_t &NI, 
        size_t &NJ, 
        size_t &NK,
        Matrix3D<FloatType> &inputDensity, 
        Matrix3D<FloatType> &inputVelX, 
        Matrix3D<FloatType> &inputVelY, 
        Matrix3D<FloatType> &inputVelZ, 
        Matrix3D<FloatType> &inputTemperature, 
        Matrix3D<Vector3D> &inputForceViscous, 
        Matrix3D<Vector3D> &inputForceInviscid);
    
    void axisymmetricRestart(
        Matrix3D<FloatType> &inputDensity, 
        Matrix3D<FloatType> &inputVelX, 
        Matrix3D<FloatType> &inputVelY, 
        Matrix3D<FloatType> &inputVelZ, 
        Matrix3D<FloatType> &inputTemperature);

    void standardRestart(
        Matrix3D<FloatType> &inputDensity, 
        Matrix3D<FloatType> &inputVelX, 
        Matrix3D<FloatType> &inputVelY, 
        Matrix3D<FloatType> &inputVelZ, 
        Matrix3D<FloatType> &inputTemperature, 
        Matrix3D<Vector3D> &inputForceViscous, 
        Matrix3D<Vector3D> &inputForceInviscid);

    void computeGradientOfField(const Matrix3D<FloatType> &var, Matrix3D<Vector3D> &grad) const;
    
    /** @brief compute the source terms related to Gong body force formulations 
     * The term computed here is simply -Omega*dudtheta. It is computed with upwinded 2nd order finite differences
     * in the direction of shaft rotation
    */
    StateVector computeGongSource(
        const FloatType& radius, 
        const FloatType& theta, 
        const FloatType& omega, 
        const size_t i, 
        const size_t j, 
        const size_t k, 
        const FloatType& volume) const;
    
    /** @brief set the momentum vector on the nodes belonging to viscous walls */
    void setMomentumSolutionOnViscousWalls(
        FlowSolution &residuals, 
        const BoundaryIndex &BoundaryIndex, 
        const Vector3D &wallVelocity) const;
    
    /** @brief perform some preprocessing of the solution */
    void preprocessSolution(FlowSolution &solution, bool updateRadialProf = true);

    StateVector computeViscousFlux(
        const StateVector& conservative, 
        const Vector3D& velXGrad, 
        const Vector3D& velYGrad, 
        const Vector3D& velZGrad, 
        const Vector3D& tempGrad, 
        const Vector3D& surface) const;

protected:
    void enforcePeriodicityOnResiduals(FlowSolution& residuals, FloatType& angleRad) const;

    void enforceNoSlipWallsOnResiduals(FlowSolution& residuals) const;

    void understandWhatSourcesAreNeeded(
        const size_t i, 
        const size_t j, 
        const size_t k, 
        bool &geometricSourceFlag, 
        bool &gongSourceFlag) const;

private:
    FlowSolution _conservativeSolution; 
    std::map<SolutionName, Matrix3D<Vector3D>> _solutionGrad;
    std::unique_ptr<OutputBase> _output;
    
    bool _isBfmActive{false}, _isTurbulenceActive{false};
    bool _isGongFormulationActive{false};
    std::unique_ptr<SourceBFMBase> _bfmSource;

    std::unique_ptr<TurbulenceModelBase> _turbulenceModel;
    
    std::vector<StateVector> _logResiduals;
    
    Matrix3D<Vector3D> _inviscidForce, _viscousForce;
    Matrix3D<FloatType> _deviationAngle;

    std::vector<unsigned int> _monitorPointsIdxI, _monitorPointsIdxJ, _monitorPointsIdxK;
    unsigned int _numberMonitorPoints = 0;

    std::unique_ptr<GreitzerModel> _greitzerModel;
    bool _isGreitzerModelingActive{false};
    
};
