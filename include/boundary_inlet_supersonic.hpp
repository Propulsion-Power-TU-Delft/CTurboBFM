#pragma once
#include "boundary_base.hpp"

/**
 * Formulation taken from 'Formulation and Implementation of Inflow/Outflow Boundary Conditions to Simulate 
 * Propulsive Effects', Rodriguez et al.
 */
class BoundaryInletSupersonic : public BoundaryBase {
    
public:
    BoundaryInletSupersonic(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType> inletValues)
        : BoundaryBase(config, mesh, fluid),
        _boundaryValues(inletValues) {}
        
    virtual ~BoundaryInletSupersonic() {}

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
    
protected:
    std::vector<FloatType> _boundaryValues;
};