#pragma once
#include "boundary_base.hpp"

/**
* @brief Compute the boundary flux.
* Formulation taken from 'Formulation and Implementation of Inflow/Outflow Boundary Conditions 
to Simulate Propulsive Effects', Rodriguez et al.
*/
class BoundaryOutlet : public BoundaryBase{
    
public:

    /**
     * @brief Constructs the object with given references.
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param boundIndex The boundary index.
     * @param bcValues The outlet values (static pressure).
     */
    BoundaryOutlet(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType> bcValues)
        : BoundaryBase(config, mesh, fluid),
        _boundaryValues(bcValues) {}
        
    virtual ~BoundaryOutlet() {}

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