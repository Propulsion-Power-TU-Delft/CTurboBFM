#pragma once
#include "boundary_base.hpp"

class BoundaryOutletSupersonic : public BoundaryBase {
    
public:

    BoundaryOutletSupersonic(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType> inletValues)
        : BoundaryBase(config, mesh, fluid), 
        _boundaryValues(inletValues) {}            

    virtual ~BoundaryOutletSupersonic() = default;

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