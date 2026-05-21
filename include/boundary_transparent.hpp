#pragma once
#include "boundary_base.hpp"
#include "advection_base.hpp"

/** 
  *  \brief     Class handling transparent boundary conditions (1D only).
  */
class BoundaryTransparent : public BoundaryBase {
    
public:

    BoundaryTransparent(const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        AdvectionBase &advScheme)
        : BoundaryBase(config, mesh, fluid), 
        _advScheme(advScheme) {}            

    virtual ~BoundaryTransparent() {}

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;

private:
    AdvectionBase &_advScheme;
        
};