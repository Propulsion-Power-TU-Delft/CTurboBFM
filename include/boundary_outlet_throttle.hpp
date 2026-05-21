#pragma once

#include "boundary_base.hpp"

/** 
  *  \brief     Class handling outlet boundary with throttle condition.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class BoundaryOutletThrottle : public BoundaryBase {
    
public:
    
    BoundaryOutletThrottle(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType>& pressure)
        : BoundaryBase(config, mesh, fluid), 
        _radialPressureProfile(pressure){}            

    virtual ~BoundaryOutletThrottle() {}
    
    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
    
protected:
    std::vector<FloatType>& _radialPressureProfile; 

};