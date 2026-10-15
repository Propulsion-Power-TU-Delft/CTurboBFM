#pragma once

#include "boundary_base.hpp"

/** 
  *  \brief     Class handling outlet boundary with prescribed mass flow rate condition.
  *  \details   It computes the boundary fluxes using characteristic outlet formulation with dynamically adjusted pressure.
  *  \author    Francesco Neri
  */
class BoundaryOutletMassFlow : public BoundaryBase {
    
public:
    
    BoundaryOutletMassFlow(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType>& pressure)
        : BoundaryBase(config, mesh, fluid), 
        _pressureProfile(pressure){}            

    virtual ~BoundaryOutletMassFlow() = default;
    
    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
    
protected:
    std::vector<FloatType>& _pressureProfile; 

};
