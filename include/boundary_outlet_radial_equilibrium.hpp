#pragma once
#include "boundary_base.hpp"

class BoundaryOutletRadialEquilibrium : public BoundaryBase {
    
public:
    
    /** 
     * @brief Constructs the object with given references.
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param boundIndex The boundary index.
     * @param radialProfile reference to the radial profile of pressure used in the boundary condition flux evaluation.
     */
    BoundaryOutletRadialEquilibrium(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::vector<FloatType>& pressure)
        : BoundaryBase(config, mesh, fluid), 
        _radialPressureProfile(pressure){}        
        
    virtual ~BoundaryOutletRadialEquilibrium() {}
    
    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
    
protected:
    std::vector<FloatType>& _radialPressureProfile; // Reference to the radial profile of the static pressure (which is updated at every iteration)
};