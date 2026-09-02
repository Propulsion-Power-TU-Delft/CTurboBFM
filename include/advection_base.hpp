#pragma once

#include "types.hpp"
#include "fluid_base.hpp"
#include "config.hpp"

class AdvectionBase {
public:

    AdvectionBase(const Config &config, const FluidBase &fluid);

    virtual ~AdvectionBase() = default;  

    /**
     * @param Ull Second cell to the left of the interface.
     * @param Ul  First cell to the left of the interface.
     * @param Ur  First cell to the right of the interface.
     * @param Urr Second cell to the right of the interface.
     * @param surface The surface normal vector at the cell interface (oriented from left to right).
     */
    virtual StateVector computeFlux(const StateVector &Ull, 
        const StateVector &Ul, 
        const StateVector &Ur, 
        const StateVector &Urr,
        const Vector3D &surface) const = 0;
    
    void musclReconstructLeftRight(
        StateVector& Wll, 
        StateVector& Wl, 
        StateVector& Wr, 
        StateVector& Wrr, 
        FluxLimiter fluxLimiter) const;

    StateVector computeReconstructionSmoothness(
        const StateVector& Wl, 
        const StateVector& Wc, 
        const StateVector& Wr) const;

    StateVector computeLimiter(const StateVector& smoothness) const;

protected:
    const FluidBase& _fluid;  
    const Config& _config; 
    bool _isMusclActive = false; 
    FluxLimiter _fluxLimiter = FluxLimiter::NONE; 

};
