#include "advection_base.hpp"

AdvectionBase::AdvectionBase(const Config &config, const FluidBase &fluid)
    : _fluid(fluid), _config(config) 
{
    _isMusclActive = _config.getMUSCLreconstruction();
    _fluxLimiter = _config.getFluxLimiter();
}

void AdvectionBase::musclReconstructLeftRight(StateVector& Wll, 
    StateVector& Wl, 
    StateVector& Wr, 
    StateVector& Wrr, 
    FluxLimiter fluxLimiter) const {
    
    StateVector smoothnessL = computeReconstructionSmoothness(Wll, Wl, Wr);
    StateVector smoothnessR = computeReconstructionSmoothness(Wl, Wr, Wrr);

    StateVector limiterL = computeLimiter(smoothnessL);
    StateVector limiterR = computeLimiter(smoothnessR);

    Wl += limiterL * (Wr - Wl) * 0.5;
    Wr -= limiterR * (Wrr - Wr) * 0.5;
}


StateVector AdvectionBase::computeReconstructionSmoothness(
    const StateVector& Wl, 
    const StateVector& Wc, 
    const StateVector& Wr) const {
    
    StateVector smoothness;
    
    for (std::size_t i = 0; i < 5; ++i) {
        FloatType dR = Wr[i] - Wc[i];
        FloatType dL = Wc[i] - Wl[i];
        smoothness[i] = (std::abs(dR) < 1E-12) ? 0.0 : dL / dR;
    }

    return smoothness;
}

StateVector AdvectionBase::computeLimiter(const StateVector& smoothness) const{
    StateVector limiter({0.0, 0.0, 0.0, 0.0, 0.0});

    for (std::size_t i = 0; i < smoothness.size(); ++i) {
        switch (_fluxLimiter) {
            case FluxLimiter::NONE:
                limiter[i] = 0.0;
                break;
            case FluxLimiter::VAN_ALBADA:
                limiter[i] = (smoothness[i] <= 0.0) ? 0.0 :
                    (smoothness[i]*smoothness[i] + smoothness[i]) / (1.0 + smoothness[i]*smoothness[i]);
                break;
            case FluxLimiter::VAN_LEER:
                limiter[i] = (smoothness[i] + std::abs(smoothness[i])) / (1.0 + std::abs(smoothness[i]));
                break;
            case FluxLimiter::MIN_MOD:
                limiter[i] = std::max(0.0, std::min(1.0, smoothness[i]));
                break;
        }
    }

    return limiter;
}