#pragma once

#include "turbulence_model_base.hpp"

class TurbulenceModelSA : public TurbulenceModelBase {
public:
    explicit TurbulenceModelSA(const Config &config, const FluidBase &fluid, const Mesh &mesh);
    
    ~TurbulenceModelSA() override = default;

private:

    const FloatType _cb1 = 0.1355;
    const FloatType _cb2 = 0.622;
    const FloatType _sigma = 2/3;
    const FloatType _kappa = 0.41;
    const FloatType _cw2 = 0.3;
    const FloatType _cw3 = 2.0;
    const FloatType _cv1 = 7.1;
    const FloatType _ct3 = 1.2;
    const FloatType _ct4 = 0.5;
    const FloatType _cw1 = _cb1/(_kappa * _kappa) + (1.0 + _cb2) / _sigma;

    FloatType _laminarKinematicViscosity;
};