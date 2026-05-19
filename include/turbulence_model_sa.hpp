#pragma once

#include "turbulence_model_base.hpp"

class TurbulenceModelSA : public TurbulenceModelBase {
public:
    explicit TurbulenceModelSA() : TurbulenceModelBase() {};
    ~TurbulenceModelSA() override = default;

    // void solve() override;
    // void initialize() override;

private:
    // SA-specific fields
    double _nu_t;   // turbulent viscosity
    double _chi;    // viscosity ratio
};