#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"

class TurbulenceModelBase {
public:
    explicit TurbulenceModelBase(const Config &config, const FluidBase &fluid, const Mesh &mesh);
    
    virtual ~TurbulenceModelBase() = default;

protected:
    const Config& _config;
    const FluidBase& _fluid;
    const Mesh& _mesh;
    size_t _ni, _nj, _nk;

    Matrix3D<FloatType> _eddyViscosity;
};
