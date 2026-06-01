#include "turbulence_model_base.hpp"

TurbulenceModelBase::TurbulenceModelBase(
    const Config &config, 
    const FluidBase &fluid, 
    const Mesh &mesh, 
    const std::vector<Boundary> &boundaries,
    const FlowSolution &initialSolution
) 
    : _config(config), 
    _fluid(fluid), 
    _mesh(mesh), 
    _boundaries(boundaries),
    _initialSolution(initialSolution) {
        
        _ni = _mesh.getNumberPointsI();
        _nj = _mesh.getNumberPointsJ();
        _nk = _mesh.getNumberPointsK();
        _topology = _config.getTopology();
    }

const std::array<int, 3> TurbulenceModelBase::getStepMask(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return {1, 0, 0};
        case FluxDirection::J: return {0, 1, 0};
        case FluxDirection::K: return {0, 0, 1};
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}