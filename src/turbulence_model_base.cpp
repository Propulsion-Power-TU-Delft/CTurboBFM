#include "turbulence_model_base.hpp"

TurbulenceModelBase::TurbulenceModelBase(
    const Config &config, 
    const FluidBase &fluid, 
    const Mesh &mesh, 
    const std::vector<Boundary> &boundaries) 
    : _config(config), 
    _fluid(fluid), 
    _mesh(mesh), 
    _boundaries(boundaries) {
        
        _ni = _mesh.getNumberPointsI();
        _nj = _mesh.getNumberPointsJ();
        _nk = _mesh.getNumberPointsK();
    }