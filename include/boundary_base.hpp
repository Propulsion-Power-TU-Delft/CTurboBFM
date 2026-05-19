#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"


class BoundaryBase {
    
public:
    BoundaryBase(const Config &config, const Mesh &mesh, const FluidBase &fluid) 
        : _config(config), _mesh(mesh), _fluid(fluid) {};
        
    virtual ~BoundaryBase() = default;

    /**
     * @brief Compute the boundary flux.
     * @param internalConservative The internal point conservative variables.
     * @param surface The surface normal vector.
     * @param midPoint The midpoint of the boundary face.
     * @param indices The indices (i,j,k) of the boundary face.
     * @param flowSolution The entire flow solution.
     */
    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) = 0;

protected:

    /** Common interface for subsonic inlet boundary types. */
    StateVector computeSubsonicInletFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const FloatType& totPressureBoundary, 
        const FloatType& totTemperatureBoundary, 
        const Vector3D& flowDirection);
    
    /** Common interface for outlet boundary types. */
    StateVector computeOutletFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const FloatType& boundaryPressure,
        const size_t& iterCounter);

protected:

    const Config& _config;
    const Mesh& _mesh;
    const FluidBase& _fluid;
        
};



struct Boundary {
    std::string name;
    BoundaryType type;
    std::vector<FloatType> values;
    size_t i_min, i_max, j_min, j_max, k_min, k_max;
    std::shared_ptr<BoundaryBase> fluxMethod;
};