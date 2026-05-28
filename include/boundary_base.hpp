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

enum class BoundaryOrientation {
    I_START,
    I_END,
    J_START,
    J_END,
    K_START,
    K_END
};

struct Boundary {
    std::string name;
    BoundaryType type;
    std::vector<FloatType> values;
    size_t i_min, i_max, j_min, j_max, k_min, k_max;
    std::shared_ptr<BoundaryBase> fluxMethod;
    BoundaryOrientation orientation;

    void computeOrientation() {
        if (i_min == i_max) {
            orientation = i_min == 0 ? BoundaryOrientation::I_START : BoundaryOrientation::I_END;
        }
        else if (j_min == j_max) {
            orientation = j_min == 0 ? BoundaryOrientation::J_START : BoundaryOrientation::J_END;
        }
        else if (k_min == k_max) {
            orientation = k_min == 0 ? BoundaryOrientation::K_START : BoundaryOrientation::K_END;
        }
    }
};

struct RadialEquilibriumProfile {
    std::vector<FloatType> pressure;
    std::vector<FloatType> radius;
    Boundary boundary;
};

// indexing referred to primary nodes
struct BoundaryNodesIndexRange
{
    size_t iStart, iLast;
    size_t jStart, jLast;
    size_t kStart, kLast;
};

BoundaryNodesIndexRange fetchBoundaryNodesIndexRange(
    const Boundary& boundary,
    const size_t& nPointsI, 
    const size_t& nPointsJ, 
    const size_t& nPointsK);