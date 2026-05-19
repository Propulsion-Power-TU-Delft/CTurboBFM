#pragma once

#include "boundary_base.hpp"

/**
 * @brief Class for fake boundaries whose flux is not needed
 */
class BoundaryFake : public BoundaryBase {
public:

    BoundaryFake(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid)
        : BoundaryBase(config, mesh, fluid) {}     

    virtual ~BoundaryFake()  = default;

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;

};