#pragma once
#include "boundary_base.hpp"

/**
 * @brief Characteristic-based farfield boundary condition.
 *
 * Non-reflecting external-flow boundary based on Riemann invariants. The
 * outgoing characteristic is taken from the interior state, the incoming one
 * from the prescribed freestream state, and the two are combined to build the
 * boundary state. Subsonic/supersonic inflow and outflow are handled
 * automatically depending on the local normal Mach number.
 *
 * Freestream values are provided as [Mach, p, T, dx, dy, dz].
 */
class BoundaryFarfield : public BoundaryBase {

public:
    /**
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param farfieldValues The freestream values (Mach, static pressure,
     *        static temperature, flow direction x/y/z).
     */
    BoundaryFarfield(
        const Config& config,
        const Mesh& mesh,
        const FluidBase& fluid,
        std::vector<FloatType> farfieldValues)
        : BoundaryBase(config, mesh, fluid),
        _boundaryValues(farfieldValues) {}

    virtual ~BoundaryFarfield() = default;

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative,
        const Vector3D& surface,
        const Vector3D& midPoint,
        const std::array<size_t, 3>& indices,
        const FlowSolution& flowSolution,
        const size_t& iterCounter) override;

protected:
    std::vector<FloatType> _boundaryValues;
};
