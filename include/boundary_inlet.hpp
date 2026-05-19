#pragma once

#include "boundary_base.hpp"

/**
 * @brief Class for inlet boundary condition flux calculation.
 * Provides a common interface for it.
 */
class BoundaryInlet : public BoundaryBase {
    
public:
    /**
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param boundIndex The boundary index.
     * @param inletValues The inlet values (total pressure, total temperature, flow direction).
     */
    BoundaryInlet(
        const Config& config,
        const Mesh& mesh,
        const FluidBase& fluid,
        const std::vector<FloatType>& inletValues)
        : BoundaryBase(config, mesh, fluid),
        _boundaryValues(inletValues),
        _referenceFrame(config.getInletReferenceFrame()) {}
        
    virtual ~BoundaryInlet() = default;

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;

protected:
    std::vector<FloatType> _boundaryValues;
    ReferenceFrame _referenceFrame;
};