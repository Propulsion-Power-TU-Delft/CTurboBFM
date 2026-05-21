#pragma once

#include "boundary_inlet.hpp"
#include "boundary_base.hpp"
#include "input.hpp"

/**
 * @brief Class for non-uniform inlets defined by 2D distributions.
 */
class BoundaryInlet2D : public BoundaryBase {
public:

    /**
     * @brief Constructs the object with given references.
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param boundIndex The boundary index.
     * @param inletFilePath The path to the inlet file.
     */
    BoundaryInlet2D(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        std::string inletFilePath);
        
    virtual ~BoundaryInlet2D() {}

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
    
protected:
        std::string _inletFilePath;
        ReferenceFrame _referenceFrame;
        Input _inputGrid;
        size_t _inputNi, _inputNj, _inputNk;
};