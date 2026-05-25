#pragma once
#include "output_base.hpp"


class OutputCSV : public OutputBase {
    
public:
OutputCSV(
        const Config &config, 
        const Mesh &mesh, 
        const FlowSolution &solution, 
        const FluidBase &fluid, 
        const Matrix3D<Vector3D> &inviscidForce, 
        const Matrix3D<Vector3D> &viscousForce, 
        const Matrix3D<FloatType> &deviationAngle,
        const Matrix3D<FloatType> &wallDistance) 
        : OutputBase(config, mesh, solution, fluid, inviscidForce, viscousForce, deviationAngle, wallDistance) {};

virtual ~OutputCSV() = default;

virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) override;

};