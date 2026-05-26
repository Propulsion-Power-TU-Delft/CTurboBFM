#pragma once
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "output_base.hpp"
#include <fstream>
#include <filesystem>


class OutputBase {
    
public:
    OutputBase(
        const Config &config, 
        const Mesh &mesh, 
        const FlowSolution &solution, 
        const FluidBase &fluid, 
        const Matrix3D<Vector3D> &inviscidForce, 
        const Matrix3D<Vector3D> &viscousForce, 
        const Matrix3D<FloatType> &deviationAngle,
        const Matrix3D<FloatType> &wallDistance);
    
    virtual ~OutputBase() = default;

    virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) = 0;

    void getOutputFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap, bool alsoGradients=false) const;

    std::string getOutputFilename(size_t iterationCounter);

protected:
    void allocateSpaceForOutput(
        std::map<std::string, Matrix3D<FloatType>>& fieldsMap, 
        bool alsoGradients) const;
    
    void storeFields(
        std::map<std::string, Matrix3D<FloatType>>& fieldsMap, 
        bool alsoGradients) const;

protected:
    const Config& _config;
    const Mesh& _mesh;
    const FlowSolution& _solution;
    const FluidBase& _fluid;
    
    const Matrix3D<Vector3D>& _inviscidForce;
    const Matrix3D<Vector3D>& _viscousForce;
    const Matrix3D<FloatType>& _deviationAngle;
    const Matrix3D<FloatType>& _wallDistance;
    
    bool _isUnsteadyOutput;
    
    std::string _outputDirectory = "Volume_CSV";
    OutputFields _outputFields;
};