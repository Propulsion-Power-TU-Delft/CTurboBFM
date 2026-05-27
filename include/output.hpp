#pragma once
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "output.hpp"
#include "math_utils.hpp"
#include <fstream>
#include <filesystem>


class Output {
    
public:
    Output(
        const Config &config, 
        const Mesh &mesh, 
        const FlowSolution &solution, 
        const FluidBase &fluid, 
        const Matrix3D<Vector3D> &inviscidForce, 
        const Matrix3D<Vector3D> &viscousForce, 
        const Matrix3D<FloatType> &deviationAngle);
    
    ~Output() = default;

    void writeSolution(size_t iterationCounter);

    void getOutputFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap) const;

    std::string getOutputFilename(size_t iterationCounter);

protected:
    void allocateOutputArrays();
    
    void updateOutputArrays();

    void writeHeader(std::ofstream& file) const;

    void writeData(std::ofstream& file) const;

protected:
    const Config& _config;
    const Mesh& _mesh;
    const FlowSolution& _solution;
    const FluidBase& _fluid;
    size_t _ni, _nj, _nk;
    std::map<std::string, Matrix3D<FloatType>> _outputFields;
    
    const Matrix3D<Vector3D>& _inviscidForce;
    const Matrix3D<Vector3D>& _viscousForce;
    const Matrix3D<FloatType>& _deviationAngle;
    
    bool _isUnsteadyOutput;
    bool _isBfmActive;
    
    std::string _outputVolumeDirectory = "Volume_CSV";
    OutputFieldsType _outputFieldsType;
};