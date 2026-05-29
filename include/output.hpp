#pragma once
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "output.hpp"
#include "math_utils.hpp"
#include "boundary_base.hpp"
#include "turbulence_model_base.hpp"
#include "turbulence_model_sa.hpp"
#include <fstream>
#include <filesystem>


class Output {
    
public:
    Output(
        const Config &config, 
        const Mesh &mesh, 
        const FlowSolution &solution, 
        const std::map<SolutionName, Matrix3D<Vector3D>>& solutionGrad,
        const FluidBase &fluid, 
        const std::vector<Boundary> &boundaries,
        const TurbulenceModelBase &turbulenceModel,
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

    void updatePrimaryFields();

    void updateSecondaryFields();

    void updateTurbulenceFields();

    Vector3D computeWallShearStress(Boundary boundary, size_t i, size_t j, size_t k);

protected:
    const Config& _config;
    const Mesh& _mesh;
    const FlowSolution& _solution;
    const std::map<SolutionName, Matrix3D<Vector3D>>& _solutionGrad;
    const FluidBase& _fluid;
    const std::vector<Boundary>& _boundaries;
    const TurbulenceModelBase& _turbulenceModel;
    const Matrix3D<Vector3D>& _surfacesI;
    const Matrix3D<Vector3D>& _surfacesJ;
    const Matrix3D<Vector3D>& _surfacesK;

    size_t _ni, _nj, _nk;
    std::map<std::string, Matrix3D<FloatType>> _outputFields;
    
    const Matrix3D<Vector3D>& _inviscidForce;
    const Matrix3D<Vector3D>& _viscousForce;
    const Matrix3D<FloatType>& _deviationAngle;
    
    bool _isUnsteadyOutput;
    bool _isBfmActive;
    bool _isTurbulenceActive;
    
    std::string _outputVolumeDirectory = "Volume_CSV";
    OutputFieldsType _outputFieldsType;
};