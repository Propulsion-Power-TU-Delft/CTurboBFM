#include "boundary_inlet_2d.hpp"
#include "math_utils.hpp"

BoundaryInlet2D::BoundaryInlet2D(
    const Config &config, 
    const Mesh &mesh, 
    const FluidBase &fluid, 
    std::string inletFile)
    : BoundaryBase(config, mesh, fluid) {
        
        _referenceFrame = config.getInletReferenceFrame();
        
        // read input file and check that is compatible with the mesh
        _inletFilePath = inletFile;
        _inputGrid = Input(_inletFilePath);
        _inputNi = _inputGrid.getNumberPointsI();
        _inputNj = _inputGrid.getNumberPointsJ();
        _inputNk = _inputGrid.getNumberPointsK();
        assert(_inputNi == 1 && "The inlet 2D file must have only one node along i.");
        assert((_inputNj == _mesh.getNumberPointsJ() && _inputNk == _mesh.getNumberPointsK()) &&
            "The inlet 2D file and the mesh must share the same number of nodes along j and k.");

    }


StateVector BoundaryInlet2D::computeBoundaryFlux(
            const StateVector& internalConservative, 
            const Vector3D& surface, 
            const Vector3D& midPoint, 
            const std::array<size_t, 3>& indices, 
            const FlowSolution& flowSolution, 
            const size_t& iterCounter) {
    
    FloatType nx = _inputGrid.getField(InputField::INLET_NX, indices[0], indices[1], indices[2]);
    FloatType ny = _inputGrid.getField(InputField::INLET_NY, indices[0], indices[1], indices[2]);
    FloatType nz = _inputGrid.getField(InputField::INLET_NZ, indices[0], indices[1], indices[2]);
    Vector3D flowDirection(nx, ny, nz);
    flowDirection /= flowDirection.magnitude();

    FloatType totalPressure = _inputGrid.getField(InputField::TOTAL_PRESSURE, indices[0], indices[1], indices[2]);
    FloatType totalTemperature = _inputGrid.getField(InputField::TOTAL_TEMPERATURE, indices[0], indices[1], indices[2]);
    
    StateVector flux = computeSubsonicInletFlux(
        internalConservative, 
        surface, 
        totalPressure, 
        totalTemperature, 
        flowDirection);
    
    return flux;
}