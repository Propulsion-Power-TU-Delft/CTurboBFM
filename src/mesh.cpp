#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include "mesh.hpp"
#include "config.hpp"
#include "math_utils.hpp"

Mesh::Mesh(Config& config) : _config(config) {
    
    readPoints();
    allocateMemory();
    auto topology = config.getTopology();
    if (topology == Topology::TWO_DIMENSIONAL || topology == Topology::AXISYMMETRIC) {
        computeDualGrid2dTopology();
    } 
    else if (topology == Topology::THREE_DIMENSIONAL) {
        computeDualGrid3dTopology();
    }
    else if (topology == Topology::ONE_DIMENSIONAL) {
        computeDualGrid1dTopology();
    }
    else {
        throw std::runtime_error("Unsupported topology.");
    }
    
    computeMeshInterfaces();
    computeFiniteVolumes();
    computeMeshQualityStatistics();
    computeBoundaryAreas();
    computeInputGradients();

    printMeshInfo();

    if (_config.saveMeshQualityStats()) {
        writeMeshQualityStatistics();
    }
}


void Mesh::readPoints() {
    _inputFile = Input(_config.gridFilePath());

    _nPointsI = _inputFile.getNumberPointsI();
    _nPointsJ = _inputFile.getNumberPointsJ();
    _nPointsK = _inputFile.getNumberPointsK();
    _nPointsTotal = _nPointsI * _nPointsJ * _nPointsK;
    if (_nPointsJ==1 && _nPointsK==1) {
        _nDimensions = 1;
    }
    else if (_nPointsK==1) {
        _nDimensions = 2;
    }
    else {
        _nDimensions = 3;
    }
    
    _nDualPointsI = _nPointsI + 1;
    _nDualPointsJ = _nPointsJ + 1;
    _nDualPointsK = _nPointsK + 1;
    _nDualPointsTotal = _nDualPointsI * _nDualPointsJ * _nDualPointsK;

    _vertices.resize(_nPointsI, _nPointsJ, _nPointsK);
    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {
                _vertices(i,j,k) = _inputFile.getCoordinates(i,j,k);
            }
        }
    }
}


void Mesh::allocateMemory() {
    _dualNodes.resize(_nDualPointsI, _nDualPointsJ, _nDualPointsK);

    _surfacesI.resize(_nDualPointsI, _nDualPointsJ-1, _nDualPointsK-1);
    _surfacesJ.resize(_nDualPointsI-1, _nDualPointsJ, _nDualPointsK-1);
    _surfacesK.resize(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK);
    
    _centersI.resize(_nDualPointsI, _nDualPointsJ-1, _nDualPointsK-1);
    _centersJ.resize(_nDualPointsI-1, _nDualPointsJ, _nDualPointsK-1);
    _centersK.resize(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK);

    _volumes.resize(_nPointsI, _nPointsJ, _nPointsK);
}



void Mesh::printMeshInfo() {
    std::cout << std::endl;
    std::cout << "=========================================\n";
    std::cout << "        INFORMATION OF GRID FILE        \n";
    std::cout << "=========================================\n";
    std::cout << std::endl;

    std::cout << "Number of Dimensions    : " << _nDimensions << "\n";
    std::cout << "Points along Axis I     : " << _nPointsI << "\n";
    std::cout << "Points along Axis J     : " << _nPointsJ << "\n";
    std::cout << "Points along Axis K     : " << _nPointsK << "\n";
    std::cout << std::endl;

    std::cout << "The grid file '" << _config.gridFilePath() << "' contains a total of " << _nPointsTotal << " points.\n";
    std::cout << "Aspect ratio ";
    _aspectRatioStats.printInfo();
    std::cout << "Skewness ";
    _skewnessStats.printInfo();
    std::cout << "Orthogonality ";
    _orthogonalityStats.printInfo();
    std::cout << "=========================================\n";
}


void Mesh::computeMeshInterfaces() {
    std::cout << "Computing cell elements interfaces\n";

    // Computing the i-interfaces
    for (size_t i = 0; i < _nDualPointsI; i++) {
        for (size_t j = 0; j < _nDualPointsJ-1; j++) {
            for (size_t k = 0; k < _nDualPointsK-1; k++) {
                computeSurfaceVectorAndCenter(
                    _dualNodes(i,j,k), 
                    _dualNodes(i,j+1,k), 
                    _dualNodes(i,j+1,k+1), 
                    _dualNodes(i,j,k+1), 
                    _surfacesI(i,j,k), 
                    _centersI(i,j,k));
            }
        }
    }

    // Computing the j-interfaces
    for (size_t i = 0; i < _nDualPointsI-1; i++) {
        for (size_t j = 0; j < _nDualPointsJ; j++) {
            for (size_t k = 0; k < _nDualPointsK-1; k++) {
                computeSurfaceVectorAndCenter(
                    _dualNodes(i,j,k), 
                    _dualNodes(i,j,k+1), 
                    _dualNodes(i+1,j,k+1), 
                    _dualNodes(i+1,j,k), 
                    _surfacesJ(i,j,k), 
                    _centersJ(i,j,k));
            }
        }
    }

    // Computing the k-interfaces
    for (size_t i = 0; i < _nDualPointsI-1; i++) {
        for (size_t j = 0; j < _nDualPointsJ-1; j++) {
            for (size_t k = 0; k < _nDualPointsK; k++) {
                computeSurfaceVectorAndCenter(
                    _dualNodes(i,j,k), 
                    _dualNodes(i+1,j,k), 
                    _dualNodes(i+1,j+1,k), 
                    _dualNodes(i,j+1,k), 
                    _surfacesK(i,j,k), 
                    _centersK(i,j,k));
            }
        }
    }
}

void Mesh::computeSurfaceVectorAndCenter(
    const Vector3D &p1, 
    const Vector3D &p2, 
    const Vector3D &p3, 
    const Vector3D &p4, 
    Vector3D &normal, 
    Vector3D &center){

    // sides used to compute areas
    Vector3D a1, b1, a2, b2;
    a1 = p2 - p1;
    b1 = p3 - p2;
    a2 = p3 - p1;
    b2 = p4 - p3;

    // compute areas dividing the rectangle in two triangles
    Vector3D surf1, surf2;
    surf1 = a1.cross(b1) / 2.0;
    surf2 = a2.cross(b2) / 2.0;
    normal = surf1 + surf2;

    Vector3D center1, center2;
    center1 = (p1 + p2 + p3) / 3.0;
    center2 = (p1 + p3 + p4) / 3.0;
    center = (center1 * surf1.magnitude() + center2 * surf2.magnitude()) / (surf1.magnitude() + surf2.magnitude());
}

FloatType Mesh::computeElementVolume(
    const std::vector<Vector3D> &boundSurfaces, 
    const std::vector<Vector3D> &boundCenters) {
        
    FloatType volume {0.0};
    for (size_t i = 0; i < boundSurfaces.size(); ++i) {
        volume += boundSurfaces[i].x() * boundCenters[i].x();
    }
    return volume;
}

void Mesh::computeFiniteVolumes() {
    std::cout << "Compute element volumes\n";

    for (size_t i = 0; i< _nPointsI; i++) {
        for (size_t j = 0; j<_nPointsJ; j++) {
            for (size_t k = 0; k<_nPointsK; k++) {
                std::vector<Vector3D> boundSurfaces (6);
                std::vector<Vector3D> boundCenters (6);
        
                boundSurfaces[0] = - _surfacesI(i,j,k);
                boundSurfaces[1] = - _surfacesJ(i,j,k);
                boundSurfaces[2] = - _surfacesK(i,j,k);
                boundSurfaces[3] = _surfacesI(i+1,j,k);
                boundSurfaces[4] = _surfacesJ(i,j+1,k);
                boundSurfaces[5] = _surfacesK(i,j,k+1);

                boundCenters[0] = _centersI(i,j,k);
                boundCenters[1] = _centersJ(i,j,k);
                boundCenters[2] = _centersK(i,j,k); 
                boundCenters[3] = _centersI(i+1,j,k);
                boundCenters[4] = _centersJ(i,j+1,k);
                boundCenters[5] = _centersK(i,j,k+1);

                _volumes(i,j,k) = computeElementVolume(boundSurfaces, boundCenters);
            }
        }
    }
}

void Mesh::computeMeshQualityStatistics() {
    _aspectRatio = computeAspectRatio();
    _aspectRatioStats = Statistics(_aspectRatio);

    computeSkewnessAndOrthogonality(_skewness, _orthogonality);
    _skewnessStats = Statistics(_skewness);
    _orthogonalityStats = Statistics(_orthogonality);
}

void Mesh::computeSkewnessAndOrthogonality(std::vector<FloatType> &skewness, std::vector<FloatType> &orthogonality){

    for (size_t i=1; i<_nPointsI-1; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _vertices(i-1,j,k);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersI(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesI(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }

    // internal faces along j
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=1; j<_nPointsJ-1; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _vertices(i,j-1,k);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersJ(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesJ(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }

    // internal faces along k
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=1; k<_nPointsK-1; k++){
                Vector3D point0 = _vertices(i,j,k-1);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersK(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesK(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }
}


Matrix3D<FloatType> Mesh::computeAspectRatio() {
    Matrix3D<FloatType> aspectRatio(_nPointsI, _nPointsJ, _nPointsK);
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _dualNodes(i,j,k);
                Vector3D pointI = _dualNodes(i+1,j,k);
                Vector3D pointJ = _dualNodes(i,j+1,k);
                Vector3D pointK = _dualNodes(i,j,k+1);
                FloatType iEdge = (pointI - point0).magnitude();
                FloatType jEdge = (pointJ - point0).magnitude();
                FloatType kEdge = (pointK - point0).magnitude();
                if (_nDimensions == 2) {
                    FloatType maxEdge = std::max(iEdge, jEdge);
                    FloatType minEdge = std::min(iEdge, jEdge);
                    aspectRatio(i, j, k) = maxEdge / minEdge;
                } else {
                    FloatType maxEdge = std::max({iEdge, jEdge, kEdge});
                    FloatType minEdge = std::min({iEdge, jEdge, kEdge});
                    aspectRatio(i, j, k) = maxEdge / minEdge;
                }
            }
        }
    }
    return aspectRatio;
}


void Mesh::computeBoundaryAreas() {

    _boundarySurfaces[BoundaryIndex::I_START] = _surfacesI.getBoundarySlice(BoundaryIndex::I_START);
    _boundarySurfaces[BoundaryIndex::I_END] = _surfacesI.getBoundarySlice(BoundaryIndex::I_END);
    _boundarySurfaces[BoundaryIndex::J_START] = _surfacesJ.getBoundarySlice(BoundaryIndex::J_START);
    _boundarySurfaces[BoundaryIndex::J_END] = _surfacesJ.getBoundarySlice(BoundaryIndex::J_END);
    _boundarySurfaces[BoundaryIndex::K_START] = _surfacesK.getBoundarySlice(BoundaryIndex::K_START);
    _boundarySurfaces[BoundaryIndex::K_END] = _surfacesK.getBoundarySlice(BoundaryIndex::K_END);

    std::array<BoundaryIndex, 6> BoundaryIndexArray = {
        BoundaryIndex::I_START,
        BoundaryIndex::I_END,
        BoundaryIndex::J_START,
        BoundaryIndex::J_END,
        BoundaryIndex::K_START,
        BoundaryIndex::K_END
    };
    for (const auto& index : BoundaryIndexArray) {
        _boundaryAreas[index] = computeTotalSurfaceArea(_boundarySurfaces[index]);
    }

}

const Matrix3D<Vector3D>& Mesh::getSurfaces(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return _surfacesI;
        case FluxDirection::J: return _surfacesJ;
        case FluxDirection::K: return _surfacesK;
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}

const Matrix3D<Vector3D>& Mesh::getMidPoints(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return _centersI;
        case FluxDirection::J: return _centersJ;
        case FluxDirection::K: return _centersK;
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}


void Mesh::computeDualGrid2dTopology() {
    assert (_nDimensions == 2 && "Can only compute a 2D dual grid for a 2D mesh.");
    std::cout << "Computing dual grid coordinates for 2D mesh\n";
    
    Matrix2D<Vector3D> nodes(_nDualPointsI, _nDualPointsJ);

    // dual nodes internal to the volume
    for (size_t i = 1; i < _nDualPointsI-1; i++) {
        for (size_t j = 1; j < _nDualPointsJ-1; j++) {
            nodes(i,j) = (_vertices(i-1,j-1,0) + _vertices(i,j-1,0) + _vertices(i,j,0) + _vertices(i-1,j,0)) / 4.0;
        }
    }
    
    // on the corners
    nodes(0,0) = _vertices(0,0,0);
    nodes(0, _nDualPointsJ-1) = _vertices(0, _nPointsJ-1, 0);
    nodes(_nDualPointsI-1, 0) = _vertices(_nPointsI-1, 0, 0);
    nodes(_nDualPointsI-1, _nDualPointsJ-1) = _vertices(_nPointsI-1, _nPointsJ-1, 0);

    // on the edges
    for (size_t i = 1; i < _nDualPointsI-1; i++) {
        nodes(i, 0) = (_vertices(i-1,0,0) + _vertices(i,0,0)) / 2.0;
        nodes(i, _nDualPointsJ-1) = (_vertices(i-1, _nPointsJ-1, 0) + _vertices(i, _nPointsJ-1, 0)) / 2.0;
    }
    for (size_t j = 1; j < _nDualPointsJ-1; j++) {
        nodes(0, j) = (_vertices(0, j-1, 0) + _vertices(0, j, 0)) / 2.0;    
        nodes(_nDualPointsI-1, j) = (_vertices(_nPointsI-1, j-1, 0) + _vertices(_nPointsI-1, j, 0)) / 2.0;    
    }

    if (_config.getTopology() == Topology::AXISYMMETRIC) { 
        // build a wedge of 1 degree in the circumferential (k) direction
        _wedgeAngle = 1.0 * M_PI / 180.0;
        for (size_t i = 0; i < _nDualPointsI; i++) {
            for (size_t j = 0; j < _nDualPointsJ; j++) {
                _dualNodes(i,j,0).x() = nodes(i,j).x();
                _dualNodes(i,j,1).x() = nodes(i,j).x();

                _dualNodes(i,j,0).y() = nodes(i,j).y() * cos(-_wedgeAngle/2.0);
                _dualNodes(i,j,1).y() = nodes(i,j).y() * cos(+_wedgeAngle/2.0);

                _dualNodes(i,j,0).z() = nodes(i,j).y() * sin(-_wedgeAngle/2.0);
                _dualNodes(i,j,1).z() = nodes(i,j).y() * sin(+_wedgeAngle/2.0);
            }
        }
    }
    else {
        // for 2D use a unitary thickness in the third direction
        _cellThickness = 1.0;
        for (size_t i = 0; i < _nDualPointsI; i++) {
            for (size_t j = 0; j < _nDualPointsJ; j++) {
                _dualNodes(i,j,0).x() = nodes(i,j).x();
                _dualNodes(i,j,1).x() = nodes(i,j).x();

                _dualNodes(i,j,0).y() = nodes(i,j).y();
                _dualNodes(i,j,1).y() = nodes(i,j).y();

                _dualNodes(i,j,0).z() = -_cellThickness / 2.0;
                _dualNodes(i,j,1).z() = +_cellThickness / 2.0;
            }
        }
            
    }
}


void Mesh::computeDualGrid1dTopology() {
    assert (_nDimensions == 1 && "Can only compute a 1D dual grid for a 1D mesh.");
    std::cout << "Computing dual grid coordinates for 1D mesh\n";

    Matrix2D<Vector3D> nodes(_nDualPointsI, 1); // use the same structure as 2D but with only one point in j direction

    // internal dual nodes
    for (size_t i = 1; i < _nDualPointsI-1; i++) {
        nodes(i,0) = (_vertices(i-1,0,0) + _vertices(i,0,0)) / 2.0;
    }
    
    // extremes
    nodes(0,0) = _vertices(0,0,0);
    nodes(_nDualPointsI-1, 0) = _vertices(_nPointsI-1, 0, 0);

    // build elements
    _cellThickness = 1.0;
    for (size_t i = 0; i < _nDualPointsI; i++) {
        // x-coords
        _dualNodes(i,0,0).x() = nodes(i,0).x();
        _dualNodes(i,1,0).x() = nodes(i,0).x();
        _dualNodes(i,0,1).x() = nodes(i,0).x();
        _dualNodes(i,1,1).x() = nodes(i,0).x();

        // y-coords
        _dualNodes(i,0,0).y() = -_cellThickness / 2.0;
        _dualNodes(i,1,0).y() = +_cellThickness / 2.0;
        _dualNodes(i,0,1).y() = -_cellThickness / 2.0;
        _dualNodes(i,1,1).y() = +_cellThickness / 2.0;

        // z-coords
        _dualNodes(i,0,0).z() = -_cellThickness / 2.0;
        _dualNodes(i,1,0).z() = -_cellThickness / 2.0;
        _dualNodes(i,0,1).z() = +_cellThickness / 2.0;
        _dualNodes(i,1,1).z() = +_cellThickness / 2.0;
    }
            
}


void Mesh::computeDualGrid3dTopology() {
    assert (_nDimensions == 3 && "Can only compute a 3D dual grid for a 3D mesh.");
    std::cout << "Computing dual grid coordinates for 3D mesh\n";

    // internal dual nodes
    for (size_t i = 1; i < _nDualPointsI-1; i++) {
        for (size_t j = 1; j < _nDualPointsJ-1; j++) {
            for (size_t k = 1; k < _nDualPointsK-1; k++){
                _dualNodes(i,j,k) = (
                    _vertices(i-1,j-1,k-1) + _vertices(i,j-1,k-1) + _vertices(i-1,j,k-1) + _vertices(i-1,j-1,k) + \
                    _vertices(i  ,j,  k-1) + _vertices(i,j-1,k)   + _vertices(i-1,j,k)   + _vertices(i  ,j  ,k) ) / 8.0;
            }
        }
    }
    
    // corners
    _dualNodes(0,0,0) = _vertices(0,0,0);
    _dualNodes(_nDualPointsI-1, 0, 0) = _vertices(_nPointsI-1, 0, 0);
    _dualNodes(0, _nDualPointsJ-1, 0) = _vertices(0, _nPointsJ-1, 0);
    _dualNodes(0, 0, _nDualPointsK-1) = _vertices(0, 0, _nPointsK-1);
    _dualNodes(_nDualPointsI-1, _nDualPointsJ-1, 0) = _vertices(_nPointsI-1, _nPointsJ-1, 0);
    _dualNodes(_nDualPointsI-1, 0, _nDualPointsK-1) = _vertices(_nPointsI-1, 0, _nPointsK-1);
    _dualNodes(0, _nDualPointsJ-1, _nDualPointsK-1) = _vertices(0, _nPointsJ-1, _nPointsK-1);
    _dualNodes(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK-1) = _vertices(_nPointsI-1, _nPointsJ-1, _nPointsK-1);

    // edges
    // i-oriented
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        _dualNodes(i,0,0) = (_vertices(i-1,0,0) + _vertices(i,0,0)) / 2.0;
        _dualNodes(i,_nDualPointsJ-1,0) = (_vertices(i-1,_nPointsJ-1,0) + _vertices(i,_nPointsJ-1,0)) / 2.0;
        _dualNodes(i,0,_nDualPointsK-1) = (_vertices(i-1,0,_nPointsK-1) + _vertices(i,0,_nPointsK-1)) / 2.0;
        _dualNodes(i,_nDualPointsJ-1,_nDualPointsK-1) = (
            _vertices(i-1,_nPointsJ-1,_nPointsK-1) + _vertices(i,_nPointsJ-1,_nPointsK-1)) / 2.0;
    }
    // j-oriented
    for (size_t j=1; j<_nDualPointsJ-1; j++) {
        _dualNodes(0,j,0) = (_vertices(0,j-1,0) + _vertices(0,j,0)) / 2.0;
        _dualNodes(_nDualPointsI-1,j,0) = (_vertices(_nPointsI-1,j-1,0) + _vertices(_nPointsI-1,j,0)) / 2.0;
        _dualNodes(0,j,_nDualPointsK-1) = (_vertices(0,j-1,_nPointsK-1) + _vertices(0,j,_nPointsK-1)) / 2.0;
        _dualNodes(_nDualPointsI-1,j,_nDualPointsK-1) = (
            _vertices(_nPointsI-1,j-1,_nPointsK-1) + _vertices(_nPointsI-1,j,_nPointsK-1)) / 2.0;
    }
    // k-oriented
    for (size_t k=1; k<_nDualPointsK-1; k++) {
        _dualNodes(0,0,k) = (_vertices(0,0,k-1) + _vertices(0,0,k)) / 2.0;
        _dualNodes(_nDualPointsI-1,0,k) = (_vertices(_nPointsI-1,0,k-1) + _vertices(_nPointsI-1,0,k)) / 2.0;
        _dualNodes(0,_nDualPointsJ-1,k) = (_vertices(0,_nPointsJ-1,k-1) + _vertices(0,_nPointsJ-1,k)) / 2.0;
        _dualNodes(_nDualPointsI-1,_nDualPointsJ-1,k) = (
            _vertices(_nPointsI-1,_nPointsJ-1,k-1) + _vertices(_nPointsI-1,_nPointsJ-1,k)) / 2.0;
    }
    
    // external boundaries
    // i-faces
    for (size_t j=1; j<_nDualPointsJ-1; j++) {
        for (size_t k=1; k<_nDualPointsK-1; k++) {
            _dualNodes(0,j,k) = (_vertices(0,j-1,k-1) + _vertices(0,j,k-1) + _vertices(0,j-1,k) + _vertices(0,j,k)) / 4.0;
            _dualNodes(_nDualPointsI-1,j,k) = (
                _vertices(_nPointsI-1,j-1,k-1) + _vertices(_nPointsI-1,j,k-1) + 
                _vertices(_nPointsI-1,j-1,k) + _vertices(_nPointsI-1,j,k)) / 4.0;
        }
    }
    // j-faces
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        for (size_t k=1; k<_nDualPointsK-1; k++) {
            _dualNodes(i,0,k) = (_vertices(i-1,0,k-1) + _vertices(i,0,k-1) + _vertices(i-1,0,k) + _vertices(i,0,k)) / 4.0;
            _dualNodes(i,_nDualPointsJ-1,k) = (
                _vertices(i-1,_nPointsJ-1,k-1) + _vertices(i,_nPointsJ-1,k-1) +
                 _vertices(i-1,_nPointsJ-1,k) + _vertices(i,_nPointsJ-1,k)) / 4.0;
        }
    }
    // k-faces
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        for (size_t j=1; j<_nDualPointsJ-1; j++) {
            _dualNodes(i,j,0) = (_vertices(i-1,j-1,0) + _vertices(i,j-1,0) + _vertices(i-1,j,0) + _vertices(i,j,0)) / 4.0;
            _dualNodes(i,j,_nDualPointsK-1) = (
                _vertices(i-1,j-1,_nPointsK-1) + _vertices(i,j-1,_nPointsK-1) +
                 _vertices(i-1,j,_nPointsK-1) + _vertices(i,j,_nPointsK-1)) / 4.0;
        }
    }
}

void Mesh::getElementEdges(size_t i, size_t j, size_t k, Vector3D &iEdge, Vector3D &jEdge, Vector3D &kEdge) const {
    Vector3D pt0, ptI, ptJ, ptK;
    pt0 = _dualNodes(i, j, k);
    ptI = _dualNodes(i+1, j, k);
    ptJ = _dualNodes(i, j+1, k);
    ptK = _dualNodes(i, j, k+1);
    iEdge = ptI - pt0;
    jEdge = ptJ - pt0;
    kEdge = ptK - pt0;
}


void Mesh::computeInputGradients() {
    if (_config.isBFMActive()){
        _gradientsMap[InputField::BLOCKAGE] = Matrix3D<Vector3D>(_nPointsI, _nPointsJ, _nPointsK);
        auto blockage = getInputFields(InputField::BLOCKAGE);
        computeGradientGreenGauss(
            _surfacesI, _surfacesJ, _surfacesK, 
            _centersI, _centersJ, _centersK, 
            _vertices, _volumes, blockage, 
            _gradientsMap[InputField::BLOCKAGE]);
    }
}


void Mesh::checkPeriodicity(FloatType& translation, FloatType& angle) {
    const FloatType tolerance = 1E-4;  
    const Vector3D& drag{0.0, 0.0, translation};
    
    for (size_t i = 0; i < _nPointsI; ++i) {
        for (size_t j = 0; j < _nPointsJ; ++j) {
            const Vector3D& pointA = _vertices(i, j, 0);
            const Vector3D& pointB = _vertices(i, j, _nPointsK - 1) - drag;

            const FloatType thetaA = atan2FromZeroTo2pi(pointA.z(), pointA.y());
            const FloatType thetaB = atan2FromZeroTo2pi(pointB.z(), pointB.y());
            FloatType deltaTheta = std::fmod(thetaB - thetaA, 2.0 * M_PI);

            if (std::abs(deltaTheta - angle) > tolerance) {
                std::ostringstream msg;
                msg << "Periodicity check failed at (i,j)=("
                    << i << "," << j << "): "
                    << "Δθ=" << deltaTheta << ", expected " << angle;
                throw std::invalid_argument(msg.str());
            }

            if (std::abs(translation) > 1e-6 && angle==0.0) {
                const FloatType dz = pointB.z() - pointA.z();
                if (std::abs(dz) > tolerance) {
                    std::ostringstream msg;
                    msg << "Axial translation mismatch at (i,j)=("
                        << i << "," << j << "): "
                        << "Δz=" << dz << ", expected " << 0.0;
                    throw std::invalid_argument(msg.str());
                }
            }
        }
    }
    
    setPeriodicMesh(angle, translation);
}



void Mesh::computeUniformFlowDirection(Vector3D initDirection, Matrix3D<Vector3D> &flowDirection) const{
    std::cout << "Uniform flow direction initialization" << std::endl;
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                flowDirection(i,j,k) = initDirection;
            }
        }
    }
}


void Mesh::computeAdaptiveFlowDirection(Matrix3D<Vector3D> &flowDirection) const{
    std::cout << "Adaptive flow direction initialization" << std::endl;

    // up to last-1 point in the streamwise direction
    for (size_t i=0; i<_nPointsI-1; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D pointA = _vertices(i,j,k);
                Vector3D pointB = _vertices(i+1,j,k);
                flowDirection(i,j,k) = (pointB - pointA).normalized();
            }
        }
    }

    // last streamwise position
    for (size_t j=0; j<_nPointsJ; j++){
        for (size_t k=0; k<_nPointsK; k++){
            Vector3D pointA = _vertices(_nPointsI-2,j,k);
            Vector3D pointB = _vertices(_nPointsI-1,j,k);
            flowDirection(_nPointsI-1,j,k) = (pointB - pointA).normalized();
        }
    }
}


void Mesh::writeMeshQualityStatistics() const {
    std::string _meshQualityFolderName = "Mesh_Quality";
    try {
        if (!std::filesystem::exists(_meshQualityFolderName))
            std::filesystem::create_directories(_meshQualityFolderName);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
        return;
    }
    writeDataToCsv(_skewness, _meshQualityFolderName + "/skewness.csv");
    writeDataToCsv(_orthogonality, _meshQualityFolderName + "/orthogonality.csv");
    writeDataToCsv(_aspectRatio.getData(), _meshQualityFolderName + "/aspect_ratio.csv");
}

void Mesh::setPeriodicMesh(FloatType periodicityAngle, FloatType periodicityTranslation) {
    _periodicityAngleRad = periodicityAngle;
    _periodicityTranslation = periodicityTranslation;
    _meshHasPeriodicity = true;
}