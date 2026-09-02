#pragma once

#include "types.hpp"
#include "config.hpp"
#include <assert.h>
#include "input.hpp"
#include <filesystem>

class Mesh {
    
public:
    Mesh(Config &config);
    
    ~Mesh() {};

    FloatType computeElementVolume(
        const std::vector<Vector3D> &boundSurfaces, 
        const std::vector<Vector3D> &boundCenters);

    const Matrix3D<Vector3D>& getSurfaces(FluxDirection direction) const;

    const Matrix3D<Vector3D>& getMidPoints(FluxDirection direction) const;

    const Matrix3D<Vector3D>& getSurfacesI() const {
        return _surfacesI;
    }

    FloatType getPeriodicityAngleRad() {
        return _periodicityAngleRad;
    }

    FloatType getPeriodicityAngleDeg() {
        return _periodicityAngleRad * 180.0 / M_PI;
    }
    
    const Matrix3D<Vector3D>& getMidPointsI() const {
        return _centersI;
    }

    const Matrix3D<Vector3D>& getSurfacesJ() const {
        return _surfacesJ;
    }
    
    const Matrix3D<Vector3D>& getMidPointsJ() const {
        return _centersJ;
    }

    const Matrix3D<Vector3D>& getSurfacesK() const {
        return _surfacesK;
    }

    const Matrix3D<Vector3D>& getMidPointsK() const {
        return _centersK;
    }
    
    const Matrix3D<Vector3D>& getVertices() const {
        return _vertices;
    }

    const Vector3D& getVertex(size_t i, size_t j, size_t k) const {
        return _vertices(i,j,k);
    }

    inline FloatType getRadius(size_t i, size_t j, size_t k) const {
        return std::sqrt(_vertices(i,j,k).y()*_vertices(i,j,k).y() + _vertices(i,j,k).z()*_vertices(i,j,k).z());
    }

    inline FloatType getTheta(size_t i, size_t j, size_t k) const {
        return std::atan2(_vertices(i,j,k).z(), _vertices(i,j,k).y());
    }
    
    const Matrix3D<Vector3D>& getDualNodes() const {
        return _dualNodes;
    }
    
    const Matrix3D<FloatType>& getVolumes() const {
        return _volumes;
    }

    const FloatType getVolume(size_t i, size_t j, size_t k) const {
        return _volumes(i,j,k);
    }

    const size_t getNumberDimensions() const {
        return _nDimensions;
    }

    const size_t getNumberPointsI() const {
        return _nPointsI;
    }
    
    const size_t getNumberPointsJ() const {
        return _nPointsJ;
    }
    
    const size_t getNumberPointsK() const {
        return _nPointsK;
    }

    void getElementEdges(size_t i, size_t j, size_t k, Vector3D &iEdge, Vector3D &jEdge, Vector3D &kEdge) const;

    FloatType getBoundaryTotalArea(BoundaryIndex boundIndex) const {
        return _boundaryAreas.at(boundIndex);
    }

    const Matrix2D<Vector3D> getMeshBoundary(BoundaryIndex boundIndex) const {
        return _boundarySurfaces.at(boundIndex); 
    }

    void computeSurfaceVectorAndCenter(
        const Vector3D &p1, 
        const Vector3D &p2, 
        const Vector3D &p3, 
        const Vector3D &p4, 
        Vector3D &normal, 
        Vector3D &center);

    FloatType getWedgeAngle() const {
        return _wedgeAngle;
    }

    FloatType getInputFields(InputField fieldName, size_t i, size_t j, size_t k) const {
        return _inputFile.getField(fieldName, i, j, k);
    }

    Matrix3D<FloatType> getInputFields(InputField fieldName) const {
        return _inputFile.getField(fieldName);
    }

    Vector3D getInputFieldsGradient(InputField fieldName, size_t i, size_t j, size_t k) const {
        return _gradientsMap.at(fieldName)(i, j, k);
    }

    void computeInputGradients();

    void checkPeriodicity(FloatType& translation, FloatType& angleRad);

    void setPeriodicMesh(FloatType angle, FloatType translation);

    void computeAdaptiveFlowDirection(Matrix3D<Vector3D> &flowDirection) const;

    void computeUniformFlowDirection(Vector3D initDirection, Matrix3D<Vector3D> &flowDirection) const;

    void writeMeshQualityStatistics() const;

    bool isPeriodicityActive() const {return _meshHasPeriodicity;}

private:

    void readPoints();

    void allocateMemory();

    void computeDualGrid2dTopology();
    
    void computeDualGrid3dTopology();

    void computeDualGrid1dTopology();

    void computeMeshInterfaces();

    void computeFiniteVolumes();
    
    void computeMeshQualityStatistics();

    void computeBoundaryAreas();
    
    void printMeshInfo();

    Matrix3D<FloatType> computeAspectRatio();

    void computeSkewnessAndOrthogonality(std::vector<FloatType> &skewness, std::vector<FloatType> &orthogonality);

private:
    const Config& _config;
    Input _inputFile;

    size_t _nPointsI, _nPointsJ, _nPointsK, _nPointsTotal;
    size_t _nDualPointsI, _nDualPointsJ, _nDualPointsK, _nDualPointsTotal;
    unsigned short int _nDimensions;
    
    Matrix3D<Vector3D> _surfacesI, _surfacesJ, _surfacesK; 
    Matrix3D<FloatType> _volumes;
    Matrix3D<Vector3D> _vertices;
    Matrix3D<Vector3D> _dualNodes;
    Matrix3D<Vector3D> _centersI, _centersJ, _centersK;
    std::map<BoundaryIndex, FloatType> _boundaryAreas;
    std::map<BoundaryIndex, Matrix2D<Vector3D>> _boundarySurfaces;
    
    FloatType _wedgeAngle {0.0}, _cellThickness {0.0}; 
    FloatType _periodicityAngleRad{0.0}, _periodicityTranslation{0.0};
    bool _meshHasPeriodicity {false};

    Matrix3D<FloatType> _aspectRatio;
    std::vector<FloatType> _skewness, _orthogonality;
    Statistics _aspectRatioStats, _skewnessStats, _orthogonalityStats;
    
    std::map<InputField, Matrix3D<Vector3D>> _gradientsMap;
};