#pragma once

#include "types.hpp"
#include "fluid_base.hpp"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

FloatType computeTotalSurfaceArea(const Matrix2D<Vector3D> &surfaces);

FloatType computeSurfaceIntegral(
    const Matrix2D<Vector3D> &surfaces, 
    const Matrix2D<FloatType> &vecX, 
    const Matrix2D<FloatType> &vecY, 
    const Matrix2D<FloatType> &vecZ);

StateVector getPrimitiveVariablesFromConservative(StateVector conservative);

StateVector getConservativeVariablesFromPrimitive(StateVector primitive);

StateVector computeAdvectionFluxFromPrimitive(StateVector primitive, Vector3D surface, const FluidBase& fluid);

StateVector computeAdvectionFluxFromConservative(StateVector conservative, Vector3D surface, const FluidBase& fluid);


/** formulation taken from Hirsch second book, as summation of A*nx + B*ny + C*nz */
Matrix2D<FloatType> computeAdvectionJacobian(
    const StateVector& primitive, 
    const Vector3D& direction, 
    const FluidBase& fluid);

FloatType computeAngleBetweenVectors(const Vector3D& v1, const Vector3D& v2);

Vector3D rotateVectorAlongXAxis(const Vector3D& v, FloatType theta);

StateVector rotateStateVectorAlongXAxis(const StateVector& v, FloatType theta);

/** Compute cylindrical components (axial, radial, tangential) from cartesian components (x, y, z) */
Vector3D computeCylindricalComponentsFromCartesian(Vector3D vec, FloatType theta);

/** Compute cartesian components (x, y, z) from cylindrical components (axial, radial, tangential) */
Vector3D computeCartesianComponentsFromCylindrical(Vector3D vec, FloatType theta);

void integrateRadialEquilibrium(
    const std::vector<FloatType>& density, 
    const std::vector<FloatType>& velocityTang, 
    const std::vector<FloatType>& radius, 
    std::vector<FloatType>& pressure,
    FloatType& newHubPressure);

void computeGradientGreenGauss(
    const Matrix3D<Vector3D>& surfacesI, 
    const Matrix3D<Vector3D>& surfacesJ, 
    const Matrix3D<Vector3D>& surfacesK, 
    const Matrix3D<Vector3D>& midPointsI, 
    const Matrix3D<Vector3D>& midPointsJ, 
    const Matrix3D<Vector3D>& midPointsK, 
    const Matrix3D<Vector3D>& nodes, 
    const Matrix3D<FloatType>& volumes, 
    const Matrix3D<FloatType>& field, 
    Matrix3D<Vector3D>& gradient);

void computeGradientsGreenGauss4(
    const Matrix3D<Vector3D>& surfacesI, 
    const Matrix3D<Vector3D>& surfacesJ, 
    const Matrix3D<Vector3D>& surfacesK, 
    const Matrix3D<Vector3D>& midPointsI, 
    const Matrix3D<Vector3D>& midPointsJ, 
    const Matrix3D<Vector3D>& midPointsK, 
    const Matrix3D<Vector3D>& nodes, 
    const Matrix3D<FloatType>& volumes, 
    const Matrix3D<FloatType>& field0, 
    const Matrix3D<FloatType>& field1, 
    const Matrix3D<FloatType>& field2, 
    const Matrix3D<FloatType>& field3, 
    Matrix3D<Vector3D>& grad0,
    Matrix3D<Vector3D>& grad1,
    Matrix3D<Vector3D>& grad2,
    Matrix3D<Vector3D>& grad3);

FloatType interpolateScalar(
    const Matrix3D<FloatType>& field, 
    size_t i, 
    size_t j, 
    size_t k, 
    const Matrix3D<Vector3D>& nodes, 
    const Vector3D& point, 
    Direction3D direction);

static inline Vector3D computeGreenGaussDivergence(
    const std::array<Vector3D,6>& surfaces, 
    const std::array<FloatType,6>& scalars, 
    const FloatType& volume){

    Vector3D gradient(0.0,0.0,0.0);
    int length = surfaces.size();
    for (int i=0; i<length; i++){
        gradient += surfaces[i] * scalars[i];
    }
    return gradient / volume;
}

FloatType atan2FromZeroTo2pi(FloatType y, FloatType x);

FloatType inverseRescalingMinMax(const FloatType& value, const FloatType& min, const FloatType& max);

FloatType rescaleMinMax(const FloatType& value, const FloatType& min, const FloatType& max);

FloatType linearInterpolation(const std::vector<double>& x, const std::vector<double>& y, const FloatType& xp);

void writeDataToCsv(const std::vector<FloatType>& data, const std::string& filename);

struct ViscousStressTensor {
    FloatType xx = 0.0;
    FloatType yy = 0.0;
    FloatType zz = 0.0;
    FloatType xy = 0.0;
    FloatType xz = 0.0;
    FloatType yz = 0.0;
};

ViscousStressTensor computeViscousStressTensor(
    const FloatType& mu, 
    const FloatType& lambda, 
    const Vector3D& gradU, 
    const Vector3D& gradV, 
    const Vector3D& gradW
);

Vector3D computeViscousStressVector(const ViscousStressTensor& viscousStressTensor, const Vector3D& surface);

FloatType computeRotationRateMagnitude(
    const Vector3D& gradU, 
    const Vector3D& gradV, 
    const Vector3D& gradW
);