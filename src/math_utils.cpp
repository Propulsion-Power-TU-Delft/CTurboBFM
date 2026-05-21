#include "math_utils.hpp"
#include <iostream>


FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    
    FloatType sum = 0.0;
    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).dot(vectors(i,j));
        }
    }
    return sum;
}

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, 
    const Matrix2D<FloatType> &vecX, 
    const Matrix2D<FloatType> &vecY, 
    const Matrix2D<FloatType> &vecZ){

    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).x()*vecX(i,j) + surfaces(i,j).y()*vecY(i,j) + surfaces(i,j).z()*vecZ(i,j);
        }
    }
    return sum;
}

FloatType computeTotalSurfaceArea(const Matrix2D<Vector3D> &surfaces){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).magnitude();
        }
    }
    return sum;
}


StateVector getPrimitiveVariablesFromConservative(StateVector conservative){
    StateVector primitive;
    primitive[0] = conservative[0];
    primitive[1] = conservative[1] / conservative[0];
    primitive[2] = conservative[2] / conservative[0];
    primitive[3] = conservative[3] / conservative[0];
    primitive[4] = conservative[4] / conservative[0];
    return primitive;
}


StateVector getConservativeVariablesFromPrimitive(StateVector primitive){
    StateVector conservative;
    conservative[0] = primitive[0];
    conservative[1] = primitive[1] * primitive[0];
    conservative[2] = primitive[2] * primitive[0];
    conservative[3] = primitive[3] * primitive[0];
    conservative[4] = primitive[4] * primitive[0];
    return conservative;
}


StateVector computeAdvectionFluxFromPrimitive(StateVector primitive, Vector3D surface, const FluidBase& fluid){
    Vector3D normal = surface / surface.magnitude();
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    FloatType normalVelocity = velocity.dot(normal);
    FloatType pressure = fluid.computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType totEnthalpy = fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);

    StateVector flux {};
    flux[0] = primitive[0] * normalVelocity;
    flux[1] = primitive[0] * normalVelocity * velocity.x() + pressure * normal.x();
    flux[2] = primitive[0] * normalVelocity * velocity.y() + pressure * normal.y();
    flux[3] = primitive[0] * normalVelocity * velocity.z() + pressure * normal.z();
    flux[4] = primitive[0] * normalVelocity * totEnthalpy;
    return flux;

}


StateVector computeAdvectionFluxFromConservative(StateVector conservative, Vector3D surface, const FluidBase& fluid){
    StateVector primitive = getPrimitiveVariablesFromConservative(conservative);
    StateVector flux = computeAdvectionFluxFromPrimitive(primitive, surface, fluid);
    return flux;
}


FloatType computeAngleBetweenVectors(const Vector3D& v1, const Vector3D& v2) {
    Vector3D v1norm = v1.normalized();
    Vector3D v2norm = v2.normalized();

    FloatType dotProduct = v1norm.dot(v2norm);

    // Clamp dot product to valid range for acos()
    dotProduct = std::clamp(dotProduct, static_cast<FloatType>(-1.0), static_cast<FloatType>(1.0));

    return std::acos(dotProduct);
}


Vector3D rotateVectorAlongXAxis(const Vector3D& v, FloatType theta){
    Vector3D newVector;
    newVector.x() = v.x();
    newVector.y() = v.y()*cos(theta) - v.z()*sin(theta);
    newVector.z() = v.y()*sin(theta) + v.z()*cos(theta);
    return newVector;
}

StateVector rotateStateVectorAlongXAxis(const StateVector& state, FloatType theta){
    Vector3D vector(state[1], state[2], state[3]);
    Vector3D rotatedVector = rotateVectorAlongXAxis(vector, theta);
    StateVector newVector({state[0], rotatedVector.x(), rotatedVector.y(), rotatedVector.z(), state[4]});
    return newVector;
}

Vector3D computeCylindricalComponentsFromCartesian(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y()*cos(theta) + vec.z()*sin(theta);
    newVector.z() = - vec.y()*sin(theta) + vec.z()*cos(theta);
    return newVector;
}

Vector3D computeCartesianComponentsFromCylindrical(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y() * cos(theta) - vec.z() * sin(theta);
    newVector.z() = vec.y() * sin(theta) + vec.z() * cos(theta);
    return newVector;
}


void integrateRadialEquilibrium(const std::vector<FloatType>& density,
    const std::vector<FloatType>& velocityTang,
    const std::vector<FloatType>& radius,
    std::vector<FloatType>& pressure,
    FloatType& newHubPressure)
{

    const size_t N = radius.size();
    pressure[0] = newHubPressure; 
    for (size_t i = 1; i < N; ++i) {
        FloatType r1 = radius[i - 1];
        FloatType r2 = radius[i];

        FloatType dr = r2 - r1;

        FloatType rho1 = density[i - 1];
        FloatType rho2 = density[i];

        FloatType ut1 = velocityTang[i - 1];
        FloatType ut2 = velocityTang[i];

        // Evaluate dp/dr = rho * u_theta^2 / r at both points
        FloatType dpdr1 = rho1 * ut1 * ut1 / r1;
        FloatType dpdr2 = rho2 * ut2 * ut2 / r2;

        // Trapezoidal integration step
        pressure[i] = pressure[i - 1] + 0.5 * (dpdr1 + dpdr2) * dr;
    }

}


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
    Matrix3D<Vector3D>& gradient){

    size_t ni = volumes.sizeI();
    size_t nj = volumes.sizeJ();
    size_t nk = volumes.sizeK();

    for (size_t i = 0; i < ni; i++){
        for (size_t j = 0; j < nj; j++){
            for (size_t k = 0; k < nk; k++){
                Vector3D surfaceWest = - surfacesI(i,j,k);
                Vector3D surfaceEast = surfacesI(i+1,j,k);
                Vector3D surfaceNorth = surfacesJ(i,j+1,k);
                Vector3D surfaceSouth = - surfacesJ(i,j,k);
                Vector3D surfaceBottom = - surfacesK(i,j,k);
                Vector3D surfaceTop = surfacesK(i,j,k+1);

                Vector3D centerWest = midPointsI(i,j,k);
                Vector3D centerEast = midPointsI(i+1,j,k);
                Vector3D centerNorth = midPointsJ(i,j+1,k);
                Vector3D centerSouth = midPointsJ(i,j,k);
                Vector3D centerBottom = midPointsK(i,j,k);
                Vector3D centerTop = midPointsK(i,j,k+1);

                FloatType scalarWest = interpolateScalar(field, i, j, k, nodes, centerWest, Direction3D::WEST);
                FloatType scalarEast = interpolateScalar(field, i, j, k, nodes, centerEast, Direction3D::EAST);
                FloatType scalarNorth = interpolateScalar(field, i, j, k, nodes, centerNorth, Direction3D::NORTH);
                FloatType scalarSouth = interpolateScalar(field, i, j, k, nodes, centerSouth, Direction3D::SOUTH);
                FloatType scalarBottom = interpolateScalar(field, i, j, k, nodes, centerBottom, Direction3D::BOTTOM);
                FloatType scalarTop = interpolateScalar(field, i, j, k, nodes, centerTop, Direction3D::TOP);

                std::array<Vector3D, 6> surfaces = {
                    surfaceWest, 
                    surfaceEast, 
                    surfaceNorth, 
                    surfaceSouth, 
                    surfaceBottom, 
                    surfaceTop};

                std::array<FloatType, 6> scalars = {
                    scalarWest, 
                    scalarEast, 
                    scalarNorth, 
                    scalarSouth, 
                    scalarBottom, 
                    scalarTop};

                gradient(i,j,k) = computeGreenGaussDivergence(surfaces, scalars, volumes(i,j,k));
            }
        }
    }
}


FloatType interpolateScalar(
    const Matrix3D<FloatType>& field, 
    size_t i, 
    size_t j, 
    size_t k, 
    const Matrix3D<Vector3D>& nodes, 
    const Vector3D& point, 
    Direction3D direction){

    Vector3D point0 = nodes(i,j,k);
    FloatType field0 = field(i,j,k);
    
    FloatType field1 {0.0};
    Vector3D point1(0.0,0.0,0.0);
    FloatType fieldInterpolated = 0.0;
    size_t ni = field.sizeI();
    size_t nj = field.sizeJ();
    size_t nk = field.sizeK();
    if (direction==Direction3D::WEST && i==0){
        return field0;
    }
    else if (direction==Direction3D::WEST){
        point1 = nodes(i-1,j,k);
        field1 = field(i-1,j,k);
    }
    else if (direction==Direction3D::EAST && i==ni-1){
        return field0;
    }
    else if (direction==Direction3D::EAST){
        point1 = nodes(i+1,j,k);
        field1 = field(i+1,j,k);
    }

    else if (direction==Direction3D::NORTH && j==nj-1){
        return field0;
    }
    else if (direction==Direction3D::NORTH){
        point1 = nodes(i,j+1,k);
        field1 = field(i,j+1,k);
    }
    else if (direction==Direction3D::SOUTH && j==0){
        return field0;
    }
    else if (direction==Direction3D::SOUTH){
        point1 = nodes(i,j-1,k);
        field1 = field(i,j-1,k);
    }

    else if (direction==Direction3D::BOTTOM && k==0){
        return field0;
    }
    else if (direction==Direction3D::BOTTOM){
        point1 = nodes(i,j,k-1);
        field1 = field(i,j,k-1);
    }
    else if (direction==Direction3D::TOP && k==nk-1){
        return field0;
    }
    else if (direction==Direction3D::TOP){
        point1 = nodes(i,j,k+1);
        field1 = field(i,j,k+1);
    }

    FloatType distance0 = (point - point0).magnitude();
    FloatType distance1 = (point - point1).magnitude();
    fieldInterpolated = field0 + (field1 - field0) * distance0 / (distance0 + distance1);
    return fieldInterpolated;
}


FloatType atan2FromZeroTo2pi(FloatType y, FloatType x){
    FloatType angle = std::atan2(y, x);
    if (angle < 0)
        angle += 2 * M_PI;
    return angle;
}


FloatType inverseRescalingMinMax(const FloatType& value, const FloatType& min, const FloatType& max){
    FloatType res = min + (max - min) * value;
    return res;
}


FloatType rescaleMinMax(const FloatType& value, const FloatType& min, const FloatType& max){
    FloatType res = (value - min) / (max - min);
    return res;
}




FloatType linearInterpolation(const std::vector<double>& x,
                            const std::vector<double>& y,
                            const FloatType& xp) {

    // Extrapolate with zero gradient if out of interval
    if (xp <= x.front()) {
        return y[0];
    }
    else if (xp >= x.back()) {
        size_t n = x.size();
        return y[n - 1];
    }
    else{
        // Find the interval [x[i], x[i+1]] such that x[i] <= xp < x[i+1]
        auto it = std::upper_bound(x.begin(), x.end(), xp);
        size_t i = std::distance(x.begin(), it) - 1;

        double x0 = x[i];
        double x1 = x[i + 1];
        double y0 = y[i];
        double y1 = y[i + 1];

        double t = (xp - x0) / (x1 - x0);
        return ((1.0 - t) * y0 + t * y1);
    }

    
}

void writeDataToCsv(const std::vector<FloatType>& data, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << " for writing.\n";
        return;
    }

    for (const auto& value : data)
        file << value << "\n";  // each value on its own line

    file.close();

    std::cout << "Data written to " << filename << " successfully.\n";
}


Matrix2D<FloatType> computeAdvectionJacobian(
    const StateVector& primitive, 
    const Vector3D& direction, 
    const FluidBase& fluid){

    Matrix2D<FloatType> jacobian(5, 5);

    FloatType u = primitive[1];
    FloatType v = primitive[2];
    FloatType w = primitive[3];
    FloatType et = primitive[4];
    FloatType nx = direction(0);
    FloatType ny = direction(1);
    FloatType nz = direction(2);

    // precompute common terms
    FloatType gmma = fluid.getGamma();
    FloatType umag2 = u*u + v*v + w*w; 
    FloatType V = u*nx + v*ny + w*nz; 
    FloatType phi = 0.5 * (gmma - 1) * umag2;

    // first row
    jacobian(0, 0) = 0.0;
    jacobian(0, 1) = nx;
    jacobian(0, 2) = ny;
    jacobian(0, 3) = nz;
    jacobian(0, 4) = 0.0;

    // second row
    jacobian(1, 0) = (phi - u*u)*nx - u*v*ny - u*w*nz;
    jacobian(1, 1) = (3.0 - gmma)*u*nx + v*ny + w*nz;
    jacobian(1, 2) = -(gmma - 1.0)*v*nx + u*ny + 0.0*nz;
    jacobian(1, 3) = -(gmma-1.0)*w*nx + 0.0*ny + u*nz;
    jacobian(1, 4) = (gmma-1.0)*nx + 0.0*ny + 0.0*nz;

    // third row
    jacobian(2, 0) = -u*v*nx + (phi-v*v)*ny - v*w*nz;
    jacobian(2, 1) = v*nx -(gmma-1.0)*u*ny + 0.0*nz;
    jacobian(2, 2) = u*nx + (3.0-gmma)*v*ny + w*nz;
    jacobian(2, 3) = 0.0*nx  - (gmma-1.0)*w*ny + v*nz;
    jacobian(2, 4) = 0.0*nx + (gmma-1.0)*ny + 0.0*nz;

    // fourth row
    jacobian(3, 0) = -u*w*nx - v*w*ny + (phi - w*w)*nz;
    jacobian(3, 1) = w*nx + 0.0*ny - (gmma-1.0)*u*nz;
    jacobian(3, 2) = 0.0*nx + w*ny - (gmma-1.0)*v*nz;
    jacobian(3, 3) = u*nx + v*ny + (3.0 - gmma)*w*nz;
    jacobian(3, 4) = 0.0*nx + 0.0*ny + (gmma-1.0)*nz;

    // fifth row
    jacobian(4, 0) = (-u*(gmma*et - 2.0*phi))*nx + (-v*(gmma*et - 2.0*phi))*ny + (-w*(gmma*et - 2.0*phi))*nz;
    jacobian(4, 1) = (gmma*et - phi - (gmma-1.0)*u*u)*nx - (gmma-1.0)*u*v*ny - (gmma-1.0)*u*w*nz;
    jacobian(4, 2) = -(gmma-1.0)*u*v*nx + (gmma*et - phi - (gmma-1.0)*v*v)*ny - (gmma-1.0)*v*w*nz;
    jacobian(4, 3) = -(gmma-1.0)*u*w*nx - (gmma-1.0)*v*w*ny + (gmma*et - phi - (gmma-1.0)*w*w)*nz;
    jacobian(4, 4) = gmma * V;

    return jacobian;

}