#pragma once
#include <iostream>
#include <vector>
#include <algorithm> 
#include <numeric>   
#include <cmath>
#include <map>
#include <cassert>
#include <array>
#include <cstring>

using FloatType = double;


/** @brief Class for 3D vectors of FloatType type elements. The 3 components are named as x,y,z. */
class Vector3D {

    private:
        std::array<FloatType, 3> _data; 

    public:
        
        Vector3D() : _data{0.0, 0.0, 0.0} {} 

        Vector3D(FloatType x, FloatType y, FloatType z) : _data{x, y, z} {} 

        FloatType& operator()(size_t index) {
            assert (index < 3 && "Index out of bounds");
            return _data[index];
        }

        const FloatType& operator()(size_t index) const {
            assert (index < 3 && "Index out of bounds");
            return _data[index];
        }

        FloatType& x() { return _data[0]; }
        FloatType& y() { return _data[1]; }
        FloatType& z() { return _data[2]; }

        const FloatType& x() const { return _data[0]; }
        const FloatType& y() const { return _data[1]; }
        const FloatType& z() const { return _data[2]; }

        Vector3D operator+(const Vector3D& other) const {
            return Vector3D(x() + other.x(), y() + other.y(), z() + other.z());
        }

        Vector3D operator-(const Vector3D& other) const {
            return Vector3D(x() - other.x(), y() - other.y(), z() - other.z());
        }

        Vector3D operator*(FloatType scalar) const {
            return Vector3D(x() * scalar, y() * scalar, z() * scalar);
        }

        Vector3D operator/(FloatType scalar) const {
            return Vector3D(x() / scalar, y() / scalar, z() / scalar);
        }

        Vector3D& operator+=(const Vector3D& other) {
            x() += other.x(); 
            y() += other.y(); 
            z() += other.z();
            return *this;
        }

        Vector3D& operator-=(const Vector3D& other) {
            x() -= other.x(); 
            y() -= other.y(); 
            z() -= other.z();
            return *this;
        }

        Vector3D operator-() const {
            return Vector3D(-x(), -y(), -z());
        }

        Vector3D& operator*=(FloatType scalar) {
            x() *= scalar; 
            y() *= scalar; 
            z() *= scalar;
            return *this;
        }

        Vector3D& operator/=(FloatType scalar) {
            x() /= scalar; 
            y() /= scalar; 
            z() /= scalar;
            return *this;
        }

        bool operator==(const Vector3D& other) const {
            return (x() == other.x() && y() == other.y() && z() == other.z());
        }

        FloatType magnitude() const {
            return std::sqrt(x()*x() + y()*y() + z()*z());
        }

        Vector3D normalized() const {
            FloatType mag = magnitude();
            return (mag == 0) ? Vector3D(0, 0, 0) : *this / mag;
        }

        Vector3D cross(const Vector3D& other) const {
            return Vector3D(
                y()*other.z() - z()*other.y(),
                z()*other.x() - x()*other.z(),
                x()*other.y() - y()*other.x()
            );
        }

        FloatType dot(const Vector3D& other) const {
            FloatType result = x()*other.x() + y()*other.y() + z()*other.z();
            return result;
        }

        friend std::ostream& operator<<(std::ostream& os, const Vector3D& v) {
        return os << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
        }
};



enum class BodyForceModel {
    NONE,
    HALL,
    HALL_THOLLET,
    CHIMA,
    ONLY_BLOCKAGE,
    CORRELATIONS,
    LIFT_DRAG,
};

enum class TurbulenceModel {
    NONE,
    SPALART_ALLMARAS,
};

enum class BoundaryIndex {
    I_START,
    I_END,
    J_START,
    J_END,
    K_START,
    K_END
};


enum class OutputFormat {
    CSV = 0,
    VTK = 1
};


enum class KindSolver {
    EULER = 0,
    RANS = 1
};


enum class Topology{
    TWO_DIMENSIONAL,
    THREE_DIMENSIONAL,
    AXISYMMETRIC,
    ONE_DIMENSIONAL
};


enum class BoundaryType {
    INLET,
    INLET_SUPERSONIC,
    OUTLET_SUPERSONIC,
    OUTLET,
    RADIAL_EQUILIBRIUM,
    INVISCID_WALL,
    EMPTY,
    WEDGE,
    PERIODIC,
    THROTTLE,
    INLET_2D,
    NO_SLIP_WALL,
    TRANSPARENT
};


enum class TimeIntegration {
    RUNGE_KUTTA_4,
    RUNGE_KUTTA_3
};


enum class TimeStepMethod {
    LOCAL,
    GLOBAL,
    FIXED
};

enum class TurboPerformance {
    MASS_FLOW,
    TOTAL_PRESSURE_RATIO,
    TOTAL_EFFICIENCY,
    TOTAL_TEMPERATURE_RATIO,
};

enum class MonitorOutputField {
    PRESSURE,
    VELOCITY_X,
    VELOCITY_Y,
    VELOCITY_Z,
    TIME,
};

enum class FluxLimiter {
    NONE,
    VAN_ALBADA,
    VAN_LEER,
    MIN_MOD,
};

enum class AdvectionScheme {
    JST,
    ROE,
};


/** @brief Class for 2D data container generic elements*/
template<typename T>
class Matrix2D {
    public:
        Matrix2D() : _ni(0), _nj(0) {}

        Matrix2D(size_t ni, size_t nj) : _ni(ni), _nj(nj), _data(ni * nj, T{}) {}

        void resize(size_t ni, size_t nj) {
            _ni = ni;
            _nj = nj;
            _data.assign(ni * nj, T{});
        }

        T& operator()(size_t i, size_t j) {
            check_bounds(i, j);
            return _data[i * _nj + j];
        }

        const T& operator()(size_t i, size_t j) const {
            check_bounds(i, j);
            return _data[i * _nj + j];
        }

        size_t sizeI() const { return _ni; }
        size_t sizeJ() const { return _nj; }

        Matrix2D<T> operator+(const Matrix2D<T>& other) const {
            assert(_ni == other._ni && _nj == other._nj && "Matrix dimensions must match for addition");
            Matrix2D<T> result(_ni, _nj);
            for (size_t n = 0; n < _data.size(); ++n)
                result._data[n] = _data[n] + other._data[n];
            return result;
        }

        // Regular subtraction (A - B)
        Matrix2D<T> operator-(const Matrix2D<T>& other) const {
            assert(_ni == other._ni && _nj == other._nj && "Matrix dimensions must match for addition");
            Matrix2D<T> result(_ni, _nj);
            for (size_t n = 0; n < _data.size(); ++n)
                result._data[n] = _data[n] - other._data[n];
            return result;
        }

        // Matrix * scalar
        Matrix2D<T> operator*(const T& scalar) const {
            Matrix2D<T> result(_ni, _nj);
            for (size_t n = 0; n < _data.size(); ++n)
                result._data[n] = _data[n] * scalar;
            return result;
        }

        // Matrix * Matrix multiplication
        Matrix2D<T> operator*(const Matrix2D<T>& other) const {
            assert(_nj == other._ni && "Matrix dimensions must be compatible for multiplication");

            Matrix2D<T> result(_ni, other._nj);

            for (size_t i = 0; i < _ni; ++i) {
                for (size_t j = 0; j < other._nj; ++j) {
                    T sum = T{};  // initialize to zero
                    for (size_t k = 0; k < _nj; ++k) {
                        sum += (*this)(i, k) * other(k, j);
                    }
                    result(i, j) = sum;
                }
            }

            return result;
        }

    private:
        size_t _ni, _nj;
        std::vector<T> _data;

        void check_bounds(size_t i, size_t j) const {
            assert(i < _ni && "i-index out of range");
            assert(j < _nj && "j-index out of range");
        }
};



/** @brief Class for 3D data container */
template <typename T>
class Matrix3D {
    public:
        Matrix3D() : _ni(0), _nj(0), _nk(0) {}

        Matrix3D(size_t ni, size_t nj, size_t nk, const std::vector<T>& initial_data = {})
            : _ni(ni), _nj(nj), _nk(nk) {
            if (!initial_data.empty()) {
                if (initial_data.size() != ni * nj * nk) {
                    throw std::invalid_argument("Initial data size does not match matrix dimensions");
                }
                _data = initial_data;
            } else {
                _data.resize(ni * nj * nk);
            }
        }

        void resize(size_t ni, size_t nj, size_t nk) {
            _ni = ni;
            _nj = nj;
            _nk = nk;
            _data.resize(ni * nj * nk);
        }

        void copyFrom(const Matrix3D<T>& other) {
            assert(_data.size() == other._data.size() && "Matrix3D::copyFrom data size mismatch");
            std::memcpy(_data.data(), other._data.data(), _data.size() * sizeof(T));
        }

        T& operator()(size_t i, size_t j, size_t k) {
            check_bounds(i, j, k);
            return _data[i * _nj * _nk + j * _nk + k];
        }

        const T& operator()(size_t i, size_t j, size_t k) const {
            check_bounds(i, j, k);
            return _data[i * _nj * _nk + j * _nk + k];
        }

        size_t sizeI() const { return _ni; }
        size_t sizeJ() const { return _nj; }
        size_t sizeK() const { return _nk; }

        // Compute L2 norm of all the matrix elements
        T norm() const {
            T sum = 0;
            for (size_t i = 0; i < _ni; ++i) {
                for (size_t j = 0; j < _nj; ++j) {
                    for (size_t k = 0; k < _nk; ++k) {
                        sum += std::pow(operator()(i, j, k), 2); 
                    }
                }
            }
            return std::sqrt(sum); 
        }

        Matrix3D& operator+=(const Matrix3D& other) {
            assert(_ni == other._ni && "Matrix3D::operator+= dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+= dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+= dimension mismatch on nk");
            for (size_t idx = 0; idx < _data.size(); ++idx) {
                _data[idx] += other._data[idx];
            }
            return *this;
        }

        void setToZero() {
            size_t n = _data.size();
            for (size_t idx = 0; idx < n; ++idx) {
                _data[idx] = T(0);
            }
        }

        Matrix3D& operator*=(const Matrix3D& other) {
            assert(_ni == other._ni && "Matrix3D::operator+= dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+= dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+= dimension mismatch on nk");
            for (size_t idx = 0; idx < _data.size(); ++idx) {
                _data[idx] *= other._data[idx];
            }
            return *this;
        }

        Matrix3D& operator/=(const Matrix3D& other) {
            assert(_ni == other._ni && "Matrix3D::operator+= dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+= dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+= dimension mismatch on nk");
            for (size_t idx = 0; idx < _data.size(); ++idx) {
                _data[idx] /= other._data[idx];
            }
            return *this;
        }

        Matrix3D operator*(const Matrix3D& other) const {
            assert(_ni == other._ni && "Matrix3D::operator* dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator* dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator* dimension mismatch on nk");

            Matrix3D result(*this); // copy current matrix
            result *= other;        // use your operator*=
            return result;
        }

        Matrix3D operator+(const Matrix3D& other) const {
            assert(_ni == other._ni && "Matrix3D::operator+ dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+ dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+ dimension mismatch on nk");

            Matrix3D result(*this); // copy current matrix
            result += other;        // use your operator*=
            return result;
        }

        // In-place division
        Matrix3D<T> operator/(const Matrix3D<FloatType>& other) const {
            assert(_ni == other._ni && "Matrix3D::operator+= dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+= dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+= dimension mismatch on nk");
            Matrix3D<FloatType> result(_ni, _nj, _nk);
            for (size_t i = 0; i < _ni; ++i) {
                for (size_t j = 0; j < _nj; ++j) {
                    for (size_t k = 0; k < _nk; ++k) {
                        if (other(i, j, k) == 0) throw std::runtime_error("Division by zero");
                        result(i, j, k) = (*this)(i, j, k) / other(i, j, k);
                    }
                }
            }
            return result;
        }

        Matrix3D& operator*=(T scalar) {
            for (size_t idx = 0; idx < _data.size(); ++idx) {
                _data[idx] *= scalar;
            }
            return *this;
        }

        Matrix3D operator*(T scalar) const {
            Matrix3D result(*this);  
            result *= scalar;        
            return result;
        }

        Matrix3D& operator-=(const Matrix3D& other) {
            assert(_ni == other._ni && "Matrix3D::operator+= dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator+= dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator+= dimension mismatch on nk");
            for (size_t idx = 0; idx < _data.size(); ++idx) {
                _data[idx] -= other._data[idx];
            }
            return *this;
        }

        Matrix3D operator-(const Matrix3D& other) const {
            assert(_ni == other._ni && "Matrix3D::operator- dimension mismatch on ni");
            assert(_nj == other._nj && "Matrix3D::operator- dimension mismatch on nj");
            assert(_nk == other._nk && "Matrix3D::operator- dimension mismatch on nk");

            Matrix3D result(*this);  
            result -= other;         
            return result;
        }

        std::vector<T> getData() const {return _data;}

        // get a slice of one of the boundary
        Matrix2D<T> getBoundarySlice(BoundaryIndex index) const {
            Matrix2D<T> slice;
            size_t ni = sizeI();
            size_t nj = sizeJ();
            size_t nk = sizeK();
            
            // i slices
            if (index==BoundaryIndex::I_START || index==BoundaryIndex::I_END){
                
                slice.resize(nj, nk);

                if (index==BoundaryIndex::I_START){
                    for (size_t j=0; j<nj; j++){
                        for (size_t k=0; k<nk; k++){
                            slice(j,k) = (*this)(0, j, k);
                        }
                    }
                } else {
                    for (size_t j=0; j<nj; j++){
                        for (size_t k=0; k<nk; k++){
                            slice(j,k) = (*this)(ni-1, j, k);
                        }
                    }
                }
            }

            // j slices
            if (index==BoundaryIndex::J_START || index==BoundaryIndex::J_END){
                slice.resize(ni, nk);

                if (index==BoundaryIndex::J_START){
                    for (size_t i=0; i<ni; i++){
                        for (size_t k=0; k<nk; k++){
                            slice(i,k) = (*this)(i, 0, k);
                        }
                    }
                } else {
                    for (size_t i=0; i<ni; i++){
                        for (size_t k=0; k<nk; k++){
                            slice(i,k) = (*this)(i, nj-1, k);
                        }
                    }
                }
            }

            // k slices
            if (index==BoundaryIndex::K_START || index==BoundaryIndex::K_END){
                slice.resize(ni, nj);

                if (index==BoundaryIndex::K_START){
                    for (size_t i=0; i<ni; i++){
                        for (size_t j=0; j<nj; j++){
                            slice(i,j) = (*this)(i, j, 0);
                        }
                    }
                } else {
                    for (size_t i=0; i<ni; i++){
                        for (size_t j=0; j<nj; j++){
                            slice(i,j) = (*this)(i, j, nk-1);
                        }
                    }
                }
            }

            return slice;
        }

        T min() const {
            if (_data.empty()) {
                throw std::runtime_error("Matrix3D::min called on empty matrix");
            }
            return *std::min_element(_data.begin(), _data.end());
        }
        
        T max() const {
            if (_data.empty()) {
                throw std::runtime_error("Matrix3D::max called on empty matrix");
            }
            return *std::max_element(_data.begin(), _data.end());
        }

    private:
        size_t _ni, _nj, _nk;
        std::vector<T> _data;

        void check_bounds(size_t i, size_t j, size_t k) const {
            assert(i < _ni && "i index out of range");
            assert(j < _nj && "j index out of range");
            assert(k < _nk && "k index out of range");
        }
};


/** \brief State vector class (conservative or primitive), holding 5 float values, with overloaded operators */
class StateVector {
    private:
        static constexpr std::size_t Size = 5; 
        std::array<FloatType, Size> _data;
    
    public:
        StateVector() = default;
    
        explicit StateVector(const std::array<FloatType, Size>& arr) : _data(arr) {}

        size_t size() const { return Size; }
    
        FloatType& operator[](std::size_t idx) { 
            return _data[idx]; }

        const FloatType& operator[](std::size_t idx) const { 
            return _data[idx]; }
    
        StateVector operator+(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] + other[i];
            return result;
        }
    
        StateVector operator-(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] - other[i];
            return result;
        }
    
        StateVector operator*(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] * other[i];
            return result;
        }
    
        StateVector operator/(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i) {
                if (other[i] == 0.0)
                    throw std::runtime_error("Division by zero in StateVector");
                result[i] = _data[i] / other[i];
            }
            return result;
        }
    
        StateVector operator*(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] * scalar;
            return result;
        }
    
        StateVector operator/(FloatType scalar) const {
            if (scalar == 0.0)
                throw std::runtime_error("Division by zero in scalar operation");
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] / scalar;
            return result;
        }
    
        StateVector operator+(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] + scalar;
            return result;
        }
    
        StateVector operator-(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] - scalar;
            return result;
        }
    
        StateVector& operator+=(const StateVector& other) {
            for (std::size_t i = 0; i < Size; ++i)
                _data[i] += other[i];
            return *this;
        }
    
        StateVector& operator-=(const StateVector& other) {
            for (std::size_t i = 0; i < Size; ++i)
                _data[i] -= other[i];
            return *this;
        }
    
        operator std::array<FloatType, Size>() const {
            return _data;
        }
    
        FloatType dot(const StateVector& other) const {
            FloatType result = 0.0;
            for (std::size_t i = 0; i < Size; ++i)
                result += _data[i] * other[i];
            return result;
        }

    };


/** \brief FlowSolution class, holding 3D arrays of conserviative variables */
struct FlowSolution {
    Matrix3D<FloatType> _rho;
    Matrix3D<FloatType> _rhoU;
    Matrix3D<FloatType> _rhoV;
    Matrix3D<FloatType> _rhoW;
    Matrix3D<FloatType> _rhoE;

    FlowSolution(size_t ni, size_t nj, size_t nk) {
        resize(ni, nj, nk);
    }

    FlowSolution() = default;

    FlowSolution(const FlowSolution& other)
        : _rho(other._rho),
        _rhoU(other._rhoU),
        _rhoV(other._rhoV),
        _rhoW(other._rhoW),
        _rhoE(other._rhoE)
    {}

    FloatType norm(size_t iVar) const {
        switch (iVar)
        {
        case 0:
            return _rho.norm();
        case 1:
            return _rhoU.norm();
        case 2:
            return _rhoV.norm();
        case 3:
            return _rhoW.norm();
        case 4:
            return _rhoE.norm();        
        default:
            throw std::out_of_range("Index out of range for norm calculation");
        }
    }

    void setToZero() {
        _rho.setToZero();
        _rhoU.setToZero();
        _rhoV.setToZero();
        _rhoW.setToZero();
        _rhoE.setToZero();
    }

    void copyFrom(const FlowSolution& other) {
        _rho = other._rho;
        _rhoU = other._rhoU;
        _rhoV = other._rhoV;
        _rhoW = other._rhoW;
        _rhoE = other._rhoE;
    }

    Matrix3D<FloatType> getSolution(size_t iVar) const {
        switch (iVar)
        {
        case 0:
            return _rho;
        case 1:
            return _rhoU;
        case 2:
            return _rhoV;
        case 3:
            return _rhoW;
        case 4:
            return _rhoE;        
        default:
            throw std::out_of_range("Index out of range");
        }
    }

    Matrix3D<FloatType> getDensity() const {
        return _rho;
    }

    Matrix3D<FloatType> getVelocityX() const {
        auto results = _rhoU / _rho;
        return results;
    }

    Matrix3D<FloatType> getVelocityY() const {
        auto results = _rhoV / _rho;
        return results;
    }

    Matrix3D<FloatType> getVelocityZ() const {
        auto results = _rhoW / _rho;
        return results;
    }

    Matrix3D<FloatType> getTotalEnergy() const {
        auto results = _rhoE / _rho;
        return results;
    }

    void add(size_t i, size_t j, size_t k, const StateVector& delta) {
        _rho(i,j,k)   += delta[0];
        _rhoU(i,j,k)  += delta[1];
        _rhoV(i,j,k)  += delta[2];
        _rhoW(i,j,k)  += delta[3];
        _rhoE(i,j,k)  += delta[4];
    }

    void subtract(size_t i, size_t j, size_t k, const StateVector& delta) {
        _rho(i,j,k)   -= delta[0];
        _rhoU(i,j,k)  -= delta[1];
        _rhoV(i,j,k)  -= delta[2];
        _rhoW(i,j,k)  -= delta[3];
        _rhoE(i,j,k)  -= delta[4];
    }

    StateVector at(size_t i, size_t j, size_t k) const {
        return StateVector({_rho(i,j,k), _rhoU(i,j,k), _rhoV(i,j,k), _rhoW(i,j,k), _rhoE(i,j,k)});
    }

    void set(size_t i, size_t j, size_t k, const StateVector& vals) {
        _rho(i,j,k) = vals[0];
        _rhoU(i,j,k) = vals[1];
        _rhoV(i,j,k) = vals[2];
        _rhoW(i,j,k) = vals[3];
        _rhoE(i,j,k) = vals[4];
    }

    void set(size_t i, size_t j, size_t k, const size_t iVar, FloatType val) {
        switch (iVar)
        {
        case 0:
            _rho(i,j,k) = val;
            break;
        case 1:
            _rhoU(i,j,k) = val;
            break;
        case 2:
            _rhoV(i,j,k) = val;
            break;
        case 3:
            _rhoW(i,j,k) = val;
            break;
        case 4:
            _rhoE(i,j,k) = val;
            break;
        default:
            break;
        }
    }

    void resize(size_t i, size_t j, size_t k){
        _rho.resize(i,j,k);
        _rhoU.resize(i,j,k);
        _rhoV.resize(i,j,k);
        _rhoW.resize(i,j,k);
        _rhoE.resize(i,j,k);
    }

    FlowSolution operator+(const FlowSolution& other) const {
        FlowSolution result = *this;
        result += other;
        return result;
    }

    FlowSolution& operator+=(const FlowSolution& other) {
        _rho  += other._rho;
        _rhoU += other._rhoU;
        _rhoV += other._rhoV;
        _rhoW += other._rhoW;
        _rhoE += other._rhoE;
        return *this;
    }

    FlowSolution& operator*(const Matrix3D<FloatType>& other) {
        _rho  *= other;
        _rhoU += other;
        _rhoV += other;
        _rhoW += other;
        _rhoE += other;
        return *this;
    }

    FlowSolution& operator/(const Matrix3D<FloatType>& other) {
        _rho  /= other;
        _rhoU /= other;
        _rhoV /= other;
        _rhoW /= other;
        _rhoE /= other;
        return *this;
    }

    FlowSolution operator-(const FlowSolution& other) const {
        FlowSolution result = *this;
        result -= other;
        return result;
    }

    FlowSolution& operator-=(const FlowSolution& other) {
        _rho  -= other._rho;
        _rhoU -= other._rhoU;
        _rhoV -= other._rhoV;
        _rhoW -= other._rhoW;
        _rhoE -= other._rhoE;
        return *this;
    }

    FlowSolution& operator*=(FloatType scalar) {
        _rho  *= scalar;
        _rhoU *= scalar;
        _rhoV *= scalar;
        _rhoW *= scalar;
        _rhoE *= scalar;
        return *this;
    }

    FlowSolution operator*(FloatType scalar) const {
        FlowSolution result = *this;
        result *= scalar;  
        return result;
    }
};


enum class FluxDirection {
    I,
    J,
    K,
};

enum class SolutionName {
    DENSITY=0,
    VELOCITY_X=1,
    VELOCITY_Y=2,
    VELOCITY_Z=3,
    TOTAL_ENERGY=4,
    TEMPERATURE=5
};

enum class InputField {
    X_COORDS,
    Y_COORDS,
    Z_COORDS,
    RHO,
    RHO_U,
    RHO_V,
    RHO_W,
    RHO_E,
    BLOCKAGE,
    NORMAL_AXIAL,
    NORMAL_RADIAL,
    NORMAL_TANGENTIAL,
    RPM,
    STREAMWISE_LENGTH,
    BLADE_PRESENT,
    NUMBER_BLADES,
    SPANWISE_LENGTH,
    DELTA_TOT_ENTHALPY_DM,
    DELTA_ENTROPY_DM,
    CHIMA_MASS_FLOW,
    CHIMA_SCALING_TURNING,
    CHIMA_SCALING_LOSS,
    CHIMA_SCALING_DEVIATION,
    BLADE_METAL_ANGLE,
    BLADE_GAS_PATH_ANGLE,
    BLADE_LEAN_ANGLE,
    BLADE_CAMBER_CURVATURE,
    TOTAL_PRESSURE,
    TOTAL_TEMPERATURE,
    INLET_NX,
    INLET_NY,
    INLET_NZ,
    BLADE_DANGLE_DMERIDIONAL,
    DEVIATION_ANGLE_PIVOT,
    DELTA_ANGULAR_MOMENTUM,
    DELTA_ENTROPY,
    DELTA_FLOW_ANGLE,
    DEVIATION_ANGLE_STALL,
    LIFT_DRAG_BFM_BETA0,
    LIFT_DRAG_BFM_BETA_ETAMAX,
    LIFT_DRAG_BFM_KP_ETAMAX,
};


class FieldNameMapper {
public:
    static const std::map<std::string, InputField>& stringToEnum() {
        static const std::map<std::string, InputField> map = {
            {"x", InputField::X_COORDS},
            {"y", InputField::Y_COORDS},
            {"z", InputField::Z_COORDS},
            {"rho", InputField::RHO},
            {"rho_u", InputField::RHO_U},
            {"rho_v", InputField::RHO_V},
            {"rho_w", InputField::RHO_W},
            {"rho_e", InputField::RHO_E},
            {"blockage", InputField::BLOCKAGE},
            {"normalAxial", InputField::NORMAL_AXIAL},
            {"normalRadial", InputField::NORMAL_RADIAL},
            {"normalTangential", InputField::NORMAL_TANGENTIAL},
            {"rpm", InputField::RPM},
            {"streamwiseLength", InputField::STREAMWISE_LENGTH},
            {"bladePresent", InputField::BLADE_PRESENT},
            {"numberBlades", InputField::NUMBER_BLADES},
            {"spanwiseLength", InputField::SPANWISE_LENGTH},
            {"Dht_Dm", InputField::DELTA_TOT_ENTHALPY_DM},
            {"Ds_Dm", InputField::DELTA_ENTROPY_DM},
            {"MassFlow", InputField::CHIMA_MASS_FLOW},
            {"PhiTurn", InputField::CHIMA_SCALING_TURNING},
            {"PhiLoss", InputField::CHIMA_SCALING_LOSS},
            {"PhiDev", InputField::CHIMA_SCALING_DEVIATION},
            {"bladeMetalAngle", InputField::BLADE_METAL_ANGLE},
            {"bladeGasPathAngle", InputField::BLADE_GAS_PATH_ANGLE},
            {"bladeLeanAngle", InputField::BLADE_LEAN_ANGLE},
            {"bladeCamberCurvature", InputField::BLADE_CAMBER_CURVATURE},
            {"Pt", InputField::TOTAL_PRESSURE},
            {"Tt", InputField::TOTAL_TEMPERATURE},
            {"nx", InputField::INLET_NX},
            {"ny", InputField::INLET_NY},
            {"nz", InputField::INLET_NZ},
            {"dbladeMetalAngle_dm", InputField::BLADE_DANGLE_DMERIDIONAL},
            {"deviationAnglePivot", InputField::DEVIATION_ANGLE_PIVOT},
            {"deltaAngularMomentum", InputField::DELTA_ANGULAR_MOMENTUM},
            {"deltaEntropy", InputField::DELTA_ENTROPY},
            {"deltaFlowAngle", InputField::DELTA_FLOW_ANGLE},
            {"deviationAngleStall", InputField::DEVIATION_ANGLE_STALL},
            {"liftDragBeta0", InputField::LIFT_DRAG_BFM_BETA0},
            {"liftDragBetaEtamax", InputField::LIFT_DRAG_BFM_BETA_ETAMAX},
            {"liftDragKpEtamax", InputField::LIFT_DRAG_BFM_KP_ETAMAX}
        };
        return map;
    }

    static const std::map<InputField, std::string>& enumToString() {
        static const std::map<InputField, std::string> map = {
            {InputField::X_COORDS, "x"},
            {InputField::Y_COORDS, "y"},
            {InputField::Z_COORDS, "z"},
            {InputField::RHO, "rho"},
            {InputField::RHO_U, "rho_u"},
            {InputField::RHO_V, "rho_v"},
            {InputField::RHO_W, "rho_w"},
            {InputField::RHO_E, "rho_e"},
            {InputField::BLOCKAGE, "blockage"},
            {InputField::NORMAL_AXIAL, "normalAxial"},
            {InputField::NORMAL_RADIAL, "normalRadial"},
            {InputField::NORMAL_TANGENTIAL, "normalTangential"},
            {InputField::RPM, "rpm"},
            {InputField::STREAMWISE_LENGTH, "streamwiseLength"},
            {InputField::BLADE_PRESENT, "bladePresent"},
            {InputField::NUMBER_BLADES, "numberBlades"},
            {InputField::SPANWISE_LENGTH, "spanwiseLength"},
            {InputField::DELTA_TOT_ENTHALPY_DM, "Dht_Dm"},
            {InputField::DELTA_ENTROPY_DM, "Ds_Dm"},
            {InputField::CHIMA_MASS_FLOW, "MassFlow"},
            {InputField::CHIMA_SCALING_TURNING, "PhiTurn"},
            {InputField::CHIMA_SCALING_LOSS, "PhiLoss"},
            {InputField::CHIMA_SCALING_DEVIATION, "PhiDev"},
            {InputField::BLADE_METAL_ANGLE, "bladeMetalAngle"},
            {InputField::BLADE_GAS_PATH_ANGLE, "bladeGasPathAngle"},
            {InputField::BLADE_LEAN_ANGLE, "bladeLeanAngle"},
            {InputField::BLADE_CAMBER_CURVATURE, "bladeCamberCurvature"},
            {InputField::TOTAL_PRESSURE, "Pt"},
            {InputField::TOTAL_TEMPERATURE, "Tt"},
            {InputField::INLET_NX, "nx"},
            {InputField::INLET_NY, "ny"},
            {InputField::INLET_NZ, "nz"},
            {InputField::BLADE_DANGLE_DMERIDIONAL, "dbladeMetalAngle_dm"},
            {InputField::DEVIATION_ANGLE_PIVOT, "deviationAnglePivot"},
            {InputField::DELTA_ANGULAR_MOMENTUM, "deltaAngularMomentum"},
            {InputField::DELTA_ENTROPY, "deltaEntropy"},
            {InputField::DELTA_FLOW_ANGLE, "deltaFlowAngle"},
            {InputField::DEVIATION_ANGLE_STALL, "deviationAngleStall"},
            {InputField::LIFT_DRAG_BFM_BETA0, "liftDragBeta0"},
            {InputField::LIFT_DRAG_BFM_BETA_ETAMAX, "liftDragBetaEtamax"},
            {InputField::LIFT_DRAG_BFM_KP_ETAMAX, "liftDragKpEtamax"}
        };
        return map;
    }

    static InputField fromString(const std::string& str) {
        const auto& map = stringToEnum();
        auto it = map.find(str);
        if (it == map.end()) {
            throw std::runtime_error("Unknown field name: " + str);
        }
        return it->second;
    }

    static std::string toString(InputField field) {
        const auto& map = enumToString();
        auto it = map.find(field);
        if (it == map.end()) {
            throw std::runtime_error("Unknown FieldNames enum value");
        }
        return it->second;
    }
};


struct Statistics {
    FloatType mean = 0;
    FloatType min = 0;
    FloatType max = 0;
    FloatType stdDev = 0;
    size_t size = 0;

    Statistics() = default;

    Statistics(const std::vector<FloatType>& data) {
        computeFromData(data);
    }

    Statistics(const Matrix3D<FloatType>& matrix) {
        if (matrix.sizeI() * matrix.sizeJ() * matrix.sizeK() == 0) {
            mean = min = max = stdDev = 0;
            size = 0;
            return;
        }
        computeFromData(matrix.getData());
    }

    void printInfo(){
        std::cout << "Stats: {\n";
        std::cout << "min: " << min << std::endl;
        std::cout << "max: " << max << std::endl;
        std::cout << "mean: " << mean << std::endl;
        std::cout << "std: " << stdDev << std::endl;
        std::cout << "size: " << size << std::endl;
        std::cout << "}\n\n";
    }

private:
    void computeFromData(const std::vector<FloatType>& data) {
        auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
        min = *minIt;
        max = *maxIt;

        FloatType sum = std::accumulate(data.begin(), data.end(), static_cast<FloatType>(0));
        mean = sum / static_cast<FloatType>(data.size());

        FloatType sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), static_cast<FloatType>(0));
        stdDev = std::sqrt(sq_sum / static_cast<FloatType>(data.size()) - mean * mean);

        size = data.size();
    }
};

enum class Direction3D {
    WEST,
    EAST,
    NORTH,
    SOUTH,
    BOTTOM,
    TOP
};

enum class FluidModel {
    IDEAL,
    REAL
};

enum class ReferenceFrame {
    CARTESIAN,
    CYLINDRICAL
};

enum class OutputFields {
    PRIMARY,
    SECONDARY,
    TURBO_BFM
};