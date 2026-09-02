#pragma once

#include "advection_base.hpp"


class AdvectionRoe : public AdvectionBase {
public:

    AdvectionRoe(const Config &config, FluidBase& fluid)  : AdvectionBase(config, fluid) {}

    StateVector computeFlux(
        const StateVector& Ull,
        const StateVector& Ul,
        const StateVector& Ur,
        const StateVector& Urr,
        const Vector3D& S) const override;

private:
    struct RoeState {
        Vector3D n1, n2, n3;
        FloatType rhoAVG{0}, u1AVG{0}, u2AVG{0}, u3AVG{0}, htAVG{0}, aAVG{0};
        StateVector eigenvalues{};
        std::array<StateVector, 5> eigenvectors{};
        StateVector waveStrengths{};
    };

    /** 
     * @brief Computes the orthonormal triad (normal, tangential, binormal) associated with the surface vector S.
     * The normal is aligned with S (left-to-right orientation).
     * */ 
    void computeNormalTriad(const Vector3D& S, RoeState& state) const;

    /**
     * @brief Rotate the primitive variables into the local coordinate system defined by the normal triad.
     */
    StateVector computeRotatedPrimitive(const StateVector& W, const RoeState& state) const;

    void computeRoeAvgVariables(const StateVector& WnormL, const StateVector& WnormR, RoeState& state) const;

    /**
     * @brief Compute Roe avg values of variable phi.
     */
    FloatType roeAverage(FloatType rhoL, FloatType rhoR, FloatType phiL, FloatType phiR) const;

    void computeEigenvalues(RoeState& state) const;

    void computeEigenvectors(RoeState& state) const;

    void computeWaveStrengths(const StateVector& WnormL, const StateVector& WnormR, RoeState& state) const;

    void assembleTotalFlux(
        const Vector3D& S, 
        const StateVector& WnormL, 
        const StateVector& WnormR, 
        StateVector& flux,
        const RoeState& state) const;

private:
    Vector3D _xVersor {1.0, 0.0, 0.0}; 
    Vector3D _yVersor {0.0, 1.0, 0.0}; 
    Vector3D _zVersor {0.0, 0.0, 1.0};
    FloatType _entropyFixCoefficient {0.01}; // default value
};
