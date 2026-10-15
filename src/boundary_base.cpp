#include "boundary_base.hpp"
#include "math_utils.hpp"
#include "fluid_ideal.hpp"
#include <algorithm>
#include <cmath>

StateVector BoundaryBase::computeSubsonicInletFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const FloatType& totPressureBoundary, 
    const FloatType& totTemperatureBoundary, 
    const Vector3D& flowDirection){

    // properties of internal point
    StateVector primitive = getPrimitiveVariablesFromConservative(internalConservative);
    Vector3D velocityInt({primitive[1], primitive[2], primitive[3]});
    FloatType soundSpeedInt = _fluid.computeSoundSpeed_rho_u_et(primitive[0], velocityInt, primitive[4]);

    Vector3D normal = surface / surface.magnitude();
    FloatType unInt = velocityInt.dot(normal);
    FloatType Z = primitive[0] * soundSpeedInt;
    FloatType pInt = _fluid.computePressure_rho_u_et(primitive[0], velocityInt, primitive[4]);

    // Outflow normal component of inflow direction (negative for inward flow)
    FloatType Cgeom = flowDirection.dot(normal);
    if (Cgeom >= 0.0) {
        Cgeom = -1.0;
    }

    // Outgoing acoustic characteristic carried from the interior towards the boundary:
    // Along dx/dt = un + a > 0, (p + Z * un) is constant.
    // Therefore at the boundary: p_b + Z * un_b = p_int + Z * un_int
    FloatType Rout = pInt + Z * unInt;

    FloatType densityBound = 0.0;
    FloatType energyBound = 0.0;
    FloatType totEnergyBound = 0.0;
    Vector3D velocityBound;

    if (auto fluidIdeal = dynamic_cast<const FluidIdeal*>(&_fluid)) {
        FloatType gamma = fluidIdeal->getGamma();
        FloatType R = fluidIdeal->getRconstant();
        FloatType cp = gamma * R / (gamma - 1.0);
        FloatType htIn = cp * totTemperatureBoundary;

        // Subsonic pressure search range [0.05*pt, pt]
        FloatType pLow = 0.05 * totPressureBoundary;
        FloatType pHigh = 0.999999 * totPressureBoundary;
        FloatType pBound = 0.5 * (pLow + pHigh);

        for (int iter = 0; iter < 30; ++iter) {
            FloatType pMid = 0.5 * (pLow + pHigh);
            FloatType TMid = totTemperatureBoundary * std::pow(pMid / totPressureBoundary, (gamma - 1.0) / gamma);
            FloatType hMid = cp * TMid;
            FloatType VMid = std::sqrt(2.0 * std::max(static_cast<FloatType>(0.0), htIn - hMid));
            FloatType residual = pMid + Z * Cgeom * VMid - Rout;

            if (std::abs(residual) < 1e-4 * totPressureBoundary) {
                pBound = pMid;
                break;
            }
            if (residual < 0.0) {
                pLow = pMid;
            } else {
                pHigh = pMid;
            }
            pBound = pMid;
        }

        FloatType TBound = totTemperatureBoundary * std::pow(pBound / totPressureBoundary, (gamma - 1.0) / gamma);
        FloatType hBound = cp * TBound;
        FloatType velocityBoundMag = std::sqrt(2.0 * std::max(static_cast<FloatType>(0.0), htIn - hBound));
        velocityBound = flowDirection * velocityBoundMag;
        densityBound = pBound / (R * TBound);
        energyBound = fluidIdeal->computeStaticEnergy_p_rho(pBound, densityBound);
        totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);
    }
    else {
        // General acoustic characteristic boundary condition for real fluids
        FloatType sIn = _fluid.computeEntropy_p_T(totPressureBoundary, totTemperatureBoundary);
        FloatType htIn = _fluid.computeEnthalpy_p_s(totPressureBoundary, sIn);

        // Subsonic pressure search range [0.05*pt, pt]
        FloatType pLow = 0.05 * totPressureBoundary;
        FloatType pHigh = 0.999999 * totPressureBoundary;
        FloatType pBound = 0.5 * (pLow + pHigh);

        for (int iter = 0; iter < 30; ++iter) {
            FloatType pMid = 0.5 * (pLow + pHigh);
            FloatType hMid = _fluid.computeEnthalpy_p_s(pMid, sIn);
            FloatType VMid = std::sqrt(2.0 * std::max(static_cast<FloatType>(0.0), htIn - hMid));
            FloatType residual = pMid + Z * Cgeom * VMid - Rout;

            if (std::abs(residual) < 1e-4 * totPressureBoundary) {
                pBound = pMid;
                break;
            }
            if (residual < 0.0) {
                pLow = pMid;
            } else {
                pHigh = pMid;
            }
            pBound = pMid;
        }

        FloatType hBound = _fluid.computeEnthalpy_p_s(pBound, sIn);
        FloatType velocityBoundMag = std::sqrt(2.0 * std::max(static_cast<FloatType>(0.0), htIn - hBound));
        velocityBound = flowDirection * velocityBoundMag;
        densityBound = _fluid.computeDensity_p_s(pBound, sIn);
        energyBound = _fluid.computeInternalEnergy_p_s(pBound, sIn);
        totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);
    }

    // compute boundary flux
    StateVector primitiveBoundary({
        densityBound, 
        velocityBound.x(), 
        velocityBound.y(), 
        velocityBound.z(), 
        totEnergyBound});
        
    StateVector flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}

StateVector BoundaryBase::computeStagnationWallFlux(
    const Vector3D& surface, 
    FloatType totPressure, 
    FloatType totTemperature, 
    const Vector3D& wallVelocity) {

    FloatType uMag = wallVelocity.magnitude();
    FloatType pBound = totPressure;
    FloatType TBound = totTemperature;

    if (auto fluidIdeal = dynamic_cast<const FluidIdeal*>(&_fluid)) {
        FloatType gamma = fluidIdeal->getGamma();
        FloatType R = fluidIdeal->getRconstant();
        FloatType cp = gamma * R / (gamma - 1.0);
        if (uMag > 0.0) {
            TBound = totTemperature - 0.5 * uMag * uMag / cp;
            if (TBound < 0.01 * totTemperature) {
                TBound = 0.01 * totTemperature;
            }
            pBound = totPressure * std::pow(TBound / totTemperature, gamma / (gamma - 1.0));
        }
    } else {
        if (uMag > 0.0) {
            FloatType sIn = _fluid.computeEntropy_p_T(totPressure, totTemperature);
            FloatType htIn = _fluid.computeEnthalpy_p_s(totPressure, sIn);
            FloatType hBound = htIn - 0.5 * uMag * uMag;
            FloatType pLow = 0.5 * totPressure;
            FloatType pHigh = totPressure;
            for (int iter = 0; iter < 20; ++iter) {
                FloatType pMid = 0.5 * (pLow + pHigh);
                FloatType hMid = _fluid.computeEnthalpy_p_s(pMid, sIn);
                if (hMid < hBound) {
                    pLow = pMid;
                } else {
                    pHigh = pMid;
                }
            }
            pBound = 0.5 * (pLow + pHigh);
        }
    }

    Vector3D normal = surface / surface.magnitude();

    // Solid wall adheres to fluid (no penetration); convective mass and energy fluxes are zero.
    // Momentum flux is the stagnation pressure force.
    StateVector flux {};
    flux[0] = 0.0;
    flux[1] = pBound * normal.x();
    flux[2] = pBound * normal.y();
    flux[3] = pBound * normal.z();
    flux[4] = 0.0;

    return flux;
}

StateVector BoundaryBase::computeOutletFlux(
    const StateVector& internalConservative, 
    const Vector3D& surface, 
    const FloatType& boundaryPressure,
    const size_t& iterCounter){

    auto primitive = getPrimitiveVariablesFromConservative(internalConservative);
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    auto density = primitive[0];
    auto pressure = _fluid.computePressure_rho_u_et(density, velocity, primitive[4]);
    auto soundSpeed = _fluid.computeSoundSpeed_p_rho(pressure, primitive[0]);
    
    if (velocity.magnitude() >= soundSpeed) {
        auto flux = computeAdvectionFluxFromPrimitive(primitive, surface, _fluid);
        return flux;
    }
    else {
        FloatType pressureBoundary = _config.computeRampedOutletPressure(iterCounter, boundaryPressure);
        FloatType densityBoundary = 0.0;
        FloatType energyBoundary = 0.0;

        if (dynamic_cast<const FluidIdeal*>(&_fluid)) {
            densityBoundary = pressureBoundary * density / pressure;
            energyBoundary = _fluid.computeStaticEnergy_p_rho(pressureBoundary, densityBoundary);
        } else {
            // Real fluid: isentropic extrapolation along outgoing entropy s_b = s_int
            FloatType e_int = primitive[4] - 0.5 * velocity.dot(velocity);
            FloatType s_int = _fluid.computeEntropy_rho_e(density, e_int);
            densityBoundary = _fluid.computeDensity_p_s(pressureBoundary, s_int);
            energyBoundary = _fluid.computeInternalEnergy_p_s(pressureBoundary, s_int);
        }

        Vector3D velocityBoundary = velocity;
        FloatType totEnergyBoundary = energyBoundary + 0.5 * velocityBoundary.dot(velocityBoundary);
        StateVector primitiveBoundary({
            densityBoundary, 
            velocityBoundary.x(), 
            velocityBoundary.y(), 
            velocityBoundary.z(), 
            totEnergyBoundary});
        auto flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
        return flux;
    }
}


BoundaryNodesIndexRange fetchBoundaryNodesIndexRange(
    const Boundary& boundary,
    const size_t& nPointsI, 
    const size_t& nPointsJ, 
    const size_t& nPointsK)
{
    BoundaryNodesIndexRange range{
        boundary.i_min, boundary.i_max, boundary.j_min, boundary.j_max, boundary.k_min, boundary.k_max
    };


    if (range.iStart == range.iLast && range.iStart == 0) {
        range.iLast = 1;
    }
    else if (range.iStart == range.iLast && range.iStart == nPointsI) {
        range.iStart = nPointsI - 1;
        range.iLast  = range.iStart + 1;
    }
    else if (range.jStart == range.jLast && range.jStart == 0) {
        range.jLast = 1;
    }
    else if (range.jStart == range.jLast && range.jStart == nPointsJ) {
        range.jStart = nPointsJ - 1;
        range.jLast  = range.jStart + 1;
    }
    else if (range.kStart == range.kLast && range.kStart == 0) {
        range.kLast = 1;
    }
    else if (range.kStart == range.kLast && range.kStart == nPointsK) {
        range.kStart = nPointsK - 1;
        range.kLast  = range.kStart + 1;
    }
    else {
        throw std::runtime_error(
            std::string("Invalid indexing on boundary condition specs for patch ")
            + boundary.name
        );
    }

    return range;
}
