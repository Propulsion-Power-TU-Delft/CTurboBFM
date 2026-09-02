#include "advection_jst.hpp"
#include "math_utils.hpp"

StateVector AdvectionJst::computeFlux(
    const StateVector& Ull,
    const StateVector& Ul,
    const StateVector& Ur,
    const StateVector& Urr,
    const Vector3D& S) const
{
    
    auto Wll = getPrimitiveVariablesFromConservative(Ull);
    auto Wl = getPrimitiveVariablesFromConservative(Ul);
    auto Wr = getPrimitiveVariablesFromConservative(Ur);
    auto Wrr = getPrimitiveVariablesFromConservative(Urr);
    auto Wavg = (Wl + Wr) * 0.5; 

    FloatType rCoeff[2] = {
        computeRFactor(Wl),
        computeRFactor(Wr)
    };

    FloatType sCoeff[2] = {
        computeSFactor(Wll, Wl, Wr),
        computeSFactor(Wl, Wr, Wrr)
    };

    FloatType r = std::max(rCoeff[0], rCoeff[1]);
    FloatType s = std::max(sCoeff[0], sCoeff[1]);

    FloatType psi2 = _kappa2 * s * r;
    FloatType psi4 = std::max(0.0, _kappa4 * r - _c4 * psi2);

    auto flux = computeAdvectionFluxFromPrimitive(Wavg, S, _fluid);
    auto dissipation = (Ur - Ul) * psi2 - ((Urr - Ur) - (Ur - Ul) * 2 + (Ul - Ull)) * psi4;
    flux -= dissipation;

    return flux;
}

FloatType AdvectionJst::computeRFactor(const StateVector primitive) const {
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]}; 
    FloatType velMag = velocity.magnitude();
    FloatType soundSpeed = _fluid.computeSoundSpeed_rho_u_et(primitive[0], velocity, primitive[4]);
    return (velMag+soundSpeed);
}

FloatType AdvectionJst::computeSFactor(
    const StateVector prim1, 
    const StateVector prim2, 
    const StateVector prim3) const {

    FloatType p1 = _fluid.computePressure_primitive(prim1);
    FloatType p2 = _fluid.computePressure_primitive(prim2);
    FloatType p3 = _fluid.computePressure_primitive(prim3);
    FloatType sFactor = std::abs((p1 - 2*p2 + p3) / (p1 + 2*p2 + p3));
    return sFactor;
}