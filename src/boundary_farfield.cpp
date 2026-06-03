#include "boundary_farfield.hpp"
#include "math_utils.hpp"

StateVector BoundaryFarfield::computeBoundaryFlux(
    const StateVector& internalConservative,
    const Vector3D& surface,
    const Vector3D& midPoint,
    const std::array<size_t, 3>& indices,
    const FlowSolution& flowSolution,
    const size_t& iterCounter) {

    const FloatType gamma = _fluid.getGamma();

    // freestream (reference) state from config: [Mach, p, T, dx, dy, dz]
    FloatType machInf = _boundaryValues[0];
    FloatType pressureInf = _boundaryValues[1];
    FloatType temperatureInf = _boundaryValues[2];
    Vector3D flowDirection({_boundaryValues[3], _boundaryValues[4], _boundaryValues[5]});
    flowDirection /= flowDirection.magnitude();

    FloatType densityInf = _fluid.computeDensity_p_T(pressureInf, temperatureInf);
    FloatType soundSpeedInf = _fluid.computeSoundSpeed_p_rho(pressureInf, densityInf);
    Vector3D velocityInf = flowDirection * (machInf * soundSpeedInf);

    // internal point state
    StateVector primitiveInt = getPrimitiveVariablesFromConservative(internalConservative);
    FloatType densityInt = primitiveInt[0];
    Vector3D velocityInt({primitiveInt[1], primitiveInt[2], primitiveInt[3]});
    FloatType pressureInt = _fluid.computePressure_rho_u_et(densityInt, velocityInt, primitiveInt[4]);
    FloatType soundSpeedInt = _fluid.computeSoundSpeed_p_rho(pressureInt, densityInt);

    // outward-pointing unit normal of the boundary face
    Vector3D normal = surface / surface.magnitude();

    // normal velocity components (positive = outflow, leaving the domain)
    FloatType normalVelInt = velocityInt.dot(normal);
    FloatType normalVelInf = velocityInf.dot(normal);

    StateVector primitiveBoundary;

    if (normalVelInt >= soundSpeedInt) {
        // supersonic outflow: boundary state fully determined by the interior
        primitiveBoundary = primitiveInt;
    }
    else if (normalVelInt <= -soundSpeedInt) {
        // supersonic inflow: boundary state fully determined by the freestream
        FloatType energyInf = _fluid.computeStaticEnergy_p_rho(pressureInf, densityInf);
        FloatType totEnergyInf = energyInf + 0.5 * velocityInf.dot(velocityInf);
        primitiveBoundary = StateVector({
            densityInf,
            velocityInf.x(),
            velocityInf.y(),
            velocityInf.z(),
            totEnergyInf});
    }
    else {
        // subsonic: combine the outgoing (interior) and incoming (freestream)
        // Riemann invariants
        FloatType riemannPlus  = normalVelInt + 2.0 * soundSpeedInt / (gamma - 1.0);
        FloatType riemannMinus = normalVelInf - 2.0 * soundSpeedInf / (gamma - 1.0);

        FloatType normalVelBound = 0.5 * (riemannPlus + riemannMinus);
        FloatType soundSpeedBound = 0.25 * (gamma - 1.0) * (riemannPlus - riemannMinus);

        // entropy and tangential velocity come from the upwind side
        FloatType densityRef, pressureRef;
        Vector3D velocityRef;
        if (normalVelBound > 0.0) {
            // outflow: take them from the interior
            densityRef = densityInt;
            pressureRef = pressureInt;
            velocityRef = velocityInt;
        }
        else {
            // inflow: take them from the freestream
            densityRef = densityInf;
            pressureRef = pressureInf;
            velocityRef = velocityInf;
        }

        FloatType entropy = pressureRef / std::pow(densityRef, gamma);
        FloatType densityBound = std::pow(
            soundSpeedBound * soundSpeedBound / (gamma * entropy), 1.0 / (gamma - 1.0));
        FloatType pressureBound = densityBound * soundSpeedBound * soundSpeedBound / gamma;

        // rebuild the velocity from the upwind tangential part and the
        // characteristic normal component
        Vector3D velocityTangential = velocityRef - normal * velocityRef.dot(normal);
        Vector3D velocityBound = velocityTangential + normal * normalVelBound;

        FloatType energyBound = _fluid.computeStaticEnergy_p_rho(pressureBound, densityBound);
        FloatType totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);
        primitiveBoundary = StateVector({
            densityBound,
            velocityBound.x(),
            velocityBound.y(),
            velocityBound.z(),
            totEnergyBound});
    }

    StateVector flux = computeAdvectionFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}
