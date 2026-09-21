#include "advection_roe.hpp"
#include "math_utils.hpp"



StateVector AdvectionRoe::computeFlux(
    const StateVector& Ull,
    const StateVector& Ul,
    const StateVector& Ur,
    const StateVector& Urr,
    const Vector3D& S) const
{
    StateVector Wl = getPrimitiveVariablesFromConservative(Ul);
    StateVector Wr = getPrimitiveVariablesFromConservative(Ur);
    StateVector Wrr = getPrimitiveVariablesFromConservative(Urr);
    StateVector Wll = getPrimitiveVariablesFromConservative(Ull);

    if (_isMusclActive){
        musclReconstructLeftRight(Wll, Wl, Wr, Wrr, _fluxLimiter);
    }

    RoeState state;

    computeNormalTriad(S, state);

    StateVector WnormL = computeRotatedPrimitive(Wl, state);
    StateVector WnormR = computeRotatedPrimitive(Wr, state);

    computeRoeAvgVariables(WnormL, WnormR, state);
    computeEigenvalues(state);
    computeEigenvectors(state);
    computeWaveStrengths(WnormL, WnormR, state);

    StateVector flux({0.0, 0.0, 0.0, 0.0, 0.0});
    assembleTotalFlux(S, WnormL, WnormR, flux, state);

    return flux;
}


void AdvectionRoe::computeNormalTriad(const Vector3D& S, RoeState& state) const {
    
    // the first is the surface normal direction
    state.n1 = S.normalized();
    
    // the second is simply taken orthogonal to the first
    if (state.n1.x() == 0 && state.n1.y() == 0){
        state.n2 = Vector3D(1, 0, 0);
    }
    else{
        state.n2 = Vector3D(-state.n1.y(), state.n1.x(), 0);
    }
    state.n2 = state.n2.normalized();
    
    // the third is the cross product
    state.n3 = state.n1.cross(state.n2);
    
}


StateVector AdvectionRoe::computeRotatedPrimitive(const StateVector& W, const RoeState& state) const {
    Vector3D velocity ({W[1], W[2], W[3]});
    StateVector rotatedState ({W[0], velocity.dot(state.n1), velocity.dot(state.n2), velocity.dot(state.n3), W[4]});
    return rotatedState;
}


void AdvectionRoe::computeRoeAvgVariables(const StateVector& Wl, const StateVector& Wr, RoeState& state) const {

    FloatType rhoL = Wl[0];
    FloatType rhoR = Wr[0];
    FloatType u1L = Wl[1];
    FloatType u1R = Wr[1];
    FloatType u2L = Wl[2];
    FloatType u2R = Wr[2];
    FloatType u3L = Wl[3];
    FloatType u3R = Wr[3];
    FloatType etL = Wl[4];
    FloatType etR = Wr[4];
    FloatType htL = _fluid.computeTotalEnthalpy_rho_u_et(rhoL, {u1L, u2L, u3L}, etL);
    FloatType htR = _fluid.computeTotalEnthalpy_rho_u_et(rhoR, {u1R, u2R, u3R}, etR);

    state.rhoAVG = std::sqrt(rhoL * rhoR);
    state.u1AVG = roeAverage(rhoL, rhoR, u1L, u1R);
    state.u2AVG = roeAverage(rhoL, rhoR, u2L, u2R);
    state.u3AVG = roeAverage(rhoL, rhoR, u3L, u3R);
    state.htAVG = roeAverage(rhoL, rhoR, htL, htR);

    FloatType eL = etL - 0.5 * (u1L * u1L + u2L * u2L + u3L * u3L);
    FloatType eR = etR - 0.5 * (u1R * u1R + u2R * u2R + u3R * u3R);
    FloatType pL = _fluid.computePressure_rho_e(rhoL, eL);
    FloatType pR = _fluid.computePressure_rho_e(rhoR, eR);

    FloatType chiL = _fluid.computeDpDrho_e(rhoL, eL);
    FloatType kappaL = _fluid.computeDpDe_rho(rhoL, eL);
    FloatType chiR = _fluid.computeDpDrho_e(rhoR, eR);
    FloatType kappaR = _fluid.computeDpDe_rho(rhoR, eR);

    FloatType barChi = 0.5 * (chiL + chiR);
    FloatType barKappa = 0.5 * (kappaL + kappaR);
    FloatType barRho = 0.5 * (rhoL + rhoR);
    FloatType barE = 0.5 * (eL + eR);

    FloatType deltaP = pR - pL;
    FloatType deltaRho = rhoR - rhoL;
    FloatType deltaE = eR - eL;

    // Vinokur-Montagne orthogonal projection to satisfy jump condition: deltaP = chi * deltaRho + kappa * deltaE
    FloatType dRhoHat = (barRho > 1e-12) ? deltaRho / barRho : 0.0;
    FloatType dEHat = (std::abs(barE) > 1e-12) ? deltaE / std::abs(barE) : 0.0;
    FloatType dNorm2 = dRhoHat * dRhoHat + dEHat * dEHat;

    FloatType chiTilde = barChi;
    FloatType kappaTilde = barKappa;

    if (dNorm2 > 1e-12) {
        FloatType residual = deltaP - (barChi * deltaRho + barKappa * deltaE);
        chiTilde += (residual * dRhoHat) / (dNorm2 * barRho);
        kappaTilde += (residual * dEHat) / (dNorm2 * std::abs(barE));
    }

    if (barKappa > 0.0 && kappaTilde <= 0.0) {
        kappaTilde = 0.01 * barKappa;
    }

    state.chiAVG = chiTilde;
    state.kappaAVG = kappaTilde;
    state.kappaPrime = (state.rhoAVG > 1e-12) ? kappaTilde / state.rhoAVG : 0.4;

    FloatType eAVG = roeAverage(rhoL, rhoR, eL, eR);
    state.chiPrime = chiTilde - state.kappaPrime * eAVG;

    FloatType q2AVG = state.u1AVG * state.u1AVG + state.u2AVG * state.u2AVG + state.u3AVG * state.u3AVG;
    FloatType hAVG = state.htAVG - 0.5 * q2AVG;

    FloatType a2 = state.chiPrime + state.kappaPrime * hAVG;
    if (a2 <= 0.0) {
        FloatType aL = _fluid.computeSoundSpeed_rho_e(rhoL, eL);
        FloatType aR = _fluid.computeSoundSpeed_rho_e(rhoR, eR);
        a2 = 0.5 * (aL * aL + aR * aR);
    }
    state.aAVG = std::sqrt(a2);
}


FloatType AdvectionRoe::roeAverage(FloatType rhoL, FloatType rhoR, FloatType phiL, FloatType phiR) const {
    FloatType avg = (std::sqrt(rhoL) * phiL + std::sqrt(rhoR) * phiR) / (std::sqrt(rhoL) + std::sqrt(rhoR));
    return avg;
}


void AdvectionRoe::computeEigenvalues(RoeState& state) const {
    state.eigenvalues[0] = state.u1AVG - state.aAVG;
    state.eigenvalues[1] = state.u1AVG;
    state.eigenvalues[2] = state.u1AVG;
    state.eigenvalues[3] = state.u1AVG;
    state.eigenvalues[4] = state.u1AVG + state.aAVG;
}

void AdvectionRoe::computeEigenvectors(RoeState& state) const {
    FloatType contactE5 = state.htAVG - (state.aAVG * state.aAVG) / state.kappaPrime;
    state.eigenvectors[0] = StateVector({1.0, state.u1AVG - state.aAVG, state.u2AVG, state.u3AVG, state.htAVG - state.aAVG * state.u1AVG});
    state.eigenvectors[1] = StateVector({1.0, state.u1AVG, state.u2AVG, state.u3AVG, contactE5});
    state.eigenvectors[2] = StateVector({0.0, 0.0, 1.0, 0.0, state.u2AVG});
    state.eigenvectors[3] = StateVector({0.0, 0.0, 0.0, 1.0, state.u3AVG});
    state.eigenvectors[4] = StateVector({1.0, state.u1AVG + state.aAVG, state.u2AVG, state.u3AVG, state.htAVG + state.aAVG * state.u1AVG});
}

void AdvectionRoe::computeWaveStrengths(const StateVector& WnormL, const StateVector& WnormR, RoeState& state) const {
    FloatType deltaRho = WnormR[0] - WnormL[0];
    FloatType deltaU1 = WnormR[1] - WnormL[1];
    FloatType deltaU2 = WnormR[2] - WnormL[2];
    FloatType deltaU3 = WnormR[3] - WnormL[3];
    FloatType pL = _fluid.computePressure_rho_u_et(WnormL[0], {WnormL[1], WnormL[2], WnormL[3]}, WnormL[4]);
    FloatType pR = _fluid.computePressure_rho_u_et(WnormR[0], {WnormR[1], WnormR[2], WnormR[3]}, WnormR[4]);
    FloatType deltaP = pR - pL;
    
    state.waveStrengths[0] = 1.0 / (2.0 * state.aAVG * state.aAVG) * (deltaP - state.rhoAVG * state.aAVG * deltaU1);
    state.waveStrengths[1] = deltaRho - deltaP / (state.aAVG * state.aAVG);
    state.waveStrengths[2] = state.rhoAVG * deltaU2;
    state.waveStrengths[3] = state.rhoAVG * deltaU3;
    state.waveStrengths[4] = 1.0 / (2.0 * state.aAVG * state.aAVG) * (deltaP + state.rhoAVG * state.aAVG * deltaU1);
}


void AdvectionRoe::assembleTotalFlux(
    const Vector3D& S, 
    const StateVector& WnormL, 
    const StateVector& WnormR, 
    StateVector& flux,
    const RoeState& state) const {

    // Compute flux in normal reference frame (n1, n2, n3).
    // The direction is (1, 0, 0) because the state has been aligned to the normal triad
    StateVector fluxL = computeAdvectionFluxFromPrimitive(WnormL, {1.0, 0.0, 0.0}, _fluid);
    StateVector fluxR = computeAdvectionFluxFromPrimitive(WnormR, {1.0, 0.0, 0.0}, _fluid);

    // Entropy fix of the eigenvalues
    StateVector fixedEigenvalues({0.0, 0.0, 0.0, 0.0, 0.0});
    for (size_t i = 0; i < fixedEigenvalues.size(); ++i) {
        if (std::abs(state.eigenvalues[i]) < _entropyFixCoefficient) {
            fixedEigenvalues[i] = _entropyFixCoefficient;
        } else {
            fixedEigenvalues[i] = std::abs(state.eigenvalues[i]);
        }
    }

    StateVector fluxRoe = (fluxL + fluxR) * 0.5;
    for (size_t i = 0; i < fixedEigenvalues.size(); ++i) {
        for (size_t j = 0; j < fixedEigenvalues.size(); ++j) {
            fluxRoe[i] -= 0.5 * state.waveStrengths[j] * fixedEigenvalues[j] * state.eigenvectors[j][i];
        }
    }

    // project the flux back to the original cartesian reference frame (x,y,z)
    flux[0] = fluxRoe[0];
    flux[1] = fluxRoe[1] * (state.n1.dot(_xVersor)) + fluxRoe[2] * (state.n2.dot(_xVersor)) + fluxRoe[3] * (state.n3.dot(_xVersor));
    flux[2] = fluxRoe[1] * (state.n1.dot(_yVersor)) + fluxRoe[2] * (state.n2.dot(_yVersor)) + fluxRoe[3] * (state.n3.dot(_yVersor));
    flux[3] = fluxRoe[1] * (state.n1.dot(_zVersor)) + fluxRoe[2] * (state.n2.dot(_zVersor)) + fluxRoe[3] * (state.n3.dot(_zVersor));
    flux[4] = fluxRoe[4];
}