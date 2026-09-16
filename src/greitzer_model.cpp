#include "greitzer_model.hpp"

GreitzerModel::GreitzerModel(
    const Config &config, 
    const FluidBase &fluid,
    FloatType throttleCoefficient)
    : _fluid(fluid), _config(config), _throttleCoefficient(throttleCoefficient)
{
    _plenumVolume = _config.getGreitzerPlenumVolume();
    _deltaTime = _config.getFixedTimeStep();
    _fluidGamma = _fluid.getGamma();
    _fluidRConstant = _fluid.getRconstant();
    _ambientPressure = _config.getGreitzerAmbientPressure();
}

void GreitzerModel::initializeState(
    FloatType plenumPressure, 
    FloatType plenumInletMassflow, 
    FloatType plenumOutletMassflow,
    FloatType plenumTemperature) {
    
    _plenumPressure.push_back(plenumPressure);
    _plenumInletMassflow.push_back(plenumInletMassflow);
    _plenumOutletMassflow.push_back(plenumOutletMassflow);
    _plenumTemperature = plenumTemperature;
    _time.push_back(0.0);
}



FloatType GreitzerModel::computePlenumPressure(FloatType massFlow) {
    _plenumInletMassflow.push_back(massFlow);
    _time.push_back(_time.back() + _deltaTime);
    
    FloatType aPlenum = std::sqrt(_fluidGamma * _fluidRConstant * _plenumTemperature); 

    FloatType newP = _plenumPressure.back() + _deltaTime * aPlenum * aPlenum / _plenumVolume * (
        _plenumInletMassflow.back() - _plenumOutletMassflow.back());
    _plenumPressure.push_back(newP);

    FloatType deltaP = _plenumPressure.back() - _ambientPressure;
    FloatType newM = (deltaP > 0.0) ? std::sqrt(deltaP / _throttleCoefficient) : 0.0;
    _plenumOutletMassflow.push_back(newM);
    
    return _plenumPressure.back();
}