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
    
    _currentPlenumPressure = plenumPressure;
    _currentPlenumInletMassflow = plenumInletMassflow;
    _currentPlenumOutletMassflow = plenumOutletMassflow;
    _plenumTemperature = plenumTemperature;
    _currentTime = 0.0;

    _plenumPressure.push_back(_currentPlenumPressure);
    _plenumInletMassflow.push_back(_currentPlenumInletMassflow);
    _plenumOutletMassflow.push_back(_currentPlenumOutletMassflow);
    _time.push_back(_currentTime);
}



FloatType GreitzerModel::computePlenumPressure(FloatType massFlow) {
    _currentPlenumInletMassflow = massFlow;
    _currentTime += _deltaTime;
    
    FloatType aPlenum = std::sqrt(_fluidGamma * _fluidRConstant * _plenumTemperature); 

    FloatType newP = _currentPlenumPressure + _deltaTime * aPlenum * aPlenum / _plenumVolume * (
        _currentPlenumInletMassflow - _currentPlenumOutletMassflow);
    _currentPlenumPressure = newP;

    FloatType deltaP = _currentPlenumPressure - _ambientPressure;
    FloatType newM = (deltaP > 0.0) ? std::sqrt(deltaP / _throttleCoefficient) : 0.0;
    _currentPlenumOutletMassflow = newM;

    _plenumInletMassflow.push_back(_currentPlenumInletMassflow);
    _time.push_back(_currentTime);
    _plenumPressure.push_back(_currentPlenumPressure);
    _plenumOutletMassflow.push_back(_currentPlenumOutletMassflow);
    
    return _currentPlenumPressure;
}

void GreitzerModel::clearBuffer() {
    _plenumInletMassflow.clear();
    _plenumOutletMassflow.clear();
    _plenumPressure.clear();
    _time.clear();
}