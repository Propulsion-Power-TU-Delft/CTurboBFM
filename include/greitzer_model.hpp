#pragma once

#include "types.hpp"
#include "fluid_base.hpp"
#include "config.hpp"

class GreitzerModel {
public:
    GreitzerModel(
        const Config &config, 
        const FluidBase &fluid,
        FloatType throttleCoefficient);

    virtual ~GreitzerModel() = default;  

    FloatType computePlenumPressure(FloatType massFlow);

    void initializeState(FloatType plenumPressure, FloatType plenumInletMassflow, FloatType plenumOutletMassflow);

    FloatType getTime(size_t i) const { return _time.at(i); }

    FloatType getPlenumPressure(size_t i) const { return _plenumPressure.at(i); }
    
    FloatType getPlenumInletMassflow(size_t i) const { return _plenumInletMassflow.at(i); }
    
    FloatType getPlenumOutletMassflow(size_t i) const { return _plenumOutletMassflow.at(i); }
    
    size_t getSize() const { return _time.size(); }

private:
    const FluidBase& _fluid;  
    const Config& _config;
    std::vector<FloatType> _plenumInletMassflow;
    std::vector<FloatType> _plenumOutletMassflow;
    std::vector<FloatType> _plenumPressure;
    std::vector<FloatType> _time;
    FloatType _plenumVolume; 
    FloatType _throttleCoefficient;
    FloatType _deltaTime = 0.0;
    FloatType _fluidGamma = 1.4; 
    FloatType _fluidRConstant = 287.0;
};
