#pragma once
#include <string>
#include <map>
#include <sstream>
#include "types.hpp"

class Config {

private:

    size_t massFlowUpdateFrequency = {100};

    std::map<std::string, std::string> _data;
    
    bool load(const std::string& filename);
    
    bool parseBool(const std::string& value, bool defaultValue) const;
    

    /**
     * @brief Return the string file corresponding to a certain entry.
     * @param value option to read from config file.
     * @param defaultNoneValue if true return a string 'None' if the option is not present. Otherwise throws an error
     */
    std::string parseString(const std::string& value, bool defaultNoneValue=false) const;
    
    FloatType parseFloat(const std::string& value) const;

    FloatType parseFloat(const std::string& value, FloatType defaultValue) const;

    int parseInt(const std::string& value) const;

    int parseInt(const std::string& value, int defaultValue) const;

    template <typename T>
    std::vector<T> parseVector(const std::string& key) const {
        if (!has(key)) {
            throw std::runtime_error("Missing value for key \"" + key + "\" in configuration.");
        }

        std::string value = get(key);
        std::vector<T> result;
        std::stringstream ss(value);
        std::string item;

        while (std::getline(ss, item, ',')) {
            // Trim leading/trailing whitespace
            item.erase(0, item.find_first_not_of(" \t\r\n"));
            item.erase(item.find_last_not_of(" \t\r\n") + 1);

            if (item.empty()) continue;

            std::stringstream itemStream(item);
            T val;
            itemStream >> val;

            if (itemStream.fail() || !itemStream.eof()) {
                throw std::runtime_error("Invalid vector element for key \"" + key + "\": \"" + item + "\"");
            }

            result.push_back(val);
        }

        return result;
    }

    template <typename T>
    std::vector<T> parseVector(const std::string& key, std::vector<T> defaultValue) const {
        if (!has(key)) {
            return defaultValue;
        }

        std::string value = get(key);
        std::vector<T> result;
        std::stringstream ss(value);
        std::string item;

        while (std::getline(ss, item, ',')) {
            // Trim leading/trailing whitespace
            item.erase(0, item.find_first_not_of(" \t\r\n"));
            item.erase(item.find_last_not_of(" \t\r\n") + 1);

            if (item.empty()) continue;

            std::stringstream itemStream(item);
            T val;
            itemStream >> val;

            if (itemStream.fail() || !itemStream.eof()) {
                throw std::runtime_error("Invalid vector element for key \"" + key + "\": \"" + item + "\"");
            }

            result.push_back(val);
        }

        return result;
    }

public:
    Config () = default;

    Config (const std::string& filename);
    
    ~Config () = default;
    
    std::string get(const std::string& key) const;
    
    bool has(const std::string& key) const;
    
    std::string gridFilePath() const;

    bool restartSolution() const;
    
    std::string restartSolutionFilepath() const;
    
    bool isBFMActive() const;

    bool isTurbulenceActive() const;
    
    bool isBlockageActive() const;
    
    BodyForceModel getBFMModel() const;

    TurbulenceModel getTurbulenceModel() const;

    KindSolver getKindSolver() const;

    Topology getTopology() const;

    FloatType getFluidGamma() const {return parseFloat("FLUID_GAMMA");} 

    FloatType getFluidGasConstant() const {return parseFloat("FLUID_R_CONSTANT");} 
    
    FloatType getFluidKinematicViscosity() const {return parseFloat("FLUID_KINEMATIC_VISCOSITY");} 

    FloatType getFluidPrandtlNumber() const {return parseFloat("FLUID_PRANDTL_NUMBER");}

    BoundaryType getBoundaryType(const BoundaryIndex) const;

    std::vector<FloatType> getInletBCValues() const { return parseVector<FloatType>("INLET_VALUE");}

    std::vector<FloatType> getOutletBCValues() const { return parseVector<FloatType>("OUTLET_VALUE");}

    FloatType getInitMachNumber() const { return parseFloat("INIT_MACH_NUMBER");}
    
    FloatType getInitPressure() const { return parseFloat("INIT_PRESSURE");}
    
    FloatType getInitTemperature() const { return parseFloat("INIT_TEMPERATURE");}

    ReferenceFrame getInletReferenceFrame() const;
    
    Vector3D getInitDirection() const ;

    bool isGongFormulationActive() const {return parseBool("GONG_MODELING_ACTIVE", false);}
    
    FloatType getCFL() const {return parseFloat("CFL");}

    long unsigned int getMaxIterations() const {return static_cast<long unsigned int>(parseFloat("N_ITERATIONS"));}

    bool saveUnsteadySolution() const {return parseBool("SAVE_UNSTEADY", true);}

    bool isSimulationSteady() const {return parseBool("SIMULATION_IS_STEADY", true);}

    int getSaveIterationsInterval() const {return static_cast<int>(parseFloat("SAVE_ITERATIONS_INTERVAL"));}

    FloatType getRotationalSpeedScalingFactor() const {return parseFloat("ROTATIONAL_SPEED_SCALING_FACTOR", 1.0);};

    FloatType getRotationalSpeedScalingFactor(FloatType currentTime) const;

    std::string getSolutionName() const {return parseString("SOLUTION_NAME");}

    std::string getInlet2DfilePath() const {return parseString("INLET_2D_FILEPATH");}

    std::string getChimaScalingFunctionsFile() const {return parseString("CHIMA_SCALING_FUNCTIONS_FILEPATH");}

    Vector3D getNoSlipWallVelocity(BoundaryIndex boundary) const;

    FloatType getChimaReferenceMassFlow() const {return parseFloat("CHIMA_REFERENCE_MASS_FLOW");}

    FloatType getBfmRelaxationFactor() const;

    std::string getRestartFilepath() const {return parseString("RESTART_SOLUTION_FILEPATH");}

    std::string getRestartType() const {return parseString("RESTART_TYPE", true);}

    bool saveTurboOutput() const {return parseBool("TURBO_OUTPUT", false);}

    bool isBfmLagActive() const {return parseBool("BFM_LAG_ACTIVE", false);}

    bool isMonitorPointsActive() const {return parseBool("MONITOR_POINTS_ACTIVE", false);}

    bool saveMeshQualityStats() const {return parseBool("SAVE_MESH_QUALITY_STATS", false);}

    size_t getResidualsDropConvergence() const {return static_cast<size_t>(parseInt("RESIDUALS_DROP_CONVERGENCE", 16));}

    size_t getTrailingEdgeIndex() const {return static_cast<size_t>(parseInt("TRAILING_EDGE_INDEX"));}
    
    size_t getLeadingEdgeIndex() const {return static_cast<size_t>(parseInt("LEADING_EDGE_INDEX"));}

    TimeIntegration getTimeIntegration() const;

    TimeStepMethod getTimeStepMethod() const;

    FluidModel getFluidModel() const;

    AdvectionScheme getAdvectionScheme() const;

    // get the coefficients for the Runge-Kutta time integration. (Taken from Simon thesis, page 64)
    std::vector<FloatType> getTimeIntegrationCoeffs() const ;

    size_t getMassFlowUpdateFrequency() const {return massFlowUpdateFrequency;}

    size_t getSolutionOutputFrequency() const {return getSaveIterationsInterval();}

    FloatType getPeriodicityAngleRad() const {return getPeriodicityAngleDeg()*M_PI/180.0;} 

    FloatType getPeriodicityAngleDeg() const {return getPeriodicityInfo()[1];} 

    FloatType getPeriodicityTranslation() const {return getPeriodicityInfo()[0];} 

    std::vector<FloatType> getPeriodicityInfo() const;

    void printAllConfigValues() const;

    FloatType computeRampedOutletPressure(const size_t iterCounter, const FloatType outletPressure) const;

    FloatType getOutletPressureRampMaxIterations() const {return parseFloat("OUTLET_PRESSURE_RAMP_ITERATIONS", 0.0);}

    FloatType getHallTholletCoefficient_KN() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KN");} 

    FloatType getHallTholletCoefficient_KF() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KF");} 

    FloatType getHallTholletCoefficient_KD() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KD");} 

    bool getHallTholletOffDesignActive() const;

    int getMonitorPointsCoordsI() const { return parseInt("MONITOR_POINTS_I_COORDS");}

    int getMonitorPointsCoordsJ() const { return parseInt("MONITOR_POINTS_J_COORDS");}

    int getCircumferentialNumberMonitorPoints() const { return parseInt("MONITOR_POINTS_NUMBER");}

    bool isPerturbationBodyForceActive() const {return parseBool("PERTURBATION_BODY_FORCE", false);}

    bool getMUSCLreconstruction() const {return parseBool("MUSCL_RECONSTRUCTION", false);}

    bool isStalledBfmActive() const {return parseBool("STALLED_BFM_ACTIVE", false);}

    bool enableGreitzerModeling() const;

    FluxLimiter getFluxLimiter() const ;

    std::vector<FloatType> getPerturbationCenter() const {return parseVector<FloatType>("PERTURBATION_CENTER");}

    FloatType getPerturbationScalingFactor() const {return parseFloat("PERTURBATION_SCALING_FACTOR");} 

    FloatType getPerturbationTimeDuration() const {return parseFloat("PERTURBATION_TIME_DURATION");} 
    
    FloatType getPerturbationStartTime() const {return parseFloat("PERTURBATION_TIME_START");} 

    FloatType getPerturbationRadialExtension() const {return parseFloat("PERTURBATION_RADIAL_EXTENSION");} 

    FloatType getFixedTimeStep() const {return parseFloat("FIXED_TIME_STEP");} 

    OutputFields getOutputFields() const ;

    bool isViscosityActive() const {return parseBool("VISCOSITY_ACTIVE", false);}

    FloatType getFluidHeatCapacity() const {return parseFloat("FLUID_CP");}

    FloatType getKnCorrelationBfmCoefficient() const {return parseFloat("KN_CORRELATION_BFM_COEFFICIENT");}

    FloatType getGreitzerPlenumVolume() const {return parseFloat("GREITZER_PENUM_VOLUME");} 

    FloatType getGreitzerThrottleCoefficient() const {return getOutletBCValues()[0];} 

    bool isAccelerationActive() const {return parseBool("ACCELERATION_ACTIVE", false);}

    FloatType getAccelerationInitialTime() const {return parseFloat("ACCELERATION_INITIAL_TIME");}

    FloatType getAccelerationFinalTime() const {return parseFloat("ACCELERATION_FINAL_TIME");}

    FloatType getRotationalSpeedScalingFactorInitial() const {
        return parseFloat("ROTATIONAL_SPEED_SCALING_FACTOR_INITIAL");
    }

    FloatType getRotationalSpeedScalingFactorFinal() const {
        return parseFloat("ROTATIONAL_SPEED_SCALING_FACTOR_FINAL");
    }
};
