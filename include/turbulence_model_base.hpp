#pragma once

#include "types.hpp"
#include "mesh.hpp"
#include "config.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "math_utils.hpp"
#include "input_table.hpp"
#include "solver_euler.hpp"


class TurbulenceModelBase {
public:
    explicit TurbulenceModelBase() {};
    virtual ~TurbulenceModelBase() = default;

    // virtual void solve() = 0;
    // virtual void initialize() = 0;

protected:
    // SolverEuler* _solver;
};
