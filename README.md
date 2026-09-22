# CTurboBFM #
BFM solver for turbomachinery flows. It couples a 3D finite volume Euler/RANS solver with body force models needed to mimic blade effects on the flow field.

### What is this repository for? ###

* Calculation of fluid dynamics problems (1D, 2D, Axisymmetric, and 3D).
* Simulation of turbomachinery flows through body force models (axisymmetric and full-annulus grids).
* Ideal and real gas thermodynamics (lookup tables generated via CoolProp).
* Laminar and turbulent viscous modeling (Spalart-Allmaras).
* Dynamic compression system surge/stall modeling via coupled Greitzer lumped-parameter models.
* The Python [Unsflow](https://github.com/Propulsion-Power-TU-Delft/unsflow) package can be used to generate grid files.

---

### Requirements
- C++17 compiler (GCC / Clang)
- Python 3.x (with `numpy` and `CoolProp` for real-gas table generation)
- GoogleTest (for unit tests)
- Make
- [Unsflow](https://github.com/Propulsion-Power-TU-Delft/unsflow)

---

### Build system ###
Based on `Makefile`, tested for macOS and Linux. Commands:
```bash
make release   # Optimized build (-O3)
make debug     # Debug build (-g -O0)
make clean     # Remove build artifacts
```

---

### Running a simulation

```bash
turbobfm input.ini
```

---

### Project structure

- `grids/` → Mesh files and generation scripts
- `include/` → C++ header files (`.hpp`)
- `python/` → Post-processing and thermodynamic table generation scripts
- `src/` → Solver, physics, numerical schemes, and boundary conditions
- `test/` → Unit tests and regression test cases
- `testcases/` → Ready-to-run example configurations

---

## Configuration Reference (`input.ini`)

CTurboBFM simulations are configured using a standard INI file (case-insensitive keys). Options are grouped logically below.

### 1. Mesh & Topology

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `GRID_FILE` | String | *Required* | Path to the grid coordinates file (`.csv`). |
| `TOPOLOGY` | String | `1D`, `2D`, `3D`, `axisymmetric` | Spatial topology of the computational domain. |
| `SAVE_MESH_QUALITY_STATS` | Boolean | `no` | If `yes`, computes and exports mesh quality metrics (skewness, non-orthogonality, aspect ratios, cell volumes) to CSV. |

---

### 2. Solver & Time-Stepping Controls

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `KIND_SOLVER` | String | `Euler`, `RANS` | Governing equations to solve: inviscid (`Euler`) or viscous/turbulent (`RANS`). |
| `CFL` | Float | *Required* | Courant-Friedrichs-Lewy (CFL) number. |
| `N_ITERATIONS` | Integer | *Required* | Maximum number of solver iterations. |
| `SIMULATION_IS_STEADY` | Boolean | `yes` | When `yes`, monitors residual drop and can terminate early upon convergence. |
| `RESIDUALS_DROP_CONVERGENCE` | Integer | `16` | Orders of magnitude drop in density residual required to consider the steady solution converged. |
| `TIME_STEP_METHOD` | String | `local`, `global`, `fixed` | Time-stepping algorithm:<br>• `local`: Cell-local $\Delta t$ for accelerated steady-state convergence.<br>• `global`: Uniform $\Delta t = \min(\Delta t_{\text{local}})$ for coherent CFL-limited time-marching.<br>• `fixed`: Constant physical time step specified by `FIXED_TIME_STEP`. |
| `FIXED_TIME_STEP` | Float | *Required if `fixed`* | Fixed time step size in seconds ($\Delta t$). |
| `TIME_INTEGRATION_TYPE` | Integer | `0`, `1` | Runge-Kutta explicit scheme:<br>• `0`: 4-stage RK4 (RK41A, higher boundary dissipation, max CFL $\approx 2.58$).<br>• `1`: 3-stage RK3 (RK31B, good transient damping, max CFL $\approx 1.53$). |

---

### 3. Spatial Discretization & Advection Schemes

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `CONVECTION_SCHEME` | String | `JST`, `Roe` | Advection numerical flux scheme: Jameson-Schmidt-Turkel (`JST`) or approximate Riemann solver (`Roe`). |
| `MUSCL_RECONSTRUCTION` | Boolean | `no` | Enables 2nd-order spatial accuracy for the `Roe` scheme via MUSCL primitive reconstruction. |
| `FLUX_LIMITER` | String | `None`, `van albada`, `van leer`, `min mod` | Slope limiter used in MUSCL reconstruction (default: `None`). |

---

### 4. Fluid Thermodynamics

#### Ideal Gas (`FLUID_MODEL = Ideal`)

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `FLUID_MODEL` | String | `Ideal` | Selects ideal gas equation of state. |
| `FLUID_NAME` | String | `Air` | Fluid substance name. |
| `FLUID_GAMMA` | Float | *Required for Ideal* | Ratio of specific heats $\gamma = c_p / c_v$ (e.g. `1.4` for air). |
| `FLUID_R_CONSTANT` | Float | *Required for Ideal* | Specific gas constant $R$ [$\text{J}/(\text{kg}\cdot\text{K})$] (e.g. `287.05` for air). |
| `FLUID_CP` | Float | Optional | Specific heat capacity at constant pressure $c_p$ [$\text{J}/(\text{kg}\cdot\text{K})$]. |
| `FLUID_PRANDTL_NUMBER` | Float | Optional | Laminar Prandtl number $Pr$ (e.g. `0.708`). |

#### Real Gas Look-Up Table (`FLUID_MODEL = Real`)

Real gas calculations use structured Look-Up Tables (`.lut`) generated via CoolProp:

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `FLUID_MODEL` | String | `Real` | Selects real gas look-up table (LUT) model. |
| `FLUID_NAME` (or `FLUID`) | String | `CO2` | Fluid name recognized by CoolProp (e.g. `CO2`, `Air`, `Nitrogen`). |
| `FLUID_TABLE_FILE` | String | `fluid_table_<FLUID>.lut` | Path to the thermodynamic look-up table file. |
| `FLUID_TABLE_FORCE_REGENERATE` | Boolean | `no` | Force re-generation of the table via CoolProp even if the file exists. |
| `FLUID_TABLE_P_MIN`, `FLUID_TABLE_P_MAX` | Float | $10^5$, $10^7$ Pa | Pressure range $[p_{\min}, p_{\max}]$ [Pa] for the table (can also default based on `INIT_PRESSURE`). |
| `FLUID_TABLE_T_MIN`, `FLUID_TABLE_T_MAX` | Float | 200.0, 600.0 K | Temperature range $[T_{\min}, T_{\max}]$ [K] for the table. |
| `FLUID_TABLE_N_RHO`, `FLUID_TABLE_N_E` | Integer | `60`, `60` | Number of grid points for density $\rho$ and internal energy $e$. |
| `FLUID_TABLE_N_P`, `FLUID_TABLE_N_T`, `FLUID_TABLE_N_S` | Integer | `40`, `40`, `40` | Inversion grid resolution for $(p, T)$ and $(p, s)$ tables. |
| `FLUID_TABLE_PHASE` | String | `""` | Phase restriction passed to CoolProp (e.g. `gas`). |
| `FLUID_TABLE_PYTHON` | String | `python3` | Python binary used to run the table generator script. |
| `FLUID_TABLE_GENERATOR` | String | Auto-detected | Path to `python/generate_fluid_table.py`. |
| `STOP_ON_TABLE_OUT_OF_BOUNDS` | Boolean | `no` | Terminate the simulation immediately if thermodynamic state leaves the table bounds. |

---

### 5. Viscosity & Turbulence Models

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `VISCOSITY_ACTIVE` | Boolean | `no` | Enables molecular viscous and diffusive fluxes. |
| `FLUID_VISCOSITY_MODEL` | String | `Sutherland`, `Constant` | Viscosity model (default: `Sutherland`). |
| `FLUID_KINEMATIC_VISCOSITY` | Float | `1e-5` | Constant kinematic viscosity $\nu$ [$\text{m}^2/\text{s}$]. |
| `FLUID_MU_CONSTANT` | Float | *Required if `Constant`* | Dynamic viscosity $\mu$ [$\text{Pa}\cdot\text{s}$] (alias: `FLUID_CONSTANT_MU`). |
| `FLUID_SUTHERLAND_MU_REF` | Float | *Required if `Sutherland`* | Sutherland reference viscosity $\mu_{\text{ref}}$ [$\text{Pa}\cdot\text{s}$]. |
| `FLUID_SUTHERLAND_T_REF` | Float | *Required if `Sutherland`* | Sutherland reference temperature $T_{\text{ref}}$ [$\text{K}$]. |
| `FLUID_SUTHERLAND_S` | Float | *Required if `Sutherland`* | Sutherland temperature constant $S$ [$\text{K}$]. |
| `TURBULENCE_ACTIVE` | Boolean | `no` | Enables turbulence closure modeling. |
| `TURBULENCE_MODEL` | String | `SA`, `None` | Turbulence model (`SA`: Spalart-Allmaras one-equation model). |
| `TURBULENT_PRANDTL_NUMBER` | Float | `0.86` | Turbulent Prandtl number $Pr_t$. |

---

### 6. Flow Initialization

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `INIT_MACH_NUMBER` | Float | *Required* | Initial freestream Mach number throughout the domain. |
| `INIT_PRESSURE` | Float | *Required* | Initial static pressure [$\text{Pa}$]. |
| `INIT_TEMPERATURE` | Float | *Required* | Initial static temperature [$\text{K}$]. |
| `INIT_DIRECTION` | 3 Floats / String | *Required* | Flow direction unit vector (`x, y, z`), or `adaptive` to follow grid channel geometry automatically. |
| `INLET_REFERENCE_FRAME` | String | `Cartesian`, `Cylindrical` | Coordinate frame for inlet velocity vectors (default: `Cartesian`). |

---

### 7. Solution Restart

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `RESTART_SOLUTION` | Boolean | `no` | Restart simulation from a previous solution CSV file. |
| `RESTART_SOLUTION_FILEPATH` | String | *Required if restart* | Path to the restart file (`.csv`). |
| `RESTART_TYPE` | String | `standard`, `axisymmetric`, `nearest_neighbor`, `nearest_neighbor_axisymmetric` | Interpolation method:<br>• `standard`: 1-to-1 node copy (mesh dimensions must match).<br>• `axisymmetric`: Extrudes a 2D axisymmetric meridional solution onto a 3D sector or full annulus.<br>• `nearest_neighbor`: 3D KD-Tree spatial interpolation.<br>• `nearest_neighbor_axisymmetric`: 2D $(x, r)$ KD-Tree interpolation. |

---

### 8. Boundary Conditions

Boundary patches are specified using:
1. `BOUNDARY_CONDITIONS_FILEPATH = <path_to_boundaries.csv>`
2. An entry in `input.ini` for each named patch matching the format:
   ```ini
   <PATCH_NAME> = <type>, <param1>, <param2>, ...
   ```

#### Supported Boundary Types

| Boundary Type | Syntax in `input.ini` | Parameters & Description |
| :--- | :--- | :--- |
| **Subsonic Inlet** | `<NAME> = inlet, Pt, Tt, nx, ny, nz` | Total pressure $P_t$ [Pa], total temperature $T_t$ [K], unit flow direction vector $(n_x, n_y, n_z)$.<br>*Note: Using `1, 0, 0` enforces axial direction. Direction `1, 1, 1` enforces flow normal to the face.* |
| **2D Non-Uniform Inlet** | `<NAME> = inlet_2d` | Reads spatial profile from `INLET_2D_FILEPATH = <path.csv>`. |
| **Supersonic Inlet** | `<NAME> = inlet_supersonic, p, T, u, v, w` | Static pressure $p$ [Pa], static temperature $T$ [K], and Cartesian velocity components $u, v, w$ [m/s]. |
| **Subsonic Outlet** | `<NAME> = outlet, p` | Static back-pressure $p$ [Pa]. |
| **Supersonic Outlet** | `<NAME> = outlet_supersonic` | Extrapolates all state variables from domain interior. |
| **Radial Equilibrium** | `<NAME> = radial_equilibrium, p_hub` | Hub static pressure $p_{\text{hub}}$ [Pa]. Solves $\frac{\partial p}{\partial r} = \frac{\rho v_\theta^2}{r}$ radially along $j$ on outlet face. |
| **Throttle Outlet** | `<NAME> = throttle, kt` | Throttle characteristic: adjusts $p_{\text{hub}} = P_{t,\text{in}} + k_t \dot{m}^2$ dynamically for compressor map generation. |
| **Inviscid Wall** | `<NAME> = inviscid_wall` | Slip wall condition (Euler wall: $\vec{v} \cdot \vec{n} = 0$). |
| **No-Slip Wall** | `<NAME> = no_slip_wall, [u, v, w]` | Viscous wall with optional velocity components. |
| **Periodic** | `<NAME> = periodic, id, translation, angle_deg` | Periodic boundary pair (`id` = 0 or 1), axial translation distance [m], and circumferential sector angle [deg] (e.g. `360.0`). |
| **Wedge** | `<NAME> = wedge` | Symmetry condition for axisymmetric 2.5D wedge grids. |
| **Farfield** | `<NAME> = farfield, M_inf, p_inf, T_inf, nx, ny, nz` | Non-reflecting characteristic farfield based on Riemann invariants. |
| **Transparent** | `<NAME> = transparent` | Open non-reflecting boundary using advection scheme extrapolation. |
| **Empty** | `<NAME> = empty` | Inactive boundary placeholder. |

#### Additional Boundary Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `OUTLET_PRESSURE_RAMP_ITERATIONS` | Float | `0.0` | Linearly ramps outlet back-pressure from `INIT_PRESSURE` to target pressure over $N$ iterations to prevent initial shockwaves. |
| `BOUNDARY_<I\|J\|K>_<START\|END>_VELOCITY` | 3 Floats | `0.0, 0.0, 0.0` | Prescribes moving wall velocity components `u, v, w` [m/s] for specific boundary faces (e.g. `BOUNDARY_J_END_VELOCITY = 100, 0, 0`). |
| `INLET_2D_FILEPATH` | String | Optional | Path to CSV file containing non-uniform 2D inlet conditions. |

---

### 9. Body Force Modeling (BFM)

| Parameter | Type | Default / Options | Description |
| :--- | :--- | :--- | :--- |
| `BFM_ACTIVE` | Boolean | `no` | Enables blade body force modeling source terms. |
| `BLOCKAGE_ACTIVE` | Boolean | `no` | Accounts for blade metal blockage and associated source terms. |
| `BFM_MODEL` | String | `Hall-Thollet`, `Hall`, `Chima`, `Blockage`, `Correlations`, `Lift-Drag`, `None` | BFM formulation to use. |
| `LEADING_EDGE_INDEX` | Integer | Optional | Grid $i$-index representing the blade leading edge. |
| `TRAILING_EDGE_INDEX` | Integer | Optional | Grid $i$-index representing the blade trailing edge. |

#### Hall / Hall-Thollet Model Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `HALL_THOLLET_COEFFICIENT_KN` | Float | `1.0` | Normal force calibration coefficient $K_n$. |
| `HALL_THOLLET_COEFFICIENT_KF` | Float | `1.0` | Parallel friction force calibration coefficient $K_f$. |
| `HALL_THOLLET_COEFFICIENT_KD` | Float | `1.0` | Off-design loss model exponent $K_d$. |
| `HALL_THOLLET_OFF_DESIGN_ACTIVE` | Boolean | `no` | Enables quadratic off-design incidence loss modeling. |
| `STALLED_BFM_ACTIVE` | Boolean | `no` | Deactivates viscous force if flow deviation angle exceeds the stall threshold. |

#### Chima Model Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `CHIMA_SCALING_FUNCTIONS_FILEPATH` | String | *Required for Chima* | CSV file containing turning and loss scaling factors as functions of mass flow. |
| `CHIMA_REFERENCE_MASS_FLOW` | Float | *Required for Chima* | Design reference mass flow rate [$\text{kg}/\text{s}$]. |
| `BFM_LAG_ACTIVE` | Boolean | `no` | Applies under-relaxation lag to body force evolution. |
| `BFM_RELAXATION_FACTOR` | Float | `0.1` | Relaxation factor used when `BFM_LAG_ACTIVE = yes`. |

#### Correlations BFM Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `KN_CORRELATION_BFM_COEFFICIENT` | Float | *Required for Correlations* | Calibration coefficient for correlation-based BFM. |

#### Rotational Speed & Dynamic Acceleration

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `ROTATIONAL_SPEED_SCALING_FACTOR` | Float | `1.0` | Multiplier for the mesh rotor blade speed $\Omega$. |
| `ACCELERATION_ACTIVE` | Boolean | `no` | Enables time-dependent rotor speed ramp (requires `TIME_STEP_METHOD = fixed`). |
| `ACCELERATION_INITIAL_TIME` | Float | *Required if active* | Physical start time [s] of acceleration ramp. |
| `ACCELERATION_FINAL_TIME` | Float | *Required if active* | Physical end time [s] of acceleration ramp. |
| `ROTATIONAL_SPEED_SCALING_FACTOR_INITIAL` | Float | *Required if active* | Initial rotational speed scale factor. |
| `ROTATIONAL_SPEED_SCALING_FACTOR_FINAL` | Float | *Required if active* | Final rotational speed scale factor. |

#### Localized Force Perturbations

Used to trigger rotating stall or unsteady instabilities:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `PERTURBATION_BODY_FORCE` | Boolean | `no` | Activates localized synthetic body force perturbation. |
| `PERTURBATION_CENTER` | 3 Floats | *Required if active* | Coordinates `x, y, z` [m] of the perturbation center. |
| `PERTURBATION_RADIAL_EXTENSION` | Float | *Required if active* | Spherical radius [m] of the perturbation bubble. |
| `PERTURBATION_SCALING_FACTOR` | Float | *Required if active* | Multiplier applied to the base body force in the zone. |
| `PERTURBATION_TIME_START` | Float | *Required if active* | Physical start time [s] of the perturbation. |
| `PERTURBATION_TIME_DURATION` | Float | *Required if active* | Duration [s] of the perturbation. |
| `GONG_MODELING_ACTIVE` | Boolean | `no` | Activates Gong BFM formulation terms. |

---

### 10. Greitzer Compression System Model

Simulates plenum and throttle dynamics coupled to the turbomachine to capture surge and stall limit cycles (requires `TIME_STEP_METHOD = fixed`):

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `ENABLE_GREITZER_MODELING` | Boolean | `no` | Activates downstream Greitzer plenum model. |
| `GREITZER_PLENUM_VOLUME` | Float | *Required if enabled* | Plenum volume $V_p$ [$\text{m}^3$] (alias: `GREITZER_PENUM_VOLUME`). |
| `GREITZER_AMBIENT_PRESSURE` | Float | `101325.0` | Ambient discharge pressure [$\text{Pa}$] outside throttle. |

---

### 11. Turbomachinery Diagnostics & Monitor Probes

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `TURBO_OUTPUT` | Boolean | `no` | Computes and logs turbomachinery performance (mass flow, total pressure ratio, total temperature ratio, isentropic efficiency). |
| `TURBO_DIRECTION` | String | `i` | Streamwise axis of the turbomachine: `i`, `-i`, `j`, `+j`, `-j`. |
| `MONITOR_POINTS_ACTIVE` | Boolean | `no` | Activates time-history probe recording ($p, u, v, w, t$). |
| `MONITOR_POINTS_I_COORDS` | Integer | *Required if active* | Grid $i$-index for the monitor point. |
| `MONITOR_POINTS_J_COORDS` | Integer | *Required if active* | Grid $j$-index for the monitor point. |
| `MONITOR_POINTS_NUMBER` | Integer | *Required if active* | Number of circumferentially distributed probe points around the annulus (3D). |

---

### 12. Output & Data Logging

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `SOLUTION_NAME` | String | `results` | Prefix name for output solution files. |
| `OUTPUT_FIELDS_TYPE` | String | `secondary` | Output variable set:<br>• `primary`: Conservative variables ($\rho, \rho u, \rho v, \rho w, \rho E$).<br>• `secondary`: Primitive and aerodynamic fields ($\rho, u, v, w, p, T, M, p_t, T_t, \mu_t$, etc.).<br>• `turbo_bfm`: Secondary fields plus BFM-specific variables (relative velocity, flow angles, deviation angles, inviscid and viscous forces, blockage). |
| `SAVE_UNSTEADY` | Boolean | `yes` | When `yes`, writes time-stamped intermediate solutions (`results_000250.csv`). |
| `SAVE_ITERATIONS_INTERVAL` | Integer | *Required* | Iteration interval for saving full volume solution files. |
| `HISTORY_BUFFER_SIZE` | Integer | `SAVE_ITERATIONS_INTERVAL` | Number of iterations buffered before flushing convergence logs and diagnostics to disk. |

---

## Example `input.ini` Configuration

Below is a complete annotated configuration for an axisymmetric turbomachinery compressor simulation using the Hall-Thollet BFM:

```ini
[CFD]
; --- Mesh and Geometry ---
GRID_FILE = mesh_compressor_axi.csv
TOPOLOGY = axisymmetric
SAVE_MESH_QUALITY_STATS = no

; --- Governing Equations and Schemes ---
KIND_SOLVER = Euler
CONVECTION_SCHEME = JST
MUSCL_RECONSTRUCTION = no
FLUX_LIMITER = None

; --- Fluid Properties (Ideal Air) ---
FLUID_MODEL = Ideal
FLUID_NAME = Air
FLUID_GAMMA = 1.4
FLUID_R_CONSTANT = 287.05
FLUID_KINEMATIC_VISCOSITY = 1.48e-5

; --- Boundary Conditions File & Patch Setup ---
BOUNDARY_CONDITIONS_FILEPATH = mesh_compressor_axi_boundaries.csv
INLET  = inlet, 101325, 288.15, 1, 0, 0
OUTLET = radial_equilibrium, 115000
HUB    = inviscid_wall
SHROUD = inviscid_wall

; --- Initialization ---
INIT_MACH_NUMBER = 0.45
INIT_PRESSURE = 101325
INIT_TEMPERATURE = 288.15
INIT_DIRECTION = 1, 0, 0
INLET_REFERENCE_FRAME = Cartesian

; --- Time Integration & Convergence ---
CFL = 1.1
TIME_INTEGRATION_TYPE = 0
TIME_STEP_METHOD = local
N_ITERATIONS = 50000
RESIDUALS_DROP_CONVERGENCE = 6
SIMULATION_IS_STEADY = yes

; --- Body Force Modeling (Hall-Thollet) ---
BFM_ACTIVE = yes
BLOCKAGE_ACTIVE = yes
BFM_MODEL = Hall-Thollet
HALL_THOLLET_COEFFICIENT_KN = 1.0
HALL_THOLLET_COEFFICIENT_KF = 1.0
HALL_THOLLET_COEFFICIENT_KD = 1.0
HALL_THOLLET_OFF_DESIGN_ACTIVE = no

; --- Turbomachinery Performance & Output ---
TURBO_OUTPUT = yes
TURBO_DIRECTION = i
OUTPUT_FIELDS_TYPE = turbo_bfm
SOLUTION_NAME = ger4_axi_solution
SAVE_UNSTEADY = no
SAVE_ITERATIONS_INTERVAL = 500
```

---

### Notes ###
* The code has been written and tested on macOS and Linux.
* Explicit Runge-Kutta time-stepping (3rd or 4th order) requires observing CFL stability limits.
* Restarting full-annulus 3D simulations from 2D axisymmetric meridional solutions is supported via `RESTART_TYPE = axisymmetric`.

---

### Contribution guidelines ###
* Validate modifications using unit and regression tests in `test/`.

---

### Authors and contacts ###
- **Francesco Neri**, TU Delft, `f.neri@tudelft.nl`
- **Matteo Pini**, TU Delft, `m.pini@tudelft.nl`

---

### References ###
[1] Hall, David Kenneth. *Analysis of civil aircraft propulsors with boundary layer ingestion*. Diss. Massachusetts Institute of Technology, 2015.  
[2] Thollet, William. *Body force modeling of fan-airframe interactions*. ISAE-SUPAERO: Toulouse, France (2017).  
[3] Chima, Rodrick V. *A three-dimensional unsteady CFD model of compressor stability*. Vol. 4241. 2006.  
[4] Tom-Robin Teschner, [CFD University](https://cfd.university/)