#!/usr/bin/env python3
"""
generate_fluid_table.py
Generates structured thermodynamic look-up tables (.lut) using CoolProp for CTurboBFM.
Generates:
  1. (rho, e) grid -> [p, T, a, s, mu, kappa, dp_drho_e, dp_de_rho]
  2. (p, T) grid   -> [rho, e, a, s]
  3. (p, s) grid   -> [rho, T, e]
"""

import argparse
import sys
import numpy as np

try:
    import CoolProp.CoolProp as CP
except ImportError:
    print("Error: CoolProp is not installed. Please install it using 'pip install CoolProp'.")
    sys.exit(1)


def generate_table(fluid_name, p_min, p_max, T_min, T_max,
                   n_rho=80, n_e=80, n_p=50, n_T=50, n_s=50,
                   output_path=None, phase=None,
                   rho_min=None, rho_max=None,
                   e_min=None, e_max=None,
                   s_min=None, s_max=None):
    print(f"Generating fluid table for {fluid_name}...")
    print(f"Pressure range:    [{p_min:.2e}, {p_max:.2e}] Pa")
    print(f"Temperature range: [{T_min:.2f}, {T_max:.2f}] K")
    if phase:
        print(f"Phase filter:      {phase}")

    # Critical and reference properties
    try:
        T_crit = CP.PropsSI("Tcrit", fluid_name)
        p_crit = CP.PropsSI("pcrit", fluid_name)
        rho_crit = CP.PropsSI("rhocrit", fluid_name)
    except Exception as ex:
        print(f"Warning: Could not fetch critical properties: {ex}")
        T_crit, p_crit, rho_crit = 300.0, 1e5, 1.0

    molar_mass = CP.PropsSI("molar_mass", fluid_name)
    R_u = CP.PropsSI("gas_constant", fluid_name)
    R_gas = R_u / molar_mass

    # Sample corners and edge points to determine rho and e bounds
    test_p = np.linspace(p_min, p_max, 15)
    test_T = np.linspace(T_min, T_max, 15)
    rhos = []
    es = []
    ss = []
    for p in test_p:
        for T in test_T:
            try:
                if phase == "gas":
                    ph = CP.PhaseSI("P", p, "T", T, fluid_name)
                    if ph not in ["gas", "supercritical", "supercritical_gas"]:
                        continue
                rho_val = CP.PropsSI("D", "P", p, "T", T, fluid_name)
                e_val = CP.PropsSI("U", "P", p, "T", T, fluid_name)
                s_val = CP.PropsSI("S", "P", p, "T", T, fluid_name)
                if not (np.isnan(rho_val) or np.isnan(e_val) or np.isnan(s_val)):
                    rhos.append(rho_val)
                    es.append(e_val)
                    ss.append(s_val)
            except Exception:
                pass

    if len(rhos) == 0:
        raise ValueError(f"No valid state points found for fluid {fluid_name} in the specified (p, T) range.")

    calc_rho_min = min(rhos) * 0.70 if rho_min is None else rho_min
    calc_rho_max = max(rhos) * 1.50 if rho_max is None else rho_max
    calc_e_min = min(es) * 0.70 if e_min is None else e_min
    calc_e_max = max(es) * 1.15 if e_max is None else e_max
    calc_s_min = min(ss) * 0.95 if s_min is None else s_min
    calc_s_max = max(ss) * 1.05 if s_max is None else s_max

    rho_min, rho_max = calc_rho_min, calc_rho_max
    e_min, e_max = calc_e_min, calc_e_max
    s_min, s_max = calc_s_min, calc_s_max

    print(f"Computed density bounds: [{rho_min:.2f}, {rho_max:.2f}] kg/m3")
    print(f"Computed internal energy bounds: [{e_min:.2f}, {e_max:.2f}] J/kg")
    print(f"Computed entropy bounds: [{s_min:.2f}, {s_max:.2f}] J/(kg K)")

    # 1. Generate (rho, e) grid
    rho_grid = np.linspace(rho_min, rho_max, n_rho)
    e_grid = np.linspace(e_min, e_max, n_e)

    data_rho_e = []
    for rho in rho_grid:
        for e in e_grid:
            # Query P and T
            try:
                p = CP.PropsSI("P", "D", rho, "U", e, fluid_name)
                T = CP.PropsSI("T", "D", rho, "U", e, fluid_name)
            except Exception:
                T = max(T_min, (e - e_min) / (e_max - e_min) * (T_max - T_min) + T_min)
                p = rho * R_gas * T

            # Query S
            try:
                s = CP.PropsSI("S", "D", rho, "U", e, fluid_name)
            except Exception:
                try:
                    s = CP.PropsSI("S", "P", p, "T", T, fluid_name)
                except Exception:
                    try:
                        s = CP.PropsSI("S", "P", p, "Q", 1.0, fluid_name)
                    except Exception:
                        s = 1500.0

            # Query A (Speed of sound)
            try:
                a = CP.PropsSI("A", "D", rho, "U", e, fluid_name)
            except Exception:
                try:
                    a = CP.PropsSI("A", "P", p, "T", T, fluid_name)
                except Exception:
                    try:
                        a = CP.PropsSI("A", "P", p, "Q", 1.0, fluid_name)
                    except Exception:
                        # Approximation: a = sqrt(gamma * p / rho)
                        a = np.sqrt(max(1.0, 1.25 * p / rho))

            # Query transport properties
            try:
                mu = CP.PropsSI("V", "D", rho, "U", e, fluid_name)
            except Exception:
                try:
                    mu = CP.PropsSI("V", "P", p, "T", T, fluid_name)
                except Exception:
                    try:
                        mu = CP.PropsSI("V", "P", p, "Q", 1.0, fluid_name)
                    except Exception:
                        mu = 2.0e-5

            try:
                kappa = CP.PropsSI("L", "D", rho, "U", e, fluid_name)
            except Exception:
                try:
                    kappa = CP.PropsSI("L", "P", p, "T", T, fluid_name)
                except Exception:
                    try:
                        kappa = CP.PropsSI("L", "P", p, "Q", 1.0, fluid_name)
                    except Exception:
                        kappa = 0.03

            data_rho_e.append([p, T, a, s, mu, kappa])

    arr_rho_e = np.array(data_rho_e).reshape((n_rho, n_e, 6))

    # Compute numerical derivatives dp/drho|e and dp/de|rho
    dp_drho = np.zeros((n_rho, n_e))
    dp_de = np.zeros((n_rho, n_e))
    d_rho = rho_grid[1] - rho_grid[0]
    d_e = e_grid[1] - e_grid[0]

    for i in range(n_rho):
        ip = min(i + 1, n_rho - 1)
        im = max(i - 1, 0)
        denom_rho = (ip - im) * d_rho if ip != im else d_rho
        for j in range(n_e):
            jp = min(j + 1, n_e - 1)
            jm = max(j - 1, 0)
            denom_e = (jp - jm) * d_e if jp != jm else d_e

            dp_drho[i, j] = (arr_rho_e[ip, j, 0] - arr_rho_e[im, j, 0]) / denom_rho
            dp_de[i, j] = (arr_rho_e[i, jp, 0] - arr_rho_e[i, jm, 0]) / denom_e

    # 2. Generate (p, T) grid
    p_grid = np.linspace(p_min, p_max, n_p)
    T_grid = np.linspace(T_min, T_max, n_T)
    data_p_T = []
    for p in p_grid:
        for T in T_grid:
            try:
                if phase == "gas":
                    ph = CP.PhaseSI("P", p, "T", T, fluid_name)
                    if ph == "liquid":
                        rho = CP.PropsSI("D", "P", p, "Q", 1.0, fluid_name)
                        e = CP.PropsSI("U", "P", p, "Q", 1.0, fluid_name)
                        a = CP.PropsSI("A", "P", p, "Q", 1.0, fluid_name)
                        s = CP.PropsSI("S", "P", p, "Q", 1.0, fluid_name)
                    else:
                        rho = CP.PropsSI("D", "P", p, "T", T, fluid_name)
                        e = CP.PropsSI("U", "P", p, "T", T, fluid_name)
                        a = CP.PropsSI("A", "P", p, "T", T, fluid_name)
                        s = CP.PropsSI("S", "P", p, "T", T, fluid_name)
                else:
                    rho = CP.PropsSI("D", "P", p, "T", T, fluid_name)
                    e = CP.PropsSI("U", "P", p, "T", T, fluid_name)
                    a = CP.PropsSI("A", "P", p, "T", T, fluid_name)
                    s = CP.PropsSI("S", "P", p, "T", T, fluid_name)
            except Exception:
                rho = p / (R_gas * T)
                e = R_gas / 0.4 * T
                a = np.sqrt(1.4 * R_gas * T)
                s = 1000.0
            data_p_T.append([rho, e, a, s])

    # 3. Generate (p, s) grid
    s_grid = np.linspace(s_min, s_max, n_s)
    data_p_s = []
    for p in p_grid:
        for s in s_grid:
            try:
                rho = CP.PropsSI("D", "P", p, "S", s, fluid_name)
                T = CP.PropsSI("T", "P", p, "S", s, fluid_name)
                e = CP.PropsSI("U", "P", p, "S", s, fluid_name)
            except Exception:
                T = T_min
                rho = p / (R_gas * T)
                e = R_gas / 0.4 * T
            data_p_s.append([rho, T, e])

    if output_path is None:
        output_path = f"fluid_table_{fluid_name}.lut"

    with open(output_path, "w") as f:
        f.write("# CTurboBFM Fluid Look-Up Table\n")
        f.write(f"FLUID_NAME: {fluid_name}\n")
        f.write(f"MOLAR_MASS: {molar_mass:.6e}\n")
        f.write(f"GAS_CONSTANT: {R_gas:.6e}\n")
        f.write(f"CRITICAL_PRESSURE: {p_crit:.6e}\n")
        f.write(f"CRITICAL_TEMPERATURE: {T_crit:.6e}\n")
        f.write(f"CRITICAL_DENSITY: {rho_crit:.6e}\n")
        f.write("NOMINAL_GAMMA: 1.28\n\n")

        # Section 1: RHO_E
        f.write("# SECTION: RHO_E\n")
        f.write(f"N_RHO: {n_rho}\n")
        f.write(f"N_E: {n_e}\n")
        f.write(f"RHO_MIN: {rho_min:.6e}\n")
        f.write(f"RHO_MAX: {rho_max:.6e}\n")
        f.write(f"E_MIN: {e_min:.6e}\n")
        f.write(f"E_MAX: {e_max:.6e}\n")
        f.write("# p(Pa) T(K) a(m/s) s(J/kgK) mu(Pas) kappa(W/mK) dp_drho(Pa.m3/kg) dp_de(Pa.kg/J)\n")
        for i in range(n_rho):
            for j in range(n_e):
                p, T, a, s, mu, kappa = arr_rho_e[i, j]
                dpr = dp_drho[i, j]
                dpe = dp_de[i, j]
                f.write(f"{p:.6e} {T:.6e} {a:.6e} {s:.6e} {mu:.6e} {kappa:.6e} {dpr:.6e} {dpe:.6e}\n")
        f.write("\n")

        # Section 2: P_T
        f.write("# SECTION: P_T\n")
        f.write(f"N_P: {n_p}\n")
        f.write(f"N_T: {n_T}\n")
        f.write(f"P_MIN: {p_min:.6e}\n")
        f.write(f"P_MAX: {p_max:.6e}\n")
        f.write(f"T_MIN: {T_min:.6e}\n")
        f.write(f"T_MAX: {T_max:.6e}\n")
        f.write("# rho(kg/m3) e(J/kg) a(m/s) s(J/kgK)\n")
        for row in data_p_T:
            f.write(f"{row[0]:.6e} {row[1]:.6e} {row[2]:.6e} {row[3]:.6e}\n")
        f.write("\n")

        # Section 3: P_S
        f.write("# SECTION: P_S\n")
        f.write(f"N_P: {n_p}\n")
        f.write(f"N_S: {n_s}\n")
        f.write(f"P_MIN: {p_min:.6e}\n")
        f.write(f"P_MAX: {p_max:.6e}\n")
        f.write(f"S_MIN: {s_min:.6e}\n")
        f.write(f"S_MAX: {s_max:.6e}\n")
        f.write("# rho(kg/m3) T(K) e(J/kg)\n")
        for row in data_p_s:
            f.write(f"{row[0]:.6e} {row[1]:.6e} {row[2]:.6e}\n")

    print(f"Table successfully written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CoolProp thermodynamic look-up table for CTurboBFM.")
    parser.add_argument("--fluid", type=str, default="CO2", help="Fluid name (e.g. CO2, Nitrogen, MDM, R134a)")
    parser.add_argument("--phase", type=str, default=None, choices=["gas", "liquid", "all"], help="Phase filter (e.g. gas)")
    parser.add_argument("--p_min", type=float, default=5e6, help="Minimum pressure [Pa]")
    parser.add_argument("--p_max", type=float, default=12e6, help="Maximum pressure [Pa]")
    parser.add_argument("--T_min", type=float, default=305.0, help="Minimum temperature [K]")
    parser.add_argument("--T_max", type=float, default=450.0, help="Maximum temperature [K]")
    parser.add_argument("--rho_min", type=float, default=None, help="Explicit minimum density [kg/m3]")
    parser.add_argument("--rho_max", type=float, default=None, help="Explicit maximum density [kg/m3]")
    parser.add_argument("--e_min", type=float, default=None, help="Explicit minimum internal energy [J/kg]")
    parser.add_argument("--e_max", type=float, default=None, help="Explicit maximum internal energy [J/kg]")
    parser.add_argument("--s_min", type=float, default=None, help="Explicit minimum entropy [J/kg-K]")
    parser.add_argument("--s_max", type=float, default=None, help="Explicit maximum entropy [J/kg-K]")
    parser.add_argument("--n_rho", type=int, default=60, help="Number of density points")
    parser.add_argument("--n_e", type=int, default=60, help="Number of internal energy points")
    parser.add_argument("--n_p", type=int, default=40, help="Number of pressure points")
    parser.add_argument("--n_T", type=int, default=40, help="Number of temperature points")
    parser.add_argument("--n_s", type=int, default=40, help="Number of entropy points")
    parser.add_argument("--output", type=str, default=None, help="Output .lut file path")

    args = parser.parse_args()
    generate_table(
        fluid_name=args.fluid,
        p_min=args.p_min,
        p_max=args.p_max,
        T_min=args.T_min,
        T_max=args.T_max,
        n_rho=args.n_rho,
        n_e=args.n_e,
        n_p=args.n_p,
        n_T=args.n_T,
        n_s=args.n_s,
        output_path=args.output,
        phase=args.phase,
        rho_min=args.rho_min,
        rho_max=args.rho_max,
        e_min=args.e_min,
        e_max=args.e_max,
        s_min=args.s_min,
        s_max=args.s_max
    )
