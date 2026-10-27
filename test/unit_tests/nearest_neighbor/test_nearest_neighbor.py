import subprocess
import shutil
import tempfile
from pathlib import Path
import pytest
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
SOLVER = BASE_DIR / "bin" / "turbobfm"


def test_couette_nearest_neighbor():
    """Test standard nearest neighbor on Couette flow case."""
    case_dir = BASE_DIR / "test" / "regression_tests" / "CouetteFlow_MovingWall"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        # Copy case files to temp directory
        for f in ["grid.csv", "boundaries.csv", "restart.csv", "input.ini"]:
            shutil.copy(case_dir / f, tmp_path / f)
        
        # Modify input.ini to use nearest_neighbor restart
        ini_content = (tmp_path / "input.ini").read_text()
        ini_content += "\nRESTART_TYPE = nearest_neighbor\nN_ITERATIONS = 50\n"
        (tmp_path / "input.ini").write_text(ini_content)
        
        # Run solver
        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0, f"Solver failed with output:\n{result.stdout}\n{result.stderr}"
        assert "Cartesian nearest-neighbor initialization done." in result.stdout
        print("Couette nearest_neighbor test PASSED.")


def test_coarse_to_fine_different_resolutions():
    """Test restarting fine grid simulation from a lower-resolution simulation."""
    case_dir = BASE_DIR / "test" / "regression_tests" / "CircularBump2D"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # 1. First obtain coarse results from CircularBump2D
        coarse_results = case_dir / "Volume_CSV" / "results.csv"
        assert coarse_results.exists(), "Coarse results.csv must exist"
        
        # 2. Setup a fine grid in tmp_path on the same domain [0, 1.5] x [0, 1]
        ni_fine, nj_fine, nk_fine = 35, 18, 1
        x = np.linspace(0.0, 1.5, ni_fine)
        y = np.linspace(0.0, 1.0, nj_fine)
        
        grid_lines = ["NDIMENSIONS=2\n", f"NI={ni_fine}\n", f"NJ={nj_fine}\n", f"NK={nk_fine}\n", "x,y,z\n"]
        for xi in x:
            for yj in y:
                grid_lines.append(f"{xi:.6f},{yj:.6f},0.000000\n")
        with open(tmp_path / "grid.csv", "w") as f:
            f.writelines(grid_lines)
        
        # Create matching boundaries.csv
        b_lines = [
            "NDIMENSIONS=2\n",
            f"NI={ni_fine}\n",
            f"NJ={nj_fine}\n",
            f"NK={nk_fine}\n",
            "NPATCHES=4\n",
            "PATCH_NAME,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX\n",
            f"INFLOW,0,0,0,{nj_fine},0,{nk_fine}\n",
            f"OUTFLOW,{ni_fine},{ni_fine},0,{nj_fine},0,{nk_fine}\n",
            f"UPPER_WALL,0,{ni_fine},{nj_fine},{nj_fine},0,{nk_fine}\n",
            f"LOWER_WALL,0,{ni_fine},0,0,0,{nk_fine}\n",
        ]
        with open(tmp_path / "boundaries.csv", "w") as f:
            f.writelines(b_lines)
        
        # Copy coarse results as restart file
        shutil.copy(coarse_results, tmp_path / "coarse_restart.csv")
        
        # Write input.ini
        input_ini = f"""[CFD]
GRID_FILE = grid.csv
SOLUTION_NAME = results

RESTART_SOLUTION = yes
RESTART_SOLUTION_FILEPATH = coarse_restart.csv
RESTART_TYPE = nearest_neighbor

KIND_SOLVER = Euler
TOPOLOGY = 2D

FLUID_GAMMA = 1.4
FLUID_NAME = Air
FLUID_MODEL = Ideal
FLUID_R_CONSTANT = 287.05 

BOUNDARY_CONDITIONS_FILEPATH = boundaries.csv
INFLOW = inlet, 140513.23, 316.224, 1, 0, 0
OUTFLOW = outlet, 101300
LOWER_WALL = inviscid_wall
UPPER_WALL = inviscid_wall

INIT_MACH_NUMBER = 0.5
INIT_TEMPERATURE = 285                   
INIT_PRESSURE = 101335                      
INIT_DIRECTION = 1, 0, 0     

CFL = 0.6
N_ITERATIONS = 20
RESIDUALS_DROP_CONVERGENCE = 4

SAVE_UNSTEADY = no
SAVE_ITERATIONS_INTERVAL = 250
TIME_INTEGRATION_TYPE = 0
TIME_STEP_METHOD = local
CONVECTION_SCHEME = Roe
MUSCL_RECONSTRUCTION = no
FLUX_LIMITER = van leer
"""
        (tmp_path / "input.ini").write_text(input_ini)
        
        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0, f"Solver failed with output:\n{result.stdout}\n{result.stderr}"
        assert "Cartesian nearest-neighbor initialization done." in result.stdout
        
        # Verify results were written and have fine grid dimensions
        res_file = tmp_path / "Volume_CSV" / "results.csv"
        assert res_file.exists(), "Results file should be generated"
        with open(res_file) as f:
            assert f.readline().strip() == f"NI={ni_fine}"
            assert f.readline().strip() == f"NJ={nj_fine}"
            assert f.readline().strip() == f"NK={nk_fine}"
            
        print("Coarse-to-fine different resolutions test PASSED.")


def test_coarse_to_fine_auto_nearest_neighbor():
    """Test restarting with different resolution when RESTART_TYPE is not explicitly specified (auto-fallback)."""
    case_dir = BASE_DIR / "test" / "regression_tests" / "CircularBump2D"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        coarse_results = case_dir / "Volume_CSV" / "results.csv"
        assert coarse_results.exists()
        
        ni_fine, nj_fine, nk_fine = 35, 18, 1
        x = np.linspace(0.0, 1.5, ni_fine)
        y = np.linspace(0.0, 1.0, nj_fine)
        
        grid_lines = ["NDIMENSIONS=2\n", f"NI={ni_fine}\n", f"NJ={nj_fine}\n", f"NK={nk_fine}\n", "x,y,z\n"]
        for xi in x:
            for yj in y:
                grid_lines.append(f"{xi:.6f},{yj:.6f},0.000000\n")
        with open(tmp_path / "grid.csv", "w") as f:
            f.writelines(grid_lines)
        
        b_lines = [
            "NDIMENSIONS=2\n",
            f"NI={ni_fine}\n",
            f"NJ={nj_fine}\n",
            f"NK={nk_fine}\n",
            "NPATCHES=4\n",
            "PATCH_NAME,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX\n",
            f"INFLOW,0,0,0,{nj_fine},0,{nk_fine}\n",
            f"OUTFLOW,{ni_fine},{ni_fine},0,{nj_fine},0,{nk_fine}\n",
            f"UPPER_WALL,0,{ni_fine},{nj_fine},{nj_fine},0,{nk_fine}\n",
            f"LOWER_WALL,0,{ni_fine},0,0,0,{nk_fine}\n",
        ]
        with open(tmp_path / "boundaries.csv", "w") as f:
            f.writelines(b_lines)
        
        shutil.copy(coarse_results, tmp_path / "coarse_restart.csv")
        
        # input.ini without RESTART_TYPE
        input_ini = f"""[CFD]
GRID_FILE = grid.csv
SOLUTION_NAME = results

RESTART_SOLUTION = yes
RESTART_SOLUTION_FILEPATH = coarse_restart.csv

KIND_SOLVER = Euler
TOPOLOGY = 2D

FLUID_GAMMA = 1.4
FLUID_NAME = Air
FLUID_MODEL = Ideal
FLUID_R_CONSTANT = 287.05 

BOUNDARY_CONDITIONS_FILEPATH = boundaries.csv
INFLOW = inlet, 140513.23, 316.224, 1, 0, 0
OUTFLOW = outlet, 101300
LOWER_WALL = inviscid_wall
UPPER_WALL = inviscid_wall

INIT_MACH_NUMBER = 0.5
INIT_TEMPERATURE = 285                   
INIT_PRESSURE = 101335                      
INIT_DIRECTION = 1, 0, 0     

CFL = 0.6
N_ITERATIONS = 20
RESIDUALS_DROP_CONVERGENCE = 4

SAVE_UNSTEADY = no
SAVE_ITERATIONS_INTERVAL = 250
TIME_INTEGRATION_TYPE = 0
TIME_STEP_METHOD = local
CONVECTION_SCHEME = Roe
MUSCL_RECONSTRUCTION = no
FLUX_LIMITER = van leer
"""
        (tmp_path / "input.ini").write_text(input_ini)
        
        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0, f"Solver failed with output:\n{result.stdout}\n{result.stderr}"
        assert "Using nearest-neighbor interpolation." in result.stdout
        assert "Cartesian nearest-neighbor initialization done." in result.stdout
        print("Coarse-to-fine auto nearest-neighbor test PASSED.")


def test_axisymmetric_2d_to_3d_restart():
    """Test restarting 3D simulation from 2D axisymmetric coarse simulation."""
    case_dir = BASE_DIR / "test" / "regression_tests" / "AxisymmetricSwirlingChannel"
    coarse_results = case_dir / "Volume_CSV" / "results.csv"
    assert coarse_results.exists()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        shutil.copy(coarse_results, tmp_path / "coarse_axisymmetric.csv")
        
        # Generate 3D annular grid with NI=20, NJ=10, NK=8
        ni, nj, nk = 20, 10, 8
        x = np.linspace(0.0, 1.0, ni)
        r = np.linspace(1.0, 2.0, nj)
        theta = np.linspace(0.0, np.pi/4, nk)  # 45 deg sector
        
        grid_lines = ["NDIMENSIONS=3\n", f"NI={ni}\n", f"NJ={nj}\n", f"NK={nk}\n", "x,y,z\n"]
        for xi in x:
            for rj in r:
                for tk in theta:
                    y = rj * np.cos(tk)
                    z = rj * np.sin(tk)
                    grid_lines.append(f"{xi:.6f},{y:.6f},{z:.6f}\n")
        with open(tmp_path / "grid.csv", "w") as f:
            f.writelines(grid_lines)
        
        # Boundary conditions
        b_lines = [
            "NDIMENSIONS=3\n",
            f"NI={ni}\n",
            f"NJ={nj}\n",
            f"NK={nk}\n",
            "NPATCHES=6\n",
            "PATCH_NAME,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX\n",
            f"INFLOW,0,0,0,{nj},0,{nk}\n",
            f"OUTFLOW,{ni},{ni},0,{nj},0,{nk}\n",
            f"LOWER_WALL,0,{ni},0,0,0,{nk}\n",
            f"UPPER_WALL,0,{ni},{nj},{nj},0,{nk}\n",
            f"PER1,0,{ni},0,{nj},0,0\n",
            f"PER2,0,{ni},0,{nj},{nk},{nk}\n",
        ]
        with open(tmp_path / "boundaries.csv", "w") as f:
            f.writelines(b_lines)
        
        input_ini = f"""[CFD]
GRID_FILE = grid.csv
SOLUTION_NAME = results

RESTART_SOLUTION = yes
RESTART_SOLUTION_FILEPATH = coarse_axisymmetric.csv
RESTART_TYPE = nearest_neighbor_axisymmetric

KIND_SOLVER = Euler
TOPOLOGY = 3D

FLUID_GAMMA = 1.4
FLUID_NAME = Air
FLUID_MODEL = Ideal
FLUID_R_CONSTANT = 287.05 

BOUNDARY_CONDITIONS_FILEPATH = boundaries.csv
INFLOW = inlet, 130000, 300, 1, 0, 0.5
OUTFLOW = radial_equilibrium, 101300
LOWER_WALL = inviscid_wall
UPPER_WALL = inviscid_wall

PER1 = periodic, 0, 0, 45
PER2 = periodic, 1, 0, 45

INIT_MACH_NUMBER = 0.5
INIT_TEMPERATURE = 300                   
INIT_PRESSURE = 101300                      
INIT_DIRECTION = 1, 0, 0.3                          

CFL = 1.0
N_ITERATIONS = 20
RESIDUALS_DROP_CONVERGENCE = 4

SAVE_UNSTEADY = no
SAVE_ITERATIONS_INTERVAL = 250
TIME_INTEGRATION_TYPE = 0
TIME_STEP_METHOD = local
CONVECTION_SCHEME = JST
"""
        (tmp_path / "input.ini").write_text(input_ini)
        
        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0, f"Solver failed with output:\n{result.stdout}\n{result.stderr}"
        assert "Axisymmetric nearest-neighbor initialization done." in result.stdout
        print("Axisymmetric 2D-to-3D restart test PASSED.")


def test_turbulent_flat_plate_restart_sa():
    """Test coarse-to-fine restart with Spalart-Allmaras turbulence model."""
    coarse_case = BASE_DIR / "testcases" / "Turbulent_FlatPlate" / "03_Run_roe1" / "Level_35x25"
    fine_case = BASE_DIR / "testcases" / "Turbulent_FlatPlate" / "03_Run_roe1" / "Level_69x49"
    coarse_results = coarse_case / "Volume_CSV" / "results.csv"

    if not coarse_results.exists():
        pytest.skip("Coarse results.csv for Turbulent Flat Plate not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        shutil.copy(coarse_results, tmp_path / "coarse_restart.csv")

        # Copy fine case input.ini and adjust
        ini_content = (fine_case / "input.ini").read_text()
        
        # Replace relative grid paths with absolute paths
        grid_file = (BASE_DIR / "testcases" / "Turbulent_FlatPlate" / "02_Convert_To_CTBFM" / "Output" / "grid_69x49_0.csv").resolve()
        bc_file = (BASE_DIR / "testcases" / "Turbulent_FlatPlate" / "02_Convert_To_CTBFM" / "Output" / "grid_boundaries_69x49_0.csv").resolve()

        lines = []
        for line in ini_content.splitlines():
            if line.startswith("GRID_FILE"):
                lines.append(f"GRID_FILE = {grid_file}")
            elif line.startswith("BOUNDARY_CONDITIONS_FILEPATH"):
                lines.append(f"BOUNDARY_CONDITIONS_FILEPATH = {bc_file}")
            elif line.startswith("RESTART_SOLUTION_FILEPATH"):
                lines.append("RESTART_SOLUTION_FILEPATH = coarse_restart.csv")
            elif line.startswith("RESTART_TYPE"):
                continue
            elif line.startswith("N_ITERATIONS"):
                lines.append("N_ITERATIONS = 10")
            else:
                lines.append(line)

        lines.append("RESTART_SOLUTION = yes")
        lines.append("RESTART_TYPE = nearest_neighbor")

        (tmp_path / "input.ini").write_text("\n".join(lines) + "\n")

        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Solver failed with output:\n{result.stdout}\n{result.stderr}"
        assert "Cartesian nearest-neighbor initialization done." in result.stdout
        assert "Turbulence model SA initialized from Cartesian nearest-neighbor restart file" in result.stdout
        print("Turbulent Flat Plate SA nearest-neighbor restart test PASSED.")


if __name__ == "__main__":
    test_couette_nearest_neighbor()
    test_coarse_to_fine_different_resolutions()
    test_coarse_to_fine_auto_nearest_neighbor()
    test_axisymmetric_2d_to_3d_restart()
    test_turbulent_flat_plate_restart_sa()
    print("\nALL NEAREST NEIGHBOR TESTS PASSED SUCCESSFULLY!")
