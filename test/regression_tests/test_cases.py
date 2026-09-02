import subprocess
from pathlib import Path
import pytest
import pandas as pd
import numpy as np


BASE_DIR = Path(__file__).resolve().parent.parent.parent
CASES_DIR = BASE_DIR / "test" / "regression_tests"
SOLVER = BASE_DIR / "bin" / "turbobfm"


def discover_cases():
    cases = []
    for d in CASES_DIR.iterdir():
        if d.is_dir() and (d / "input.ini").exists():
            cases.append(d)
    
    print(f"Discovered {len(cases)} test cases:")
    for case in cases:
        print(f" - {case.name}")
    return cases


def run_solver(case_dir):
    log_file = case_dir / "log_test.txt"

    result = subprocess.run(
        [str(SOLVER), "input.ini"],
        cwd=case_dir,
        stdout=open(log_file, "w"),
        stderr=subprocess.STDOUT
    )
    
    assert result.returncode == 0, f"Case {case_dir.name} did not complete successfully. Check {log_file} for details."

def read_structured_csv(filename):
    with open(filename, 'r') as f:
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])

    # Skip only the first 3 lines (dimensions), keep header line
    df = pd.read_csv(filename, skiprows=3)

    # Check consistency
    n_points = ni * nj * nk
    assert len(df) == n_points, f"Expected {n_points} points but found {len(df)}"

    shape = (ni, nj, nk)

    data = {
        col: df[col].to_numpy().reshape(shape, order='C')
        for col in df.columns
    }

    return data

def result_arrays_agree(reference_data, result_data, tol=1e-3):
    err = np.linalg.norm(reference_data - result_data) / np.linalg.norm(reference_data)
    if err<tol:
        return True, err
    else:
        return False, err


def compare_results_with_reference(case_dir):
    reference_file = case_dir / "reference.csv"
    result_file = case_dir / "Volume_CSV" / "results.csv"

    # check they exist
    assert reference_file.exists(), f"Reference file {reference_file} does not exist."
    assert result_file.exists(), f"Result file {result_file} does not exist."

    reference_data = read_structured_csv(reference_file)
    result_data = read_structured_csv(result_file)
    
    for idim in range(3):
        assert reference_data['x'].shape[idim] == result_data['x'].shape[idim], \
            f"Dimension {idim} size mismatch: reference {reference_data['x'].shape[idim]} vs result {result_data['x'].shape[idim]}"

    keys_to_check = ['Density', 'Pressure', 'Temperature']
    max_tols = [2e-2, 2e-2, 2e-2]
    for ikey, key in enumerate(keys_to_check):
        ref_value = reference_data[key]
        res_value = result_data[key]
        they_agree, err = result_arrays_agree(ref_value, res_value, tol=max_tols[ikey])
        assert they_agree, f"Mismatch in {key} with relative error {err:.9f}"
        print(f"{key} matches with relative error {err:.9f}")


@pytest.mark.parametrize("case_dir", discover_cases())
def test_run_case(case_dir):
    print(f"\n=== Running test in {case_dir.name} ===")

    log_file = case_dir / "log_test.txt"

    run_solver(case_dir)
    
    compare_results_with_reference(case_dir)