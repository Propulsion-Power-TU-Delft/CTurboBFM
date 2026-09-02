#!/usr/bin/env python3
"""
Standalone regression test runner for CTurboBFM.
Runs all test cases in test/regression_tests and validates against reference.csv.
Requires only standard Python 3 (no third-party dependencies).
"""

import os
import sys
import subprocess
import math
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
CASES_DIR = BASE_DIR / "test" / "regression_tests"
SOLVER = BASE_DIR / "bin" / "turbobfm"


def read_structured_csv(filename):
    """Read a structured CSV file with dimension headers."""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    start_idx = 0
    while start_idx < len(lines) and '=' in lines[start_idx]:
        start_idx += 1

    if start_idx >= len(lines):
        raise ValueError(f"No header found in {filename}")

    header = [h.strip() for h in lines[start_idx].split(',')]
    data = {h: [] for h in header}

    for line in lines[start_idx + 1:]:
        parts = [p.strip() for p in line.split(',')]
        for h, v in zip(header, parts):
            data[h].append(float(v))

    return data


def compute_relative_l2_error(ref_vals, res_vals):
    """Compute relative L2 norm difference: ||ref - res|| / ||ref||"""
    if len(ref_vals) != len(res_vals):
        raise ValueError(f"Length mismatch: {len(ref_vals)} vs {len(res_vals)}")

    diff_sq_sum = sum((a - b) ** 2 for a, b in zip(ref_vals, res_vals))
    ref_sq_sum = sum(a ** 2 for a in ref_vals)

    if ref_sq_sum == 0:
        return math.sqrt(diff_sq_sum)

    return math.sqrt(diff_sq_sum) / math.sqrt(ref_sq_sum)


def run_test_case(case_dir, max_tol=0.02):
    """Run solver on a test case and compare against reference."""
    case_name = case_dir.name
    log_file = case_dir / "log_test.txt"

    print(f"\n[{case_name}] Running solver...")
    with open(log_file, "w") as log_out:
        result = subprocess.run(
            [str(SOLVER), "input.ini"],
            cwd=case_dir,
            stdout=log_out,
            stderr=subprocess.STDOUT
        )

    if result.returncode != 0:
        print(f"  ❌ FAILED: Solver returned non-zero exit code ({result.returncode}). See {log_file}")
        return False

    ref_file = case_dir / "reference.csv"
    res_file = case_dir / "Volume_CSV" / "results.csv"

    if not ref_file.exists():
        print(f"  ⚠️ SKIPPED: Reference file {ref_file} not found.")
        return True

    if not res_file.exists():
        print(f"  ❌ FAILED: Result file {res_file} was not generated.")
        return False

    ref_data = read_structured_csv(ref_file)
    res_data = read_structured_csv(res_file)

    keys_to_check = ['Density', 'Pressure', 'Temperature']
    all_passed = True

    for key in keys_to_check:
        if key not in ref_data or key not in res_data:
            print(f"  ⚠️ Key '{key}' missing from data files.")
            continue

        err = compute_relative_l2_error(ref_data[key], res_data[key])
        if err <= max_tol:
            print(f"  ✅ {key:<12}: L2 Rel Error = {err:.4e} (tol <= {max_tol:.2e})")
        else:
            print(f"  ❌ {key:<12}: L2 Rel Error = {err:.4e} (EXCEEDS tol {max_tol:.2e})")
            all_passed = False

    return all_passed


def main():
    if not SOLVER.exists():
        print(f"Error: Solver executable not found at {SOLVER}. Run 'make release' first.")
        sys.exit(1)

    cases = sorted([d for d in CASES_DIR.iterdir() if d.is_dir() and (d / "input.ini").exists()])
    print(f"Discovered {len(cases)} regression test cases in {CASES_DIR}:")
    for c in cases:
        print(f"  - {c.name}")

    results = {}
    for case_dir in cases:
        passed = run_test_case(case_dir)
        results[case_dir.name] = passed

    print("\n" + "=" * 50)
    print("REGRESSION TEST SUMMARY:")
    print("=" * 50)
    all_passed = True
    for name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        icon = "✅" if passed else "❌"
        print(f" {icon} {name:<42} : {status}")
        if not passed:
            all_passed = False

    print("=" * 50)
    if all_passed:
        print("🎉 ALL REGRESSION TESTS PASSED!")
        sys.exit(0)
    else:
        print("❌ SOME REGRESSION TESTS FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
