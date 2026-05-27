#!/usr/bin/env python3

#### Convert structured CSV data to VTS format. Given the path to the csv results file (obtained with CTurboBFM), outputs a VTS file with the same name, in the same location.

import numpy as np
import pandas as pd
from pyevtk.hl import gridToVTK
import argparse
import os
import sys

def read_structured_csv(filename):
    with open(filename, 'r') as f:
        # Read grid dimensions
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
        col: df[col].to_numpy().reshape(shape, order='C')  # i-fastest
        for col in df.columns
    }

    return data


def writeVTK(data, filename):
        
    pointsData = {
        "Velocity": (
            np.ascontiguousarray(data['Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Velocity Z'], dtype=np.float64)
        )
    }

    # BFM output vector fields
    try:
        pointsData["Grid Velocity"] = (
            np.ascontiguousarray(data['Grid Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Grid Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Grid Velocity Z'], dtype=np.float64)
        )


        pointsData["Relative Velocity"] = (
            np.ascontiguousarray(data['Relative Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Relative Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Relative Velocity Z'], dtype=np.float64)
        )
        
        pointsData["Inviscid Body Force"] = (
            np.ascontiguousarray(data['Inviscid Body Force X'], dtype=np.float64),
            np.ascontiguousarray(data['Inviscid Body Force Y'], dtype=np.float64),
            np.ascontiguousarray(data['Inviscid Body Force Z'], dtype=np.float64)
        )
        
        pointsData["Viscous Body Force"] = (
            np.ascontiguousarray(data['Viscous Body Force X'], dtype=np.float64),
            np.ascontiguousarray(data['Viscous Body Force Y'], dtype=np.float64),
            np.ascontiguousarray(data['Viscous Body Force Z'], dtype=np.float64)
        )
    except:
        pass
    
    # Turbulence/Viscous output vector fields
    try:
        pointsData["Wall Shear Stress"] = (
            np.ascontiguousarray(data['Wall Shear Stress X'], dtype=np.float64),
            np.ascontiguousarray(data['Wall Shear Stress Y'], dtype=np.float64),
            np.ascontiguousarray(data['Wall Shear Stress Z'], dtype=np.float64)
        )
    except:
        pass
        
    
    for key in data.keys():
        if key not in ['Velocity X', 'Velocity Y', 'Velocity Z', 'x', 'y', 'z',
                       'Grid Velocity X', 'Grid Velocity Y', 'Grid Velocity Z',
                       'Relative Velocity X', 'Relative Velocity Y', 'Relative Velocity Z',
                       'Inviscid Body Force X', 'Inviscid Body Force Y', 'Inviscid Body Force Z',
                       'Viscous Body Force X', 'Viscous Body Force Y', 'Viscous Body Force Z',
                       'Wall Shear Stress X', 'Wall Shear Stress Y', 'Wall Shear Stress Z']:
            pointsData[key] = np.ascontiguousarray(data[key], dtype=np.float64)

    gridToVTK(
        filename,
        np.ascontiguousarray(data['x']),
        np.ascontiguousarray(data['y']),
        np.ascontiguousarray(data['z']),
        pointData=pointsData
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert structured CSV data to VTK format.")
    parser.add_argument("input_path", help="Path to the input CSV file or folder")
    args = parser.parse_args()

    input_path = args.input_path
    vts_output_dir = "Volume_VTS"
    os.makedirs(vts_output_dir, exist_ok=True)

    if os.path.isfile(input_path) and input_path.endswith(".csv"):
        # Handle single CSV file
        data = read_structured_csv(input_path)
        base_name = os.path.splitext(os.path.basename(input_path))[0]
        output_filename = os.path.join(vts_output_dir, base_name)
        writeVTK(data, output_filename)
        print(f"Written VTS file to {output_filename}.vts")

    elif os.path.isdir(input_path):
        # Handle directory: convert all .csv files inside
        csv_files = [f for f in os.listdir(input_path) if f.endswith(".csv")]
        if not csv_files:
            print("No CSV files found in directory:", input_path)
            sys.exit(1)

        for csv_file in csv_files:
            full_csv_path = os.path.join(input_path, csv_file)
            data = read_structured_csv(full_csv_path)
            base_name = os.path.splitext(csv_file)[0]
            output_filename = os.path.join(vts_output_dir, base_name)
            writeVTK(data, output_filename)
            print(f"Written VTS file to {output_filename}.vts")
    else:
        print("Invalid input. Please provide a CSV file or a directory containing CSV files.")
        sys.exit(1)