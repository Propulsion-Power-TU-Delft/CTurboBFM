#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import argparse
from styles import *


def read_monitor_points_data():
    
    # --- File Setup ---
    os.makedirs("Pictures", exist_ok=True)
    cwd = os.getcwd()
    files = os.listdir(cwd + '/Monitor_Points')
    inputFiles = sorted([cwd + '/Monitor_Points/' + file for file in files if file.endswith('.csv')])

    # --- Read & Process Data ---
    dfs = []
    for inputFile in inputFiles:
        df = pd.read_csv(inputFile)
        df["Velocity_Magnitude[m/s]"] = np.sqrt(df["Velocity_X[m/s]"]**2 + df["Velocity_Y[m/s]"]**2 + df["Velocity_Z[m/s]"]**2)
        dfs.append(df)
    
    return dfs



def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Plot monitor point data with a given time per revolution. Use 1 to plot in seconds.")
    parser.add_argument("time_per_revolution", type=float, help="Time per revolution in seconds")
    args = parser.parse_args()
    time_per_revolution = args.time_per_revolution

    # read data
    dfs = read_monitor_points_data()

    # --- Fields to Plot ---
    fieldList = [
        'Pressure[Pa]',
        'Velocity_Magnitude[m/s]',
    ]
    
    fieldNames = [
        'Pressure',
        'Velocity',
    ]

    # --- Plotting ---
    for ii, field in enumerate(fieldList):
        plt.figure()
        
        if len(dfs) > 1:
            for i in range(1, len(dfs)):
                plt.plot(dfs[i]["Time[s]"]/time_per_revolution, dfs[i][field], label=f"Point {i}")
        else:
            plt.plot(dfs[0]["Time[s]"]/time_per_revolution, dfs[0][field])

        if time_per_revolution!=1:
            plt.xlabel("Revs.")
        else:
            plt.xlabel("Time [s]") # to give availability to plit in seconds if needed
        plt.ylabel(field)

        # plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"Pictures/Monitor_Points_{fieldNames[ii]}.pdf", bbox_inches='tight')
    
if __name__ == "__main__":
    main()
    
plt.show()

