#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
from styles import *
import argparse
from Plot_MonitorPoints import read_monitor_points_data



def scale_oscillations(pressure, pressureRef, deltaTheta, reduction=1.5):
    """scale the oscillations to be nicely presented in the plot

    Args:
        pressure (array): values to scale
        pressureRef (array): reference values used to scale equivalently all the signals
        deltaTheta (float): offset to give in the value
        reduction (float, optional): factor that divides deltaTheta, to scale the peak-to-peak variations

    Returns:
        array: scaled signal
    """
    y = (pressure-pressureRef.min())/(pressureRef.max()-pressureRef.min())*deltaTheta/reduction
    y = y-y[0]
    return y


def plot_normalized_sensors_oscillation(sensors, time, time_per_revolution, name):
    """Plot the traces

    Args:
        sensors (list): list of arrays contaning the signal values
        time (array): time values
        time_per_revolution (float): used to find the equivalent revs
        name (string): name of the signal quantities
    """
    nSensors = len(sensors)
    plt.figure()
    for i in range(0, len(sensors)):
        deltaTheta = 360 / (nSensors)
        theta = i * 360 / (nSensors)
        plt.plot(time/time_per_revolution, scale_oscillations(sensors[i], sensors[-2], deltaTheta) + theta, label=f"probe {i}")

    if time_per_revolution != 1:
        plt.xlabel("Revs.")
    else:
        plt.xlabel(r"$t$ [s]")
    plt.ylabel(r"$\theta$ [deg]")
    plt.grid(alpha=.3)
    plt.title(name)
    plt.tight_layout()
    try:
        plt.savefig(f"Pictures/{name}_sensors.pdf", bbox_inches='tight')
    except:
        print("Could not save figure. Pictures folder not existing.")


def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Plot monitor point data with a given time per revolution. Use 1 to plot in seconds.")
    parser.add_argument("time_per_revolution", type=float, help="Time per revolution in seconds")
    args = parser.parse_args()
    time_per_revolution = args.time_per_revolution

    dfs = read_monitor_points_data()

    pressureSensors = []
    for df in dfs:
        pressureSensors.append(df["Pressure[Pa]"])
    time = dfs[0]["Time[s]"]
    
    plot_normalized_sensors_oscillation(pressureSensors, time, time_per_revolution, "Pressure")


if __name__ == "__main__":
    main()

plt.show()