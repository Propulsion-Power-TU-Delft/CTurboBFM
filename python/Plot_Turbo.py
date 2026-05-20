#!/usr/bin/env python3

#### Plot the turbo performance of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from styles import *

def read_file(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def shift_to_zero(res):
    return res - res[0]

def plot_turbo(data):
    
    # beta-m
    fig, ax1 = plt.subplots()
    ax1.plot(data["Massflow[kg/s]"], color='C0')
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel(r"$\dot{m}$ [kg/s]", color='C0')
    ax1.tick_params(axis='y', labelcolor='C0')
    ax1.grid(alpha=0.2)
    ax2 = ax1.twinx()
    ax2.plot(data["PRtt"], color='C1')
    ax2.set_ylabel(r"$\beta_{\rm tt}$", color='C1')
    ax2.tick_params(axis='y', labelcolor='C1')
    plt.tight_layout()
    plt.savefig("Pictures/massflow_betatt.pdf", bbox_inches='tight')
    
    # eta-m plot
    fig, ax1 = plt.subplots()
    ax1.plot(data["Massflow[kg/s]"], color='C0')
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel(r"$\dot{m}$ [kg/s]", color='C0')
    ax1.tick_params(axis='y', labelcolor='C0')
    ax1.grid(alpha=0.2)
    ax2 = ax1.twinx()
    ax2.plot(data["ETAtt"], color='C2')
    ax2.set_ylabel(r"$\eta_{\rm tt}$", color='C2')
    ax2.tick_params(axis='y', labelcolor='C2')
    plt.tight_layout()
    plt.savefig("Pictures/massflow_etatt.pdf", bbox_inches='tight')
    
    # beta-eta plot
    fig, ax1 = plt.subplots()
    ax1.plot(data["PRtt"], color='C1')
    ax1.set_xlabel("Iteration [-]")
    ax1.set_ylabel(r"$\beta_{\rm tt}$ [-]", color='C1')
    ax1.tick_params(axis='y', labelcolor='C1')
    ax1.grid(alpha=0.2)
    ax2 = ax1.twinx()
    ax2.plot(data["ETAtt"], color='C2')
    ax2.set_ylabel(r"$\eta_{\rm tt}$ [-]", color='C2')
    ax2.tick_params(axis='y', labelcolor='C2')
    plt.tight_layout()
    plt.savefig("Pictures/betatt_etatt.pdf", bbox_inches='tight')
    
    # beta evolution plot
    plt.figure()
    plt.plot(data["Massflow[kg/s]"], data["PRtt"], '-k', lw=0.5)
    plt.scatter(data["Massflow[kg/s]"], data["PRtt"], c=np.linspace(1, len(data["Massflow[kg/s]"]), len(data["Massflow[kg/s]"])), s=5, cmap='turbo')
    plt.grid(alpha=0.2)
    plt.xlabel(r"$\dot{m}$ [kg/s]")
    plt.ylabel(r"$\beta_{\rm tt}$")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("Pictures/betaTrajectory.pdf", bbox_inches='tight')
    
    # eta evolution plot
    plt.figure()
    plt.plot(data["Massflow[kg/s]"], data["ETAtt"], '-k', lw=0.5)
    plt.scatter(data["Massflow[kg/s]"], data["ETAtt"], c=np.linspace(1, len(data["Massflow[kg/s]"]), len(data["Massflow[kg/s]"])), s=5, cmap='turbo')
    plt.grid(alpha=0.2)
    plt.xlabel(r"$\dot{m}$ [kg/s]")
    plt.ylabel(r"$\eta_{\rm tt}$")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("Pictures/etaTrajectory.pdf", bbox_inches='tight')
        
    plt.show()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the performance of a CTurboBFMsimulation.")
    parser.add_argument("input_file", help="Path to the turbo performance CSV file")
    args = parser.parse_args()

    input_filename = args.input_file
    
    data = read_file(input_filename)
    os.makedirs("Pictures", exist_ok=True)
    plot_turbo(data)
