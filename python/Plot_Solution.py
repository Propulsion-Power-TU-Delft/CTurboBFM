#!/usr/bin/env python3

#### Plot the turbo performance of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from styles import *

NLEVELS = 30

def GetSolutionDict(input_filename):
    # Read the first three lines to extract grid sizes
    with open(input_filename, 'r') as f:
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])
    
    df = pd.read_csv(input_filename, skiprows=3)
    data_dict = {col: df[col].to_numpy().reshape((ni, nj, nk)) for col in df.columns}
    return data_dict



def GetDataDict(input_filename, field):
    
    data_dict = GetSolutionDict(input_filename)
    
    if field not in data_dict.keys():
        # check if the field name with spaces exists in the dataset
        try:
            fieldName = field.replace('_', ' ')
            if fieldName not in data_dict.keys():
                raise ValueError
        except:
            raise ValueError(f"Field '{field}' not found in the data columns. Available fields: {list(data_dict.keys())}. For composed names, instead of spaces use underscores '_' in the field name.")
    else:
        fieldName = field
    new_dict = {}
    new_dict['x'] = data_dict['x']
    new_dict['y'] = data_dict['y']
    new_dict['z'] = data_dict['z']
    new_dict['axial'] = data_dict['x']
    new_dict['radial'] = np.sqrt(data_dict['y']**2 + data_dict['z']**2)
    new_dict['theta'] = np.mod(np.arctan2(data_dict['z'], data_dict['y']), 2*np.pi)
    new_dict[field] = data_dict[fieldName]
    return new_dict



def PlotData(data, iSlice, jSlice, kSlice, fieldName, gridType):
    
    # get the dimensionality
    nDim = 0
    
    if iSlice==':':
        nDim += 1
    if jSlice==':':
        nDim += 1
    if kSlice==':':
        nDim += 1
    
    if (nDim !=1 and nDim !=2):
        raise ValueError("Only 1D or 2D plots are supported. Please provide exactly one ':' slice.")
    
    if nDim==1:
        Plot1D(data, iSlice, jSlice, kSlice, fieldName)
    elif nDim==2 and gridType.lower() == "cylindrical":
        Plot2DCylindrical(data, iSlice, jSlice, kSlice, fieldName)
    elif nDim==2 and gridType.lower() == "computational":
        Plot2DComputational(data, iSlice, jSlice, kSlice, fieldName)




def Plot1D(data, iSlice, jSlice, kSlice, fieldName):
    if iSlice ==':':
        j = int(jSlice)
        k = int(kSlice)
        xlabel = r'i'
        values = data[fieldName][:,j,k]
        title = f"j-{j}_k-{k}"
    elif jSlice ==':':
        i = int(iSlice)
        k = int(kSlice)
        xlabel = r'j'
        values = data[fieldName][i,:,k]
        title = f"i-{i}_k-{k}"
    elif kSlice ==':':
        i = int(iSlice)
        j = int(jSlice)
        xlabel = r'k'
        values = data[fieldName][i,j,:]
        title = f"i-{i}_j-{j}"
    
    plt.figure()
    plt.plot(values, '-o', mfc='none')
    plt.xlabel(xlabel)
    plt.ylabel(fieldName)
    plt.grid(alpha=0.2)
    plt.title(f"{title}")
    plt.tight_layout()
    plt.savefig(f"Pictures/{fieldName}_1D_plot_{title}.pdf", bbox_inches='tight')




def Plot2DCylindrical(data, iSlice, jSlice, kSlice, fieldName):
    if iSlice == ':' and jSlice == ':':
        k = int(kSlice)
        x = data['axial'][:, :, k]
        y = data['radial'][:, :, k]
        values = data[fieldName][:, :, k]
        xlabel, ylabel = r'$x$ [m]', r'$r$ [m]'
        title = f"KPlane_{k}"

    elif iSlice == ':' and kSlice == ':':
        j = int(jSlice)
        x = data['axial'][:, j, :]
        y = data['radial'][:, j, :] * data['theta'][:, j, :]
        values = data[fieldName][:, j, :]
        xlabel, ylabel = r'$x$ [m]', r'$r\theta$ [m]'
        title = f"JPlane_{j}"

    elif jSlice == ':' and kSlice == ':':
        i = int(iSlice)
        x = data['y'][i, :, :]
        y = data['z'][i, :, :]
        values = data[fieldName][i, :, :]
        xlabel, ylabel = r'$y$ [m]', r'$z$ [m]'
        title = f"IPlane_{i}"

    else:
        raise ValueError("Exactly one slice index should be numeric, others ':'")
    
    plt.figure()
    plt.contourf(x, y, values, cmap = 'coolwarm', levels=NLEVELS)
    plt.colorbar()
    plt.title(f"{fieldName} at {title}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(f"Pictures/{fieldName}_2D_plot_{title}.pdf", bbox_inches='tight')



def Plot2DComputational(data, iSlice, jSlice, kSlice, fieldName):
    if iSlice == ':' and jSlice == ':':
        k = int(kSlice)
        values = data[fieldName][:, :, k]
        xlabel, ylabel = r'I', r'J'
        title = f"KPlane_{k}"

    elif iSlice == ':' and kSlice == ':':
        j = int(jSlice)
        values = data[fieldName][:, j, :]
        xlabel, ylabel = r'I', r'K'
        title = f"JPlane_{j}"

    elif jSlice == ':' and kSlice == ':':
        i = int(iSlice)
        values = data[fieldName][i, :, :]
        xlabel, ylabel = r'J', r'K'
        title = f"IPlane_{i}"

    else:
        raise ValueError("Exactly one slice index should be numeric, others ':'")
    
    plt.figure()
    plt.imshow(values.T, cmap='turbo', interpolation='bicubic', origin='lower')
    plt.colorbar()
    plt.title(f"{fieldName} at {title}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(f"Pictures/{fieldName}_2D_plot_{title}.pdf", bbox_inches='tight')
        
    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the solution field of a CTurboBFM result.")
    
    parser.add_argument("inputFile", type=str, help="Path to the CSV file containing the simulation results.")
    parser.add_argument("iSlice", type=str, help="Slice along i-index")
    parser.add_argument("jSlice", type=str, help="Slice along j-index")
    parser.add_argument("kSlice", type=str, help="Slice along k-index")
    parser.add_argument("fieldName", type=str, help="Field name to plot")
    parser.add_argument("gridType", type=str, nargs='?', default="cylindrical",
                        help="Grid type for 2D plots, either 'cylindrical' or 'computational' (default: cylindrical)")
    
    args = parser.parse_args()

    inputFile = args.inputFile
    iSlice = args.iSlice
    jSlice = args.jSlice
    kSlice = args.kSlice
    fieldName = args.fieldName
    gridType = args.gridType
    
    os.makedirs("Pictures", exist_ok=True)
    data = GetDataDict(inputFile, fieldName)
    PlotData(data, iSlice, jSlice, kSlice, fieldName, gridType)
    plt.show()
    