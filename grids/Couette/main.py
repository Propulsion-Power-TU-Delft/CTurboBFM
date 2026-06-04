import numpy as np
import matplotlib.pyplot as plt
import os
import pickle




"""
Generate a 2D rectangular geometry, that will be used as verification
"""
OUTPUT_FOLDER = 'Grid'

LZ = 1     # flow direction
LY = LZ/5 
LX = LZ/5
NX = 6          # not needed direction
NY = 65
NZ = 6 

x = np.linspace(0, LX, NX)
y = np.linspace(0, LY, NY)
z = np.linspace(0, LZ, NZ)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
NI, NJ, NK = X.shape

# Create output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_%02i.csv' %(NI, NJ, NK), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK={NK}\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                file.write(f"{X[i,j,k]:.6f},{Y[i,j,k]:.6f},{Z[i,j,k]:.6f}\n")






