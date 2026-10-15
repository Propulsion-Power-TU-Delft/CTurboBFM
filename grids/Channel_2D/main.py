import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pickle
from unsflow.grid.functions import transfinite_grid_generation, plot_twodimensional_grid




"""
Generate a 2D rectangular geometry, that will be used as verification
"""
OUTPUT_FOLDER = 'Grid'
LX = 3
R_IN = 1
R_OUT = 2
NX = 128
NY = 128

x = np.linspace(0, LX, NX)
y = np.linspace(R_IN, R_OUT, NY)



X, Y = np.meshgrid(x, y, indexing='ij')
NI, NJ = X.shape


# Create a 3D scatter plots
mesh = pv.StructuredGrid(X, Y, np.zeros_like(X))
plotter = pv.Plotter()
plotter.add_mesh(mesh, cmap='viridis', show_edges=True)
plotter.show_axes()


grid = {'X': X, 'Y': Y}
# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i.csv' %(NX, NY), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK=1\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(1):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f}\n")


plotter.show()
plt.show()







