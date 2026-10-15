import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pickle
from unsflow.grid.functions import transfinite_grid_generation, plot_twodimensional_grid


OUTPUT_FOLDER = 'Grid'
NI = 75
NJ = 50
NK = 1
STREAMWISE_COEFF = 1
SPANWISE_COEFF = 1
ANGLE_DEG = 10

angle_rad = 10*np.pi/180

x1 = np.linspace(0, 0.5, NI//3)
y1 = np.zeros(NI//3)

x2 = np.linspace(0.5, 1.5, NI-NI//3)
y2 = (x2-0.5)*np.tan(angle_rad)

x_low = np.concatenate((x1, x2[1:]), axis=0)
y_low = np.concatenate((y1, y2[1:]), axis=0)

NI = len(x_low)

x_inlet = np.zeros(NJ)
y_inlet = np.linspace(0, 1, NJ)

x_outlet = np.zeros(NJ)
y_outlet = np.linspace(y_low[-1], 1, NJ)

x_top = x_low.copy()
y_top = np.zeros(NI)+1


X, Y = transfinite_grid_generation(
    np.vstack((x_inlet, y_inlet)), 
    np.vstack((x_low, y_low)), 
    np.vstack((x_outlet, y_outlet)), 
    np.vstack((x_top, y_top)),
    block_topology='internal',
    streamwise_coeff=1.0, 
    spanwise_coeff=1.0,
    )


# Create a 3D scatter plots
mesh = pv.StructuredGrid(X, Y, np.zeros_like(X))
plotter = pv.Plotter()
plotter.add_mesh(mesh, cmap='viridis', show_edges=True)
plotter.show_axes()


# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_%02i.csv' %(NI, NJ, NK), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK={NK}\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f}\n")
    print(f"Grid saved to {OUTPUT_FOLDER}/grid_%02i_%02i_%02i.csv" %(NI, NJ, NK))


plotter.show()
plt.show()







