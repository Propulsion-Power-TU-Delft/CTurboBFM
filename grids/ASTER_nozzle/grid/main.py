import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from unsflow.grid.functions import plot_twodimensional_grid

xThroat = 3.005                         
xStart = xThroat - (7.188)*1E-3         
xEnd = xThroat + 23.412*1E-3            
radiusThroat = 0.022/2
radiusPipe = 0.0303/2
radiusEnd = radiusPipe
rescaleAreaRatios = 0.527
ni = 100
nj = 50


x = np.linspace(xStart, xEnd, ni)
y = np.interp(x, [xStart, xThroat, xEnd], [radiusPipe, radiusThroat, radiusEnd])
X, Y = np.zeros((ni, nj)), np.zeros((ni, nj))
for i in range(ni):
    X[i, :] = x[i]
    Y[i, :] = np.linspace(0, y[i], nj)


# Create a 3D scatter plots
plot_twodimensional_grid(X, Y, frame='cartesian')


# Create output directory
OUTPUT_FOLDER = 'Grid'
os.makedirs(OUTPUT_FOLDER, exist_ok=True)
with open(OUTPUT_FOLDER + '/grid_%02ix%02i.csv' %(ni, nj), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={ni}\n")
    file.write(f"NJ={nj}\n")
    file.write(f"NK=1\n")
    file.write("x,y,z\n")
    for i in range(ni):
        for j in range(nj):
            for k in range(1):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f}\n")

plt.show()







