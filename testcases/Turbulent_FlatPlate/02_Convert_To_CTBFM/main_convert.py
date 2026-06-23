import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from utils import read_grid_file, plot_grid, save_pickle, save_csv_grid, save_csv_boundaries

os.makedirs('Output', exist_ok=True)
input_file = "../01_Nasa_Grids/flatplate_clust2_4levelsdown_35x25.p2dfmt"

X, Y = read_grid_file(input_file)
for iblock in range(len(X)):
    plot_grid(X, Y, iblock)
    # save_pickle(X, Y, iblock)
    save_csv_grid(X, Y, iblock)
    save_csv_boundaries(X, Y, iblock)
    

plt.show()