import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from Utils.styles import *


files = [
    '../Run_Roe1/Volume_CSV/results_250_roe1.csv',
    '../Run_Roe2/Volume_CSV/results_250_roe2.csv',
]

labels = [
    'Roe 1',
    'Roe 2',
]

fields = [
    'Density',
    'Velocity X',
    'Pressure',
]

for field in fields:
    plt.figure()
    for file in files:
        df = pd.read_csv(file, skiprows=3)
        plt.plot(df['x'], df[field], label=labels[files.index(file)])
    plt.xlabel(r'$x$')
    plt.ylabel(r'%s' % field)
    plt.grid(alpha=0.2)
    plt.legend()
    plt.savefig('%s.pdf' % field, bbox_inches='tight')
plt.show()