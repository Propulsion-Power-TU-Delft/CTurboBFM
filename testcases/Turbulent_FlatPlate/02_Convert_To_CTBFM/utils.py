import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

def read_grid_file(filename):
    with open(filename, 'r') as f:
        tokens = f.read().split()
    
    idx = 0
    def next_val():
        nonlocal idx
        val = tokens[idx]
        idx += 1
        return val

    # Read number of blocks
    nbl = int(next_val())
    print(f"Number of blocks: {nbl}")

    # Read dimensions for each block (2D: only idim and jdim)
    idim = np.zeros(nbl, dtype=int)
    jdim = np.zeros(nbl, dtype=int)
    for n in range(nbl):
        idim[n] = int(next_val())
        jdim[n] = int(next_val())
        print(f"Block {n}: {idim[n]} x {jdim[n]}")

    # Read x, y for each block
    X = []
    Y = []
    for n in range(nbl):
        ni, nj = idim[n], jdim[n]

        x = np.zeros((ni, nj))
        y = np.zeros((ni, nj))

        # i is innermost, j is outermost
        for j in range(nj):
            for i in range(ni):
                x[i, j] = float(next_val())
        for j in range(nj):
            for i in range(ni):
                y[i, j] = float(next_val())

        X.append(x)
        Y.append(y)

    return X, Y

def plot_grid(X, Y, iblock):
    plt.figure(figsize=(8, 6))
    ni,nj = X[iblock].shape
    for i in range(ni):
        plt.plot(X[iblock][i, :], Y[iblock][i, :], 'k-', lw=0.5)  
    for j in range(nj):
        plt.plot(X[iblock][:, j], Y[iblock][:, j], 'k-', lw=0.5)  
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.title(f'Grid Block {iblock}')
    plt.savefig(f'Output/grid_plot_block_{iblock}.png', dpi=300)

def save_pickle(X, Y, iblock):
    ni, nj = X[iblock].shape
    with open(f"Output/grid_block_{ni}x{nj}_{iblock}.pkl", 'wb') as f:
        pickle.dump((X[iblock], Y[iblock]), f)
    print(f"Saved grid block {iblock} to Output/grid_block_{ni}x{nj}_{iblock}.pkl")

def save_csv_grid(X, Y, iblock):
    ni,nj = X[iblock].shape
    filename = f'Output/grid_{ni}x{nj}_{iblock}.csv'
    with open(filename, 'w') as file:
        file.write(f"NDIMENSIONS=2\n")
        file.write(f"NI={ni}\n")
        file.write(f"NJ={nj}\n")
        file.write(f"NK=1\n")
        file.write("x,y,z\n")
        for i in range(ni):
            for j in range(nj):
                for k in range(1):
                    file.write(f"{X[iblock][i,j]:.17g},{Y[iblock][i,j]:.17g},{0:.17g}\n")
    print(f"Saved grid block {iblock} to {filename}")
    
def save_csv_boundaries(X, Y, iblock):
    """Specialized bcs for the NASA 2D bump
    """
    ni,nj = X[iblock].shape
    filename = f'Output/grid_boundaries_{ni}x{nj}_{iblock}.csv'
    
    i_start_plate = np.argmin(np.abs(X[iblock][:,0] - 0.0)) # this is the primary node of the bump start
    
    # define the patches and their corresponding i,j,k ranges of the dual grid nodes (n+1 nodes)
    patches = {}
    patches["INFLOW"] = {"I_MIN": 0, "I_MAX": 0, "J_MIN": 0, "J_MAX": nj, "K_MIN": 0, "K_MAX": 1}
    patches["OUTFLOW"] = {"I_MIN": ni, "I_MAX": ni, "J_MIN": 0, "J_MAX": nj, "K_MIN": 0, "K_MAX": 1}
    patches["UPPER_WALL"] = {"I_MIN": 0, "I_MAX": ni, "J_MIN": nj, "J_MAX": nj, "K_MIN": 0, "K_MAX": 1}
    patches["SYMMETRIC"] = {"I_MIN": 0, "I_MAX": i_start_plate+1, "J_MIN": 0, "J_MAX": 0, "K_MIN": 0, "K_MAX": 1}
    patches["FLAT_PLATE"] = {"I_MIN": i_start_plate+1, "I_MAX": ni, "J_MIN": 0, "J_MAX": 0, "K_MIN": 0, "K_MAX": 1}
    
    with open(filename, 'w') as file:
        file.write(f"NDIMENSIONS=2\n")
        file.write(f"NI={ni}\n")
        file.write(f"NJ={nj}\n")
        file.write(f"NK=1\n")
        file.write(f"NPATCHES={len(patches)}\n")
        file.write("PATCH_NAME,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX\n")
        for patch_name, bounds in patches.items():
            file.write(f"{patch_name},\
                {bounds['I_MIN']},{bounds['I_MAX']},\
                {bounds['J_MIN']},{bounds['J_MAX']},\
                {bounds['K_MIN']},{bounds['K_MAX']}\n")
        
    print(f"Saved grid block {iblock} to {filename}")

