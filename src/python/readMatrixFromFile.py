import numpy as np
import os
import matplotlib.pyplot as plt
def readMatrixFromFile(file_name):
    with open(file_name, 'r') as file:
        matrix = []
        for line in file:
            # Strip newline and split by '; '
            row = line.strip().split('; ')
            # Convert strings to floats and append to matrix
            matrix.append([float(element) for element in row if element])

    # Convert the list of lists to a NumPy array
    return np.array(matrix)

def plotMatrix(II,MM,matrix):



    plt.imshow(matrix,cmap='turbo', aspect='auto')  # Use colormap 'viridis' for grid visualization
    plt.colorbar()  # Add a color bar to indicate the scale
    plt.title("Probability density")
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.show()

def pdToIIMMgrid(pd, fullM):
    MIN_II = 1
    MIN_JJ = 1
    MIN_KK = 1
    MIN_LL = 1
    MIN_MM = -40
    MAX_II = 53
    MAX_JJ = 15
    MAX_KK = 15
    MAX_LL = 15
    MAX_MM = 40

    dimSize = np.array([1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM])
    MIN_MAX_DIM = np.array([[MIN_II, MIN_JJ, MIN_KK, MIN_LL, MIN_MM],[ MAX_II, MAX_JJ, MAX_KK, MAX_LL, MAX_MM]])
    if( not fullM):
        MIN_MAX_DIM[0,4] = 0
        dimSize[4] = 1 + MAX_MM

    A = np.zeros(dimSize) #Initialize zero array
    sumA = np.zeros((dimSize[0],dimSize[4]))
    for i in range(1, len(pd) + 1):
        # Get the coordinates based on the function
        c = coordFromLinearIdx(i, dimSize, MIN_MAX_DIM)
        
        # Assign the value from `pd` to `A` at the calculated coordinates
        # Adjust for zero-based indexing in Python (MATLAB is 1-based)
        A[c[0] - 1, c[1] - 1, c[2] - 1, c[3] - 1, c[4] - 1] = pd[i - 1]


    for j in range(MIN_JJ, MAX_JJ):
        for k in range(MIN_KK, MAX_KK):
            for l in range(MIN_LL, MAX_LL):
                sumA[:,:] = sumA[:,:] + A[:,j - 1,k - 1,l - 1,:]

    return sumA
    



def get_q2_val(AZ,AA,I_val):

    Aref = 240.0
    Zref = 94.0
    
    Q2vals = [4.81253,   6.82805,   8.90879,  11.06605,  13.31188,  15.65945,  18.12340,  20.72011,  23.46814,  26.38871,  29.50621,  32.84902,  36.45030,  38.39975,  40.34920,  42.47073,  44.59226,  46.91380,  49.23533,  51.79066,  54.34598,  57.17638,  60.00677,  63.12872,  66.25067,  69.48992,  72.72917,  76.04370,  79.35822,  86.12054,  92.99791,  99.97944, 107.05082, 114.19157, 121.42655, 128.69212, 136.03868, 143.41393, 150.81228, 158.20718, 165.78249, 173.09174, 180.81321, 188.18743, 195.65203, 203.06961, 210.66009, 218.12864, 225.62603, 233.36657, 240.74626, 247.79498, 254.57527]

    Q2factor = (AZ/Zref*AA/Aref)**(2/3)
    q2_scale = 0.01*3*AZ*1.2**2* AA**(2/3)/(4.0*np.pi)
    
    q2 = Q2factor*Q2vals[I_val - 1]/q2_scale

    return q2

def coordFromLinearIdx(idx, dimSize, MIN_MAX_DIM):
    tempidx = idx
    multfactor = 1
    
    # Calculate the initial multfactor
    for i in range(5):
        multfactor *= dimSize[i]
    
    # Initialize coordinate array to zeros
    coord = np.zeros(len(dimSize), dtype=int)
    
    # Loop from the last dimension to the first
    for i in range(len(dimSize) - 1, -1, -1):
        if i == 0:
            coord[i] = tempidx
        else:
            multfactor = multfactor // dimSize[i]  # Integer division
            coord[i] = (tempidx - 1) // multfactor + 1
            tempidx = (tempidx - 1) % multfactor + 1
    
    # Adjust coordinates based on MIN_MAX_DIM
    #for i in range(len(dimSize)):
    #    coord[i] = coord[i] - 1 + MIN_MAX_DIM[i, 0]
    
    return coord

def KDE(x1, x2, probability, grid_x1, grid_x2, bandwidth = 0.02):
    kde = np.zeros_like(grid_x2)
    print(len(x1))
    print(len(x2))
    print(np.shape(probability))
    for i in range(len(x1)):
        for j in range(len(x2)):
            xi = x1[i]
            yi = x2[j]
            pi = probability[i, j]  # Access the weight corresponding to (xi, yi)

            kde += pi * np.exp(-((grid_x1 - xi)**2 + (grid_x2 - yi)**2) / (2 * bandwidth**2))
    
    return kde
    
def plotKDE(X,Y,W,ax, title, markercords = None):
    #This takes into account points density when calculating probability.
    Xi=np.linspace(X.min() - 0.1, X.max()+0.1, 300)
    Yi=np.linspace(Y.min()- 0.1, Y.max()+0.1, 300)
    x_grid, y_grid = np.array(np.meshgrid(Xi, Yi))
    # Evaluate the kernel in grid positions
    kernel = KDE(X,Y,W,x_grid,y_grid,0.3)
    im = ax.contourf(x_grid, y_grid, kernel, cmap='hot', levels=100)
    ax.contour(x_grid,y_grid,kernel, colors='black',levels=10, linewidths = 0.3)

    plt.colorbar(im, label='Probability density')

    #ax.scatter(X, Y, edgecolors='black', marker='.', s=3, label='Data Points')
    ax.set_xlabel(r'$Mass assymetry \alpha$')  # Correct method for setting x-axis label
    ax.set_ylabel(r'$Elongation Q_2$')  # Correct method for setting y-axis label
    ax.set_title(title)
    #ax.set_aspect('equal','box')

    if markercords is not None: #Plot star at selected coordinates
        ax.scatter(markercords[1], markercords[0], marker='*', s = 15, color='red')


MIN_II = 1
MIN_JJ = 1
MIN_KK = 1
MIN_LL = 1
MIN_MM = -40
MAX_II = 53
MAX_JJ = 15
MAX_KK = 15
MAX_LL = 15
MAX_MM = 40

dimSize = np.array([1 + MAX_II-MIN_II, 1 + MAX_JJ-MIN_JJ, 1 + MAX_KK-MIN_KK, 1 + MAX_LL-MIN_LL, 1 + MAX_MM - MIN_MM])
MIN_MAX_DIM = np.array([[MIN_II, MIN_JJ, MIN_KK, MIN_LL, MIN_MM],[ MAX_II, MAX_JJ, MAX_KK, MAX_LL, MAX_MM]])
MIN_MAX_DIM = np.transpose(MIN_MAX_DIM)
fullM = True
AZ = 102.0
AA = 256.0

if( not fullM):
    MIN_MAX_DIM[0,4] = 0
    dimSize[4] = 1 + MAX_MM
MCoords = np.arange(MIN_MAX_DIM[4, 0], MIN_MAX_DIM[4, 1] + 1) * 0.02  # Using np.arange for inclusive range
ICoords = np.arange(MIN_MAX_DIM[0, 0], MIN_MAX_DIM[0, 1] + 1)  # Include endpoint in the range
QCoords = np.zeros(dimSize[0],dtype=float)
# Apply get_q2_val function to each element in yrang
for i in range(len(ICoords)):
    val = get_q2_val(AZ, AA, ICoords[i])  # Assuming get_q2_val is defined
    QCoords[i] = val  # Update yrang with the returned value

pd = readMatrixFromFile("data/PD-15")
grid = pdToIIMMgrid(pd, True)




fig, ax = plt.subplots(1,1)
plotKDE(MCoords,QCoords,np.transpose(grid),ax, "E = 15 MeV")
plt.show()