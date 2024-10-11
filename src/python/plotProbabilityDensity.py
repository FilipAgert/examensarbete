import matplotlib.pyplot as plt
from readMatrixFromFile import readMatrixFromFile

def plotMatrix(matrix):
    plt.imshow(matrix,cmap='turbo', aspect='auto')  # Use colormap 'viridis' for grid visualization
    plt.colorbar()  # Add a color bar to indicate the scale
    plt.title("Probability density")
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.show()

fileName = r'data/Fullsize'
A = readMatrixFromFile(fileName)
plotMatrix(A)