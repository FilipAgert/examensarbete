import numpy as np
import os

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