#%% !/usr/bin/env python3

Description = """
Read results from homogenisation and fit
them to model proposed by Cowin and Yang (1999)
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

np.set_printoptions(formatter={'float_kind':'{:3}'.format}, suppress=True)

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Tensor

#%% Main

def Main():

    # Read dat file
    File  = open('Test.dat', 'r')
    Text  = File.read()
    Lines = Text.split('\n')

    # Read compliance matrix
    Index = 0
    while 'COMPLIANCE Matrix' not in Lines[Index]:
        Index += 1

    SymmetrizationFactor = float(Lines[Index+2][0])

    ComplianceMatrix = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            ComplianceMatrix[i,j] = float(Lines[Index + 4 + i].split()[j])

    # Symmetrize the matrix
    ComplianceMatrix = (ComplianceMatrix+ComplianceMatrix.T)/2

    # Read stiffness matrix
    Index = 0
    while 'STIFFNESS Matrix' not in Lines[Index]:
        Index += 1

    SymmetrizationFactor = float(Lines[Index+2][0])

    StiffnessMatrix = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            StiffnessMatrix[i,j] = float(Lines[Index + 4 + i].split()[j])

    np.round(StiffnessMatrix)
    InverseCompliance = np.linalg.inv(ComplianceMatrix)

    # Print differences
    print('Components relative difference (%)')
    print(np.round((StiffnessMatrix - InverseCompliance)/StiffnessMatrix*100,3))

    N1, N2 = np.linalg.norm(StiffnessMatrix), np.linalg.norm(InverseCompliance)
    print('Norm relative difference (%)')
    print((N1-N2)/N1*100)

    # Convert matrix to tensor
    StiffnessTensor = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i < 3 and j >= 3:
                StiffnessTensor[i,j] = StiffnessMatrix[i,j] * np.sqrt(2)
            elif i >= 3 and j < 3:
                StiffnessTensor[i,j] = StiffnessMatrix[i,j] * np.sqrt(2)
            elif i >= 3 and j >= 3:
                StiffnessTensor[i,j] = StiffnessMatrix[i,j] * 2
            else:
                StiffnessTensor[i, j] = StiffnessMatrix[i, j]

    # Symmetrize the matrix
    StiffnessTensor = (StiffnessTensor+StiffnessTensor.T)/2

    # Rewrite 6D tensor into 3D tensor
    Tensor3D = Tensor.IsoMorphism66_3333(StiffnessTensor)

    # Read fabric file
    Text = open('Test.fab','r').readlines()
    BVTV = float(Text[12].split('=')[1])
    eValues = np.array(Text[18].split(':')[1].split(),float)
    eVectors = np.zeros((3,3))
    for i in range(3):
        eVectors[i] = Text[19+i].split(':')[1].split()

    # Sort eigen values and eigen vectors
    Args = np.argsort(eValues)
    eValues = eValues[Args]
    eVectors = eVectors[Args]

    # Transform tensor into fabric coordinate system
    TransformedStiffness = np.zeros((3, 3, 3, 3))
    Q = eVectors
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    TransformedStiffness[i, j, k, l] += Q[m,i]*Q[n,j]*Q[o,k]*Q[p,l] * Tensor3D[m, n, o, p]

    # Convert back transformed stiffness into 6D matrix
    Matrix = Tensor.IsoMorphism3333_66(TransformedStiffness)

    # Project onto orthotropy
    Orthotropic = np.zeros((6,6))
    for i in range(Orthotropic.shape[0]):
        for j in range(Orthotropic.shape[1]):
            if i < 3 and j < 3:
                Orthotropic[i, j] = Matrix[i, j]
            elif i == j:
                Orthotropic[i, j] = Matrix[i, j]

    # Convert stiffness tensor to matrix
    Matrix3D = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i < 3 and j >= 3:
                Matrix3D[i,j] = Orthotropic[i,j] / np.sqrt(2)
            elif i >= 3 and j < 3:
                Matrix3D[i,j] = Orthotropic[i,j] / np.sqrt(2)
            elif i >= 3 and j >= 3:
                Matrix3D[i,j] = Orthotropic[i,j] / 2
            else:
                Matrix3D[i, j] = Orthotropic[i, j]

    # Print stiffness components
    print(f'S11 :{Matrix3D[0, 0]}')
    print(f'S22 :{Matrix3D[1, 1]}')
    print(f'S33 :{Matrix3D[2, 2]}')
    print(f'S44 :{Matrix3D[3, 3]}')
    print(f'S55 :{Matrix3D[4, 4]}')
    print(f'S66 :{Matrix3D[5, 5]}')
    print(f'S12 :{Matrix3D[0, 1]}')
    print(f'S13 :{Matrix3D[0, 2]}')
    print(f'S23 :{Matrix3D[1, 2]}')

    # Compare with data frame
    Data = pd.read_csv('Data.csv')
    S = np.zeros((6,6))
    S[0,0] = Data.loc[0,'S11']
    S[1,1] = Data.loc[0,'S22']
    S[2,2] = Data.loc[0,'S33']
    S[3,3] = Data.loc[0,'S44']
    S[4,4] = Data.loc[0,'S55']
    S[5,5] = Data.loc[0,'S66']
    S[0,1] = Data.loc[0,'S12']
    S[1,2] = Data.loc[0,'S23']
    S[2,0] = Data.loc[0,'S31']
    S[1,0] = Data.loc[0,'S21']
    S[2,1] = Data.loc[0,'S32']
    S[0,2] = Data.loc[0,'S13']
    print(S/Matrix3D)

    print(S/Orthotropic)


    # Compute corresponding compliance matrix
    Compliance = np.linalg.inv(Matrix3D)

    # Print engineering constants
    print(f'E1 :{1 / Compliance[0, 0]}')
    print(f'E2 :{1 / Compliance[1, 1]}')
    print(f'E3 :{1 / Compliance[2, 2]}')
    print(f'Mu23 :{1 / Compliance[3, 3]}')
    print(f'Mu31 :{1 / Compliance[4, 4]}')
    print(f'Mu12 :{1 / Compliance[5, 5]}')
    print(f'Nu12 :{-Compliance[0, 1] / Compliance[0, 0]}')
    print(f'Nu13 :{-Compliance[0, 2] / Compliance[0, 0]}')
    print(f'Nu23 :{-Compliance[1, 2] / Compliance[1, 1]}')

    return

if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main()

#%%
