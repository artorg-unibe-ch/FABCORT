#%% !/usr/bin/env python3

Description = """
Read ISQ files and plot them in 3D using pyvista
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '15-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import pickle
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.cm import winter
from scipy.optimize import curve_fit
from scipy.stats.distributions import t

np.set_printoptions(formatter={'float': '{: 0.1f}'.format})

#%% Functions

def GetFabric(FileName):

    Text = open(FileName,'r').readlines()
    BVTV = float(Text[12].split('=')[1])
    eValues = np.array(Text[18].split(':')[1].split(),float)
    eVectors = np.zeros((3,3))
    for i in range(3):
        eVectors[i] = Text[19+i].split(':')[1].split()

    return eValues, eVectors, BVTV

def Engineering2MandelNotation(A):

    B = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i < 3 and j >= 3:
                B[i,j] = A[i,j] * np.sqrt(2)
            elif i >= 3 and j < 3:
                B[i,j] = A[i,j] * np.sqrt(2)
            elif i >= 3 and j >= 3:
                B[i,j] = A[i,j] * 2
            else:
                B[i, j] = A[i, j]

    return B

def IsoMorphism66_3333(A):

    # Check symmetry
    Symmetry = True
    for i in range(6):
        for j in range(6):
            if not A[i,j] == A[j,i]:
                Symmetry = False
                break
    if Symmetry == False:
        print('Matrix is not symmetric!')
        return

    B = np.zeros((3,3,3,3))

    # Build 4th tensor
    B[0, 0, 0, 0] = A[0, 0]
    B[1, 1, 0, 0] = A[1, 0]
    B[2, 2, 0, 0] = A[2, 0]
    B[1, 2, 0, 0] = A[3, 0] / np.sqrt(2)
    B[2, 0, 0, 0] = A[4, 0] / np.sqrt(2)
    B[0, 1, 0, 0] = A[5, 0] / np.sqrt(2)

    B[0, 0, 1, 1] = A[0, 1]
    B[1, 1, 1, 1] = A[1, 1]
    B[2, 2, 1, 1] = A[2, 1]
    B[1, 2, 1, 1] = A[3, 1] / np.sqrt(2)
    B[2, 0, 2, 1] = A[4, 1] / np.sqrt(2)
    B[0, 1, 2, 1] = A[5, 1] / np.sqrt(2)

    B[0, 0, 2, 2] = A[0, 2]
    B[1, 1, 2, 2] = A[1, 2]
    B[2, 2, 2, 2] = A[2, 2]
    B[1, 2, 2, 2] = A[3, 2] / np.sqrt(2)
    B[2, 0, 2, 2] = A[4, 2] / np.sqrt(2)
    B[0, 1, 2, 2] = A[5, 2] / np.sqrt(2)

    B[0, 0, 1, 2] = A[0, 3] / np.sqrt(2)
    B[1, 1, 1, 2] = A[1, 3] / np.sqrt(2)
    B[2, 2, 1, 2] = A[2, 3] / np.sqrt(2)
    B[1, 2, 1, 2] = A[3, 3] / 2
    B[2, 0, 1, 2] = A[4, 3] / 2
    B[0, 1, 1, 2] = A[5, 3] / 2

    B[0, 0, 2, 0] = A[0, 4] / np.sqrt(2)
    B[1, 1, 2, 0] = A[1, 4] / np.sqrt(2)
    B[2, 2, 2, 0] = A[2, 4] / np.sqrt(2)
    B[1, 2, 2, 0] = A[3, 4] / 2
    B[2, 0, 2, 0] = A[4, 4] / 2
    B[0, 1, 2, 0] = A[5, 4] / 2

    B[0, 0, 0, 1] = A[0, 5] / np.sqrt(2)
    B[1, 1, 0, 1] = A[1, 5] / np.sqrt(2)
    B[2, 2, 0, 1] = A[2, 5] / np.sqrt(2)
    B[1, 2, 0, 1] = A[3, 5] / 2
    B[2, 0, 0, 1] = A[4, 5] / 2
    B[0, 1, 0, 1] = A[5, 5] / 2



    # Add minor symmetries ijkl = ijlk and ijkl = jikl

    B[0, 0, 0, 0] = B[0, 0, 0, 0]
    B[0, 0, 0, 0] = B[0, 0, 0, 0]

    B[0, 0, 1, 0] = B[0, 0, 0, 1]
    B[0, 0, 0, 1] = B[0, 0, 0, 1]

    B[0, 0, 1, 1] = B[0, 0, 1, 1]
    B[0, 0, 1, 1] = B[0, 0, 1, 1]

    B[0, 0, 2, 1] = B[0, 0, 1, 2]
    B[0, 0, 1, 2] = B[0, 0, 1, 2]

    B[0, 0, 2, 2] = B[0, 0, 2, 2]
    B[0, 0, 2, 2] = B[0, 0, 2, 2]

    B[0, 0, 0, 2] = B[0, 0, 2, 0]
    B[0, 0, 2, 0] = B[0, 0, 2, 0]



    B[0, 1, 0, 0] = B[0, 1, 0, 0]
    B[1, 0, 0, 0] = B[0, 1, 0, 0]

    B[0, 1, 1, 0] = B[0, 1, 0, 1]
    B[1, 0, 0, 1] = B[0, 1, 0, 1]

    B[0, 1, 1, 1] = B[0, 1, 1, 1]
    B[1, 0, 1, 1] = B[0, 1, 1, 1]

    B[0, 1, 2, 1] = B[0, 1, 1, 2]
    B[1, 0, 1, 2] = B[0, 1, 1, 2]

    B[0, 1, 2, 2] = B[0, 1, 2, 2]
    B[1, 0, 2, 2] = B[0, 1, 2, 2]

    B[0, 1, 0, 2] = B[0, 1, 2, 0]
    B[1, 0, 2, 0] = B[0, 1, 2, 0]



    B[1, 1, 0, 0] = B[1, 1, 0, 0]
    B[1, 1, 0, 0] = B[1, 1, 0, 0]

    B[1, 1, 1, 0] = B[1, 1, 0, 1]
    B[1, 1, 0, 1] = B[1, 1, 0, 1]

    B[1, 1, 1, 1] = B[1, 1, 1, 1]
    B[1, 1, 1, 1] = B[1, 1, 1, 1]

    B[1, 1, 2, 1] = B[1, 1, 1, 2]
    B[1, 1, 1, 2] = B[1, 1, 1, 2]

    B[1, 1, 2, 2] = B[1, 1, 2, 2]
    B[1, 1, 2, 2] = B[1, 1, 2, 2]

    B[1, 1, 0, 2] = B[1, 1, 2, 0]
    B[1, 1, 2, 0] = B[1, 1, 2, 0]



    B[1, 2, 0, 0] = B[1, 2, 0, 0]
    B[2, 1, 0, 0] = B[1, 2, 0, 0]

    B[1, 2, 1, 0] = B[1, 2, 0, 1]
    B[2, 1, 0, 1] = B[1, 2, 0, 1]

    B[1, 2, 1, 1] = B[1, 2, 1, 1]
    B[2, 1, 1, 1] = B[1, 2, 1, 1]

    B[1, 2, 2, 1] = B[1, 2, 1, 2]
    B[2, 1, 1, 2] = B[1, 2, 1, 2]

    B[1, 2, 2, 2] = B[1, 2, 2, 2]
    B[2, 1, 2, 2] = B[1, 2, 2, 2]

    B[1, 2, 0, 2] = B[1, 2, 2, 0]
    B[2, 1, 2, 0] = B[1, 2, 2, 0]



    B[2, 2, 0, 0] = B[2, 2, 0, 0]
    B[2, 2, 0, 0] = B[2, 2, 0, 0]

    B[2, 2, 1, 0] = B[2, 2, 0, 1]
    B[2, 2, 0, 1] = B[2, 2, 0, 1]

    B[2, 2, 1, 1] = B[2, 2, 1, 1]
    B[2, 2, 1, 1] = B[2, 2, 1, 1]

    B[2, 2, 2, 1] = B[2, 2, 1, 2]
    B[2, 2, 1, 2] = B[2, 2, 1, 2]

    B[2, 2, 2, 2] = B[2, 2, 2, 2]
    B[2, 2, 2, 2] = B[2, 2, 2, 2]

    B[2, 2, 0, 2] = B[2, 2, 2, 0]
    B[2, 2, 2, 0] = B[2, 2, 2, 0]



    B[2, 0, 0, 0] = B[2, 0, 0, 0]
    B[0, 2, 0, 0] = B[2, 0, 0, 0]

    B[2, 0, 1, 0] = B[2, 0, 0, 1]
    B[0, 2, 0, 1] = B[2, 0, 0, 1]

    B[2, 0, 1, 1] = B[2, 0, 1, 1]
    B[0, 2, 1, 1] = B[2, 0, 1, 1]

    B[2, 0, 2, 1] = B[2, 0, 1, 2]
    B[0, 2, 1, 2] = B[2, 0, 1, 2]

    B[2, 0, 2, 2] = B[2, 0, 2, 2]
    B[0, 2, 2, 2] = B[2, 0, 2, 2]

    B[2, 0, 0, 2] = B[2, 0, 2, 0]
    B[0, 2, 2, 0] = B[2, 0, 2, 0]


    # Complete minor symmetries
    B[0, 2, 1, 0] = B[0, 2, 0, 1]
    B[0, 2, 0, 2] = B[0, 2, 2, 0]
    B[0, 2, 2, 1] = B[0, 2, 1, 2]

    B[1, 0, 1, 0] = B[1, 0, 0, 1]
    B[1, 0, 0, 2] = B[1, 0, 2, 0]
    B[1, 0, 2, 1] = B[1, 0, 1, 2]

    B[2, 1, 1, 0] = B[2, 1, 0, 1]
    B[2, 1, 0, 2] = B[2, 1, 2, 0]
    B[2, 1, 2, 1] = B[2, 1, 1, 2]


    # Add major symmetries ijkl = klij
    B[0, 1, 1, 1] = B[1, 1, 0, 1]
    B[1, 0, 1, 1] = B[1, 1, 1, 0]

    B[0, 2, 1, 1] = B[1, 1, 0, 2]
    B[2, 0, 1, 1] = B[1, 1, 2, 0]


    return B

def CheckMinorSymmetry(A):
    MinorSymmetry = True
    for i in range(3):
        for j in range(3):
            PartialTensor = A[:,:, i, j]
            if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                MinorSymmetry = True
            else:
                MinorSymmetry = False
                break

    if MinorSymmetry == True:
        for i in range(3):
            for j in range(3):
                PartialTensor = np.squeeze(A[i, j,:,:])
                if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                    MinorSymmetry = True
                else:
                    MinorSymmetry = False
                    break

    return MinorSymmetry

def IsoMorphism3333_66(A):

    if CheckMinorSymmetry == False:
        print('Tensor does not present minor symmetry')
    else:

        B = np.zeros((6,6))

        B[0, 0] = A[0, 0, 0, 0]
        B[0, 1] = A[0, 0, 1, 1]
        B[0, 2] = A[0, 0, 2, 2]
        B[0, 3] = np.sqrt(2) * A[0, 0, 1, 2]
        B[0, 4] = np.sqrt(2) * A[0, 0, 2, 0]
        B[0, 5] = np.sqrt(2) * A[0, 0, 0, 1]

        B[1, 0] = A[1, 1, 0, 0]
        B[1, 1] = A[1, 1, 1, 1]
        B[1, 2] = A[1, 1, 2, 2]
        B[1, 3] = np.sqrt(2) * A[1, 1, 1, 2]
        B[1, 4] = np.sqrt(2) * A[1, 1, 2, 0]
        B[1, 5] = np.sqrt(2) * A[1, 1, 0, 1]

        B[2, 0] = A[2, 2, 0, 0]
        B[2, 1] = A[2, 2, 1, 1]
        B[2, 2] = A[2, 2, 2, 2]
        B[2, 3] = np.sqrt(2) * A[2, 2, 1, 2]
        B[2, 4] = np.sqrt(2) * A[2, 2, 2, 0]
        B[2, 5] = np.sqrt(2) * A[2, 2, 0, 1]

        B[3, 0] = np.sqrt(2) * A[1, 2, 0, 0]
        B[3, 1] = np.sqrt(2) * A[1, 2, 1, 1]
        B[3, 2] = np.sqrt(2) * A[1, 2, 2, 2]
        B[3, 3] = 2 * A[1, 2, 1, 2]
        B[3, 4] = 2 * A[1, 2, 2, 0]
        B[3, 5] = 2 * A[1, 2, 0, 1]

        B[4, 0] = np.sqrt(2) * A[2, 0, 0, 0]
        B[4, 1] = np.sqrt(2) * A[2, 0, 1, 1]
        B[4, 2] = np.sqrt(2) * A[2, 0, 2, 2]
        B[4, 3] = 2 * A[2, 0, 1, 2]
        B[4, 4] = 2 * A[2, 0, 2, 0]
        B[4, 5] = 2 * A[2, 0, 0, 1]

        B[5, 0] = np.sqrt(2) * A[0, 1, 0, 0]
        B[5, 1] = np.sqrt(2) * A[0, 1, 1, 1]
        B[5, 2] = np.sqrt(2) * A[0, 1, 2, 2]
        B[5, 3] = 2 * A[0, 1, 1, 2]
        B[5, 4] = 2 * A[0, 1, 2, 0]
        B[5, 5] = 2 * A[0, 1, 0, 1]

        return B
    
def TransformTensor(A,OriginalBasis,NewBasis):

    # Build change of coordinate matrix
    O = OriginalBasis
    N = NewBasis

    Q = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            Q[i,j] = np.dot(O[i,:],N[j,:])

    if A.size == 36:
        A4 = IsoMorphism66_3333(A)

    elif A.size == 81 and A.shape == (3,3,3,3):
        A4 = A

    TransformedA = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    TransformedA[i, j, k, l] += Q[m,i]*Q[n,j]*Q[o,k]*Q[p,l] * A4[m, n, o, p]
    if A.size == 36:
        TransformedA = IsoMorphism3333_66(TransformedA)

    return TransformedA

def Mandel2EngineeringNotation(A):

    B = np.zeros((6,6))

    for i in range(6):
        for j in range(6):

            if i < 3 and j >= 3:
                B[i,j] = A[i,j] / np.sqrt(2)

            elif i >= 3 and j < 3:
                B[i,j] = A[i,j] / np.sqrt(2)

            elif i >= 3 and j >= 3:
                B[i,j] = A[i,j] / 2

            else:
                B[i, j] = A[i, j]

    return B

def Function(x, a, b, c):
    return a * np.exp(-b * x) + c

#%% Main

def Main():

    ResultsPath =Path(__file__).parent / 'Results'
    ROIs = sorted([F.name[:-4] for F in ResultsPath.iterdir() if F.name.endswith('.fab')])

    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    S = np.zeros((len(ROIs),6,6))
    eValues = np.zeros((len(ROIs),3))
    eVectors = np.zeros((len(ROIs),3,3))
    BVTV = np.zeros(len(ROIs))

    for r, ROI in enumerate(ROIs):

        # Get fabric info
        Fabric = GetFabric(ResultsPath / (ROI + '.fab'))
        eValues[r] = Fabric[0][::-1]
        eVectors[r] = Fabric[1][::-1]
        BVTV[r] = Fabric[2]

        # Get homogenization stress results
        try:
            Abaqus = open(ResultsPath / (ROI + '.out'), 'r').readlines()
        except:
            continue
        Stress = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                Stress[i,j] = float(Abaqus[i+4].split()[j+1])

        # Compute stiffness
        Stiffness = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                Stiffness[i,j] = Stress[i,j] / Strain[i]

        # Symetrize matrix
        Stiffness = 1/2 * (Stiffness + Stiffness.T)

        # Transform stiffness to fabric coordinate system
        Mandel = Engineering2MandelNotation(Stiffness)
        FStiffness = TransformTensor(Mandel, np.eye(3), eVectors[r])
        Stiffness = Mandel2EngineeringNotation(FStiffness)

        # Project onto orthotropy
        for i in range(6):
            for j in range(6):
                if i >= 3 or j >= 3:
                    if i != j:
                        Stiffness[i,j] = 0

        # Store stiffness results
        S[r] = Stiffness

    # Get average stiffness tensor
    Ref = S.mean(axis=0)

    # Load map
    with open(Path(__file__).parent / 'ROIMap', 'rb') as F:
        ROIMap = pickle.load(F)

    # Compute norm errors
    Errors = []
    Denominator = np.sum(Ref**2)
    for Stiffness in S:
        Nominator = np.sum((Stiffness - Ref)**2)
        Errors.append(np.sqrt(Nominator/Denominator))

    # Plot norm errors
    Figure, Axis = plt.subplots(1,1)
    Mean = []
    for i, ROI in enumerate(ROIMap):
        Axis.plot(np.zeros(len(ROI))+i+1, [Errors[r-1] for r in ROI], color=(0,0,1), linestyle='none',
                           marker='o', fillstyle='none')
        Mean.append(np.mean([Errors[r-1] for r in ROI]))
    Axis.plot(np.arange(16)+1, Mean, color=(0,0,1))
    Axis.set_ylabel('Norm error (-)')
    Axis.set_xlabel('Number of ROIs (-)')
    plt.show(Figure)

    # Analyse stiffness
    Indices = [[[0,0],[1,1],[2,2]], [[0,1],[0,2],[1,2]], [[3,3],[4,4],[5,5]]]
    Labels = [['$S_{11}$','$S_{22}$','$S_{33}$'], ['$S_{12}$','$S_{13}$','$S_{23}$'], ['$S_{44}$','$S_{55}$','$S_{66}$']]

    Figure, Axis = plt.subplots(3,3, sharex=True, sharey=True, dpi=200, figsize=(9,7))
    for i in range(3):
        for j in range(3):
            I, J = Indices[i][j]
            for iROI, ROI in enumerate(ROIMap):
                for r in ROI:
                    Axis[i,j].plot(iROI+1, S[r-1][I,J] / Ref[I,J], color=(0,0,1),
                        linestyle='none', marker='o', fillstyle='none')
            
            Means, Stds = np.zeros(16), np.zeros(16)
            for iROI, ROI in enumerate(ROIMap):
                Means[iROI] = np.mean([S[r-1][I,J] / Ref[I,J] for r in ROI])
                Stds[iROI] = np.std([S[r-1][I,J] / Ref[I,J] for r in ROI])

            Axis[i,j].plot(np.arange(16)+1, Means, color=(0,0,1))
            Axis[i,j].fill_between(np.arange(16)+1, Means+Stds, Means-Stds, color=(0,0,1,0.2))
            Axis[i,j].plot([], color=(0,0,1), linestyle='none',
                           marker='o', fillstyle='none', label=Labels[i][j])
            Axis[i,j].legend(loc='upper left')
    # Axis[0,0].set_ylim([0.5,1.5])
    Axis[1,0].set_ylabel('Ratio over mean (-)')
    Axis[2,1].set_xlabel('Number of ROIs (-)')
    Axis[2,1].set_xticks(np.arange(2,16,2))
    plt.show(Figure)

    AbsErrors = []
    Figure, Axis = plt.subplots(1,1)
    for i in range(3):
        for j in range(3):
            I, J = Indices[i][j]
            Color = winter((i+3*j)/8)
            MeanErrors = np.zeros(16)
            for iROI, ROI in enumerate(ROIMap):
                Error = [abs(1-S[r-1,I,J] / Ref[I,J]) for r in ROI]
                MeanErrors[iROI] = np.mean(Error)
            AbsError = np.abs(1-MeanErrors)
            AbsErrors.append(AbsError)
            Axis.plot(np.arange(16)+1, AbsError, color=Color, label=Labels[i][j])
    
    # Combine all data for a single fit
    X = np.tile(np.arange(16)+1, len(AbsErrors))  # Repeat x for each row
    Y = np.array(AbsErrors).flatten()  # Flatten the 2D array into a single 1D array

    # Fit the combined data
    Guess = [1, 1, 0]
    Params, _ = curve_fit(Function, X, Y, p0=Guess)

    Axis.plot(np.arange(16)+1, Function(np.arange(16)+1, *Params), label='Fit', color=(1,0,0), linewidth=2)
    Axis.legend(loc='right', bbox_to_anchor=(1.2,0.5))
    # Axis.set_ylim([-0.01,1])
    Axis.set_ylabel('Ratio over mean (-)')
    Axis.set_xlabel('Number of ROIs (-)')
    plt.show(Figure)

    # Read original fabric
    DataPath = Path(__file__).parents[1] / '01_Fabric/Results'
    Files = [F for F in DataPath.iterdir()]
    BVTVs = []
    for F in Files:
        _, _, BVTV = GetFabric(F)
        BVTVs.append(BVTV)
    FabricFile = Files[BVTVs.index(min(BVTVs))].name[:-4]
    Fabric = GetFabric(DataPath / (FabricFile + '.fab'))

    Figure, Axis = plt.subplots(1,1)
    Mean1, Mean2, Mean3 = [], [], []
    for r, ROI in enumerate(ROIMap):
        X = np.ones(len(ROI)) + r
        Y = np.array([eValues[R-1][2] / Fabric[0][0] for R in ROI])
        Mean1.append(np.mean(Y))
        Axis.plot(X, Y, linestyle='none', marker='o', fillstyle='none', color=(0,0,1))
        Y = np.array([eValues[R-1][1] / Fabric[0][1] for R in ROI])
        Mean2.append(np.mean(Y))
        Axis.plot(X, Y, linestyle='none', marker='o', fillstyle='none', color=(1,0,1))
        Y = np.array([eValues[R-1][0] / Fabric[0][2] for R in ROI])
        Mean3.append(np.mean(Y))
        Axis.plot(X, Y, linestyle='none', marker='o', fillstyle='none', color=(1,0,0))
    for Mean, C, L in zip([Mean1, Mean2, Mean3], [(0,0,1), (1,0,1), (1,0,0)], ['m$_1$', 'm$_2$', 'm$_3$']):
        Axis.plot(np.arange(1,17), Mean, color=C, label=L)
    plt.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.1))
    plt.show()


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
