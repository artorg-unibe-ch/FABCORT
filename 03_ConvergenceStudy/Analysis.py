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

import sys
import pickle
import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.cm import winter
from itertools import combinations
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time

#%% Functions

def Function(x, a, b, c):
    return a * np.exp(-b * x) + c

#%% Main

def Main():

    ResultsPath =Path(__file__).parent / 'Results'
    ROIs = sorted([F.name[:-4] for F in ResultsPath.iterdir() if F.name.endswith('.fab')])

    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Stiffness = np.zeros((len(ROIs),6,6))

    for r, ROI in enumerate(ROIs):

        # Get homogenization stress results
        Abaqus = open(ResultsPath / (ROI + '.out'), 'r').readlines()
        Stress = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                Stress[i,j] = float(Abaqus[i+4].split()[j+1])

        # Compute stiffness
        for i in range(6):
            for j in range(6):
                Stiffness[r,i,j] = Stress[i,j] / Strain[i]

        # Symetrize matrix
        Stiffness[r] = 1/2 * (Stiffness[r] + Stiffness[r].T)

    # Get average stiffness tensor
    Ref = np.mean(Stiffness, axis=0)

    # Load map
    with open(Path(__file__).parent / 'ROIMap', 'rb') as F:
        ROIMap = pickle.load(F)

    # Compute norm errors
    Errors = []
    Denominator = np.sum(Ref**2)
    for S in Stiffness:
        Nominator = np.sum((S - Ref)**2)
        Errors.append(np.sqrt(Nominator/Denominator))

    # Plot norm errors
    # Figure, Axis = plt.subplots(1,1)
    # Mean = []
    # for i, ROI in enumerate(ROIMap):
    #     Axis.plot(np.zeros(len(ROI))+i+1, [Errors[r-1] for r in ROI], color=(0,0,1), linestyle='none',
    #                        marker='o', fillstyle='none')
    #     Mean.append(np.mean([Errors[r-1] for r in ROI]))
    # Axis.plot(np.arange(16)+1, Mean, color=(0,0,1))
    # Axis.set_ylabel('Norm error (-)')
    # Axis.set_xlabel('Number of ROIs (-)')
    # plt.show()

    # Analyse stiffness
    Indices = [[[0,0],[1,1],[2,2]], [[0,1],[0,2],[1,2]], [[3,3],[4,4],[5,5]]]
    Labels = [['$S_{11}$','$S_{22}$','$S_{33}$'], ['$S_{12}$','$S_{13}$','$S_{23}$'], ['$S_{44}$','$S_{55}$','$S_{66}$']]

    # Figure, Axis = plt.subplots(3,3, sharex=True, sharey=True, dpi=200, figsize=(9,7))
    # for i in range(3):
    #     for j in range(3):
    #         I, J = Indices[i][j]
    #         for iROI, ROI in enumerate(ROIMap):
    #             for r in ROI:
    #                 Axis[i,j].plot(iROI+1, Stiffness[r-1][I,J] / Ref[I,J], color=(0,0,1),
    #                     linestyle='none', marker='o', fillstyle='none')
            
    #         Means, Stds = np.zeros(16), np.zeros(16)
    #         for iROI, ROI in enumerate(ROIMap):
    #             Means[iROI] = np.mean([Stiffness[r-1][I,J] / Ref[I,J] for r in ROI])
    #             Stds[iROI] = np.std([Stiffness[r-1][I,J] / Ref[I,J] for r in ROI])

    #         Axis[i,j].plot(np.arange(16)+1, Means, color=(0,0,1))
    #         Axis[i,j].fill_between(np.arange(16)+1, Means+Stds, Means-Stds, color=(0,0,1,0.2))
    #         Axis[i,j].plot([], color=(0,0,1), linestyle='none',
    #                        marker='o', fillstyle='none', label=Labels[i][j])
    #         Axis[i,j].legend(loc='upper left')
    # Axis[1,0].set_ylabel('Ratio over mean (-)')
    # Axis[2,1].set_xlabel('Number of ROIs (-)')
    # Axis[2,1].set_xticks(np.arange(2,16,2))
    # plt.show()

    AbsErrors = []
    Figure, Axis = plt.subplots(1,1)
    for i in range(3):
        for j in range(3):
            I, J = Indices[i][j]
            Color = winter((i+3*j)/8)
            MeanErrors = np.zeros(16)
            for iROI, ROI in enumerate(ROIMap):
                Error = [abs(1-Stiffness[r-1,I,J] / Ref[I,J]) for r in ROI]
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
    Axis.set_ylabel('Ratio over mean (-)')
    Axis.set_xlabel('Number of ROIs (-)')
    plt.show()

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
