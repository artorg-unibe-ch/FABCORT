#%% !/usr/bin/env python3

Description = """
This script performs a computational analysis to estimate
the minimum number of regions of interest (ROIs) needed to
achieve a desired accuracy threshold in homogenization
stiffness tensor computations.
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '09-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pyvista as pv
from pathlib import Path
from itertools import product
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment, curve_fit


sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time


#%% Functions

def BalancedClusters(X, NClusters, MaxIt=100, tol=1e-4):

    """
    Perform balanced clustering using linear sum assignment to enforce balanced clusters.

    Args:
        X (np.ndarray): Dataset, shape (n_samples, n_features).
        NClusters (int): Number of clusters.
        MaxIt (int): Maximum number of iterations.
        tol (float): Tolerance for convergence.
    
    Returns:
        labels (np.ndarray): Cluster assignments for each point, shape (n_samples,).
        centers (np.ndarray): Final cluster centers, shape (n_clusters, n_features).
    """

    # Step 1: Initialize cluster centers randomly
    NSamples = X.shape[0]
    Centers = X[np.random.choice(NSamples, NClusters, replace=False)]
    
    for _ in range(MaxIt):

        # Step 2: Assign points to clusters, ensuring balance
        rCenter = np.tile(Centers, (NSamples // NClusters + 1, 1))[:NSamples]
        distances = np.linalg.norm(X[:, None, :] - rCenter[None, :, :], axis=2)
        _, Indices = linear_sum_assignment(distances)
        Labels = Indices % NClusters

        # Step 3: Update cluster centers
        NewCenters = np.copy(Centers)
        for i in range(NClusters):
            NewCenters[i] = X[Labels == i].mean(axis=0)

        # Step 4: Check for convergence
        if np.all(np.linalg.norm(Centers - NewCenters, axis=1) < tol):
            break

        # Step 5: Update centers
        Centers = NewCenters

    # Return clusters
    Clusters = [np.where(Labels == i)[0] for i in range(NClusters)]

    return Clusters

def Fit(x, a, b):
    return a * np.exp(-b * x)

#%% Main

def Main():

    # Print time
    Time.Process(1, 'Get stiffness')

    # List simulations output files
    ResultsPath =Path(__file__).parent / 'Results'
    ROIs = sorted([F.name[:-4] for F in ResultsPath.iterdir() if F.name.endswith('.out')])

    # Compute stiffness of each ROI
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

    # Update time
    Time.Update(0.2, 'Compute errors')

    # Load map
    ROIMap = np.load(Path(__file__).parent / 'ROIMap.npy')

    # Iterate for each possibility of number of ROIs combined
    NSamples = 1000
    MeanErrors = []
    Sum = 0
    for NCluster in range(1,len(ROIMap)+1):

        # Cluster grid points
        Clusters = BalancedClusters(ROIMap, NCluster)

        # Compute maximum number of combinations
        MaxComb = np.prod([len(C) for C in Clusters])
        
        # If numerous possibilities, draw random ones
        if MaxComb / NSamples > 10:
            # Generate N random combitations of these clustes
            Combinations = set()
            while len(Combinations) < NSamples:
                Combination = tuple(np.random.choice(C) for C in Clusters)
                Combinations.add(Combination)

        # Else, compute all of them
        else:
            Combinations = set(product(*Clusters))
        
        # Compute mean norm error for N random combination
        Means = []
        Denominator = np.sum(Ref**2)
        for Indices in Combinations:

            # Compute average stiffness
            AStiffness = np.mean([Stiffness[i] for i in Indices], axis=0)
            
            # Compute norm errors
            Nominator = np.sum((AStiffness - Ref)**2)
            NormError = np.sqrt(Nominator/Denominator)
            
            # Store norm error for this combination
            Means.append(NormError)
        
        # Store all combinations norm errors
        MeanErrors.append(Means)

    # Update time
    Time.Update(0.8, 'Fit data')

    # Combine all data for a single fit
    X, Y = [], []
    for i, Error in enumerate(MeanErrors):
        for E in Error:
            X.append(i + 1)
            Y.append(E)

    # Fit the combined data
    Guess = [1, 1]
    Params, _ = curve_fit(Fit, X, Y, p0=Guess)

    # Plot norm errors
    Threshold = 0.05
    Figure, Axis = plt.subplots(1,1, dpi=200)
    for i, Error in enumerate(MeanErrors):
        X = np.zeros(len(Error)) + i + 1
        Axis.plot(X, Error, color=(0,0,0,0.1), linestyle='none', marker='.')
    Axis.plot([], color=(0,0,0), linestyle='none', marker='.', label='Samples')
    X = np.arange(len(MeanErrors))+1
    Y = Fit(X, *Params)
    Axis.plot(X, Y, color=(1,0,0), linewidth=2, label='Fit')
    Axis.plot(X, np.tile(Threshold,len(X)), color=(0,0,1), linestyle='--', label='Threshold')
    Axis.set_ylabel('Norm error (-)')
    Axis.set_xlabel('Number of ROIs (-)')
    Axis.set_xlim([0,len(ROIs)+1])
    Axis.set_ylim([0,1])
    plt.legend()
    plt.savefig(ResultsPath / 'NormError.png')
    plt.close(Figure)

    # Print time
    Time.Process(0, 'Analysis done')

    # Compute number of ROIs necessary
    NROIs = min(np.argwhere(Y<Threshold)[0])
    print(f'\nMinimum number of ROIs needed: {NROIs}')

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
