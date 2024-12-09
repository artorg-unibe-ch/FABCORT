#%% !/usr/bin/env python3

Description = """
Investigate the effect of the resolution on
homogenization resutls
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '07-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Read, Tensor, Plot


#%% Main

def Main():

    # Print time
    Time.Process(1, 'Analysis')

    # Define paths and list files
    AbaqusPath =Path(__file__).parent / 'Abaqus'
    FabricPath =Path(__file__).parent / 'Fabric'
    Files = sorted([F.name[:-4] for F in AbaqusPath.iterdir() if F.name.endswith('.out')])

    # Collect data
    Factors = [1,2,4]
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Stiffness = np.zeros((len(Files)//len(Factors),len(Factors),6,6))
    BVTV = np.zeros((len(Files)//len(Factors),len(Factors)))
    for f, File in enumerate(Files):

        # Get factor index
        r = Factors.index(int(File[-1]))
        f = f//len(Factors)

        # Get homogenization stress results
        Abaqus = open(AbaqusPath / (File + '.out'), 'r').readlines()
        Stress = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                Stress[i,j] = float(Abaqus[i+4].split()[j+1])

        # Compute stiffness
        for i in range(6):
            for j in range(6):
                Stiffness[f,r,i,j] = Stress[i,j] / Strain[i]

        # Symetrize matrix
        Stiffness[f,r] = 1/2 * (Stiffness[f,r] + Stiffness[f,r].T)

        # Read fabric
        Fabric = Read.Fabric(FabricPath / (File + '.fab'))
        BVTV[f,r] = Fabric[2]

    # Plot stiffness constants
    Args = np.argsort(BVTV[:,0])
    Indices = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2),(3,3),(4,4),(5,5)]
    Figure, Axis = plt.subplots(3,3, figsize=(15,12), dpi=200, sharex=True)
    for i in range(3):
        for j in range(3):
            ii, jj = Indices[i*3+j]
            for Arg, Color in zip(Args, [(1,0,0),(0,1,0),(0,0,1)]):
                Axis[i,j].plot(Factors, Stiffness[Arg,:,ii,jj], color=Color)
            Axis[i,j].set_title(f'S$_{ii+1}$$_{jj+1}$')
            Axis[i,j].set_xticks(Factors)
    Axis[2,1].set_xlabel('Downscale factor (-)')
    Axis[1,0].set_ylabel('Stiffness value (MPa)')
    plt.savefig(Path(__file__).parent / 'Results' / 'Stiffness.png')
    plt.close(Figure)

    # Plot relationship with bone volume fraction
    S = np.shape(BVTV)
    X, Y = [], []
    for i in range(S[0]):
        for j in range(S[1]):
            for k in range(3):
                for l in range(3):
                    kk, ll = Indices[k*3+l]
                    X.append(BVTV[i,j] / BVTV[i,0])
                    Y.append(Stiffness[i,j,kk,ll] / Stiffness[i,0,kk,ll])

    FileName = str(Path(__file__).parent / 'Results' / 'BVTV_Stiffness.png')
    Parameters, R2, SE = Plot.OLS(X,Y,XLabel='BV/TV ratio (-)',
                                  YLabel='Stiffness ratio (-)',
                                  FileName=FileName, Show=False)
    
    # Print elapsed time
    Time.Process(0)
    print('\nRegression values')
    print(Parameters)

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
