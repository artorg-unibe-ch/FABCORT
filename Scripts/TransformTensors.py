#%% !/usr/bin/env python3

Description = """
Read and process homogenization results to extract stiffness tensors and transform them into fabric coordinate system.

This script processes simulation outputs from Abaqus homogenization studies on cortical bone ROIs.
It reads homogenized stiffness tensors (from `.out` files) and fabric tensor information (from `.fab` files),
converts the data into Mandel notation, and transforms each stiffness tensor into the fabric coordinate system.
Results are saved for both Isotropic and Transverse material assumptions in a format suitable for model fitting.

Functionality:
- Collects stiffness tensors from Abaqus output files (`*.out`)
- Computes engineering stiffness matrices from stress and predefined strain values
- Converts engineering matrices to Mandel notation and rotates them into fabric-aligned coordinates
- Extracts corresponding fabric eigenvalues and apparent density from `.fab` files
- Assembles and exports two datasets:
  - `Isotropic.csv` — assuming isotropic material inputs
  - `Transverse.csv` — assuming transversely isotropic inputs

Inputs:
    - `Results/Cortical/Homogenisation/<Sample>/<ROI>_Isotropic.out`
    - `Results/Cortical/Homogenisation/<Sample>/<ROI>_Transverse.out`
    - `Results/Cortical/Fabric/<Sample>/<ROI>.fab`

Outputs:
    - `Results/Cortical/Isotropic.csv`: 9x9 stiffness tensor components, density, and fabric eigenvalues
    - `Results/Cortical/Transverse.csv`: Same format, for transverse case

Example:
    python ExtractStiffness.py
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

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor, Time

#%% Main

def Main():

    # Start timer
    Time.Process(1, 'Collect data')

    # Define paths
    ResultsPath =Path(__file__).parents[1] / 'Results/Cortical/Homogenisation'
    FabricPath = Path(__file__).parents[1] / 'Results/Cortical/Fabric'
    Folders = [Folder.name for Folder in FabricPath.iterdir() if Folder.is_dir()]

    # List simulations output files
    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    # Collect data
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Isotropic = np.zeros((len(Folders), len(ROIs), 6, 6))
    Transverse = np.zeros((len(Folders), len(ROIs), 6, 6))
    eVectors = np.zeros((len(Folders), len(ROIs), 3, 3))
    eValues = np.zeros((len(Folders), len(ROIs), 3))
    Rho = np.zeros((len(Folders), len(ROIs)))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get MIL results
            FabricData = Read.Fabric(FabricPath / Folder / (ROI+'.fab'))
            eVectors[f,r] = FabricData[1]
            eValues[f,r] = FabricData[0]
            Rho[f,r] = FabricData[2]

            # Get simulation results
            for S, sType in zip([Isotropic, Transverse],['Isotropic','Transverse']):

                # Get homogenization stress results
                File = open(ResultsPath / Folder / (ROI + f'_{sType}.out'), 'r').readlines()
                Stress = np.zeros((6,6))
                for i in range(6):
                    for j in range(6):
                        Stress[i,j] = float(File[i+4].split()[j+1])

                # Compute stiffness
                for i in range(6):
                    for j in range(6):
                        S[f,r,i,j] = Stress[i,j] / Strain[i]

                # Symetrize matrix
                S[f,r] = 1/2 * (S[f,r] + S[f,r].T)

        # Update timer
        Time.Update((f+1)/len(Folders)*0.5)

    # Reshape arrays
    Isotropic = np.reshape(Isotropic, (-1, 6, 6))
    Transverse = np.reshape(Transverse, (-1, 6, 6))
    eVectors = np.reshape(eVectors, (-1, 3, 3))
    eValues = np.reshape(eValues, (-1, 3))
    Rho = np.reshape(Rho, (-1))

    # Convert to Mandel notation and transform into fabric coordinate system
    I = np.eye(3)
    for i, (Vectors, Stiffness) in enumerate(zip(eVectors, Isotropic)):
        Q = np.array(Vectors)
        Mandel = Tensor.Engineering2MandelNotation(Stiffness)
        Mandel4 = Tensor.IsoMorphism66_3333(Mandel)
        TS4 = Tensor.TransformTensor(Mandel4, I, Q)
        Isotropic[i] = Tensor.IsoMorphism3333_66(TS4)

        # Update timer
        Time.Update((i+1)/len(Rho)*0.25+0.5)

    for i, (Vectors, Stiffness) in enumerate(zip(eVectors, Transverse)):
        Q = np.array(Vectors)
        Mandel = Tensor.Engineering2MandelNotation(Stiffness)
        Mandel4 = Tensor.IsoMorphism66_3333(Mandel)
        TS4 = Tensor.TransformTensor(Mandel4, I, Q)
        Transverse[i] = Tensor.IsoMorphism3333_66(TS4)

        # Update timer
        Time.Update((i+1)/len(Rho)*0.25+0.75)

    # Store results to into dataframe
    Index = np.arange(len(Rho))
    Columns = [f'S{i+1}{j+1}' for i in range(3) for j in range(3)]
    Columns.extend([f'S{i}{i}' for i in [4,5,6]])
    Columns.extend(['BV/TV','m1','m2','m3'])
    IsoData = pd.DataFrame(index=Index, columns=Columns)
    TraData = pd.DataFrame(index=Index, columns=Columns)
    for Data, S in zip([IsoData, TraData], [Isotropic, Transverse]):
        for i in range(3):
            for j in range(3):
                Data[f'S{i+1}{j+1}'] = S[:,i,j]
        for i in [3,4,5]:
            Data[f'S{i+1}{i+1}'] = S[:,i,i]
        Data['BV/TV'] = Rho
        for i in range(3):
            Data[f'm{i+1}'] = eValues[:,i]

    # Save to CSV
    IsoData.to_csv(Path(__file__).parents[1] / 'Results/Cortical/Isotropic.csv', index=False)
    TraData.to_csv(Path(__file__).parents[1] / 'Results/Cortical/Transverse.csv', index=False)
    
    # Stop timer
    Time.Process(0, 'Saved to csv')

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
