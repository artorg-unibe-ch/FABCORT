#%% !/usr/bin/env python3

Description = """
Read and plot the projection of the fabric tensor
of the selected ROIs
"""

__author__ = ['Mathieu Simon']
__date_created__ = '14-03-2025'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions


#%% Main

def Main():

    # Define path to data
    DataPath =Path(__file__).parents[1] / '_Results/Cortical/ROIs'
    Folders = [F for F in DataPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    Idx = pd.MultiIndex.from_product([[F.name for F in Folders], ROIs], names=['Sample', 'ROI'])
    CV = pd.DataFrame(index=Idx, columns=['CV','Rho'])
    for Folder in Folders:
        for ROI in ROIs:
             
            # Read image
            Image = sitk.ReadImage(str(Folder / f'{ROI}.mhd'))
            Array = sitk.GetArrayFromImage(Image)
            Shape = np.array(Array.shape) // 2

            # Collect bone volume fractions
            Rho = np.zeros((2,2,2))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        Sub = Array[Shape[0]*i:Shape[0]*(i+1), Shape[1]*j:Shape[1]*(j+1), Shape[2]*k:Shape[2]*(k+1)]
                        Rho[i,j,k] = np.sum(Sub-1) / Sub.size

            # Compute coefficient of variation
            CV.loc[(Folder.name, ROI),'CV'] = np.std(Rho) / np.mean(Rho)
            CV.loc[(Folder.name, ROI),'Rho'] = np.sum(Array-1) / Array.size

    # Save results
    CV.to_csv(Path(__file__).parents[1] / '_Results/Cortical/CV.csv')

    # Plot cv as function of density
    Figure, Axis = plt.subplots(1,1, figsize=(5,5), dpi=200)
    Axis.plot(1-CV['Rho'], CV['CV'], marker='o', fillstyle='none', color=(1,0,0), linestyle='none', label='Data')
    Axis.plot([1-CV['Rho'].min(), 1-CV['Rho'].max()], [0.263,0.263], color=(0,0,0), linestyle='--')
    Axis.plot([0.5, 0.5], [CV['CV'].min(), CV['CV'].max()], color=(0,0,0), linestyle='--', label='Thresholds')
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Coefficient of variation (-)')
    plt.legend()
    plt.savefig(Path(__file__).parents[1] / '_Results/Cortical/CV_Rho.png')
    plt.show(Figure)
    
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
