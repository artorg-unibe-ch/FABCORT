#%% !/usr/bin/env python3

Description = """
Compute and visualize the coefficient of variation (CV) of bone volume fraction for cortical and trabecular ROIs.

This script computes the spatial heterogeneity of bone volume in cortical regions by calculating the
coefficient of variation (CV) of subregion densities within each ROI. It compares these values with
existing trabecular ROI data and visualizes the relationship between CV and mean density.

Functionality:
- For each cortical ROI:
  - Reads the binarized `.mhd` image.
  - Splits the ROI into 8 subregions (2x2x2).
  - Computes the bone volume fraction (BV/TV) in each subregion.
  - Calculates the CV across subregions and stores it alongside global density.
- Loads trabecular ROI data (from `ROIsData.csv`).
- Plots CV versus density for both cortical and trabecular ROIs.
- Saves results to:
  - `Results/Cortical/CV.csv` (cortical CV and density values)
  - `Results/CV_Rho.png` (comparison plot)

Inputs:
    - Cortical ROI images: `Results/Cortical/ROIs/<Sample>/<ROI>.mhd`
    - Trabecular data CSV: `Data_Trabecular/ROIsData.csv`

Outputs:
    - `Results/Cortical/CV.csv`: Per-ROI CV and density
    - `Results/CV_Rho.png`: Plot comparing cortical and trabecular variability

Example:
    python ComputeCV.py
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
from Utils import Read, Tensor, Time

#%% Main

def Main():

    # Start timer
    Time.Process(1, 'Compute CV')

    # Define path to data
    DataPath =Path(__file__).parents[1] / 'Results/Cortical/ROIs'
    Folders = [F for F in DataPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    Idx = pd.MultiIndex.from_product([[F.name for F in Folders], ROIs], names=['Sample', 'ROI'])
    CVCort = pd.DataFrame(index=Idx, columns=['CV','Rho'])
    for f, Folder in enumerate(Folders):
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
            CVCort.loc[(Folder.name, ROI),'CV'] = np.std(Rho) / np.mean(Rho)
            CVCort.loc[(Folder.name, ROI),'Rho'] = np.sum(Array-1) / Array.size

        # Update timer
        Time.Update((f+1)/len(Folders)*0.8)

    # Save results
    CVCort.to_csv(Path(__file__).parents[1] / 'Results/Cortical/CV.csv')

    # Read trabecular CV
    DataTrab = pd.read_csv(Path(__file__).parents[1] / 'Data_Trabecular/ROIsData.csv')

    # Update timer
    Time.Update(0.8, 'Plot results')

    # Plot cv as function of density
    Figure, Axis = plt.subplots(1,1, figsize=(5,5), dpi=200)
    Axis.plot(DataTrab['BV/TV'], DataTrab['Variation Coefficient'], marker='o', fillstyle='none', color=(0,0,1), linestyle='none', label='Trabecular')
    Axis.plot(CVCort['Rho'], CVCort['CV'], marker='o', fillstyle='none', color=(1,0,0), linestyle='none', label='Cortical')
    Axis.plot([0, 1], [0.263,0.263], color=(0,0,0), linestyle='--', label='Threshold')
    Axis.set_xlabel(r'$\rho$ (-)')
    Axis.set_ylabel('Coefficient of variation (-)')
    plt.legend()
    plt.savefig(Path(__file__).parents[1] / 'Results/CV_Rho.png')
    plt.show()

    # Stop timer
    Time.Process(0, 'Done!')
    
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
