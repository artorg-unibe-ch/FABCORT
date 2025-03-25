#%% !/usr/bin/env python3

Description = """
Description
"""

__author__ = ['Mathieu Simon']
__date_created__ = '24-03-2025'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import SimpleITK as sitk
from skimage import measure
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read

#%% Main

def Main():

    # Read ROIs data
    Data = pd.read_csv('ROIsData.csv')

    # Keep ROIs with CV < 0.263
    F = Data['Variation Coefficient'] < 0.263
    Files = []
    for I, Row in Data[F].iterrows():
        N = Row['ROI Number']
        S = Row['Scan']
        Files.append(f'{N}_{S}')

    # Define paths
    ROIPath = Path(__file__).parent / 'ROIs'
    FabPath = Path(__file__).parent / 'Fabric'

    # Read fabric data
    eValues = np.zeros((len(Files),3))
    BVTV = np.zeros((len(Files)))
    eVectors = []
    for i, F in enumerate(Files):
        Fabric = Read.Fabric(FabPath / (F+'.fab'))
        eValues[i] = Fabric[0]
        eVectors.append(Fabric[1])
        BVTV[i] = Fabric[2]

    # Read simulations results
    Rho = np.zeros(len(Files))
    Phi = np.zeros(len(Files))
    Sp = np.zeros(len(Files))
    for i, F in enumerate(Files):
        Image = sitk.ReadImage(str(ROIPath / (F+'_Cleaned.mhd')))
        Array = sitk.GetArrayFromImage(Image)
        Bone = Array == 1

        # Compute surface
        Vertices, Faces, _, _ = measure.marching_cubes(Bone)
        Area = measure.mesh_surface_area(Vertices, Faces)

        # Store values
        Rho[i] = Bone.sum() / Bone.size
        Phi[i] = Bone.sum() / Area
        Sp[i] = Area / Bone.size
        
    # Check bone volume fractions
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, Rho, linestyle='none', marker='o', color=(1,0,0), fillstyle='none')
    Axis.set_xlabel('Medtool')
    Axis.set_ylabel('Python')
    plt.show(Figure)


    # Read cortical bone values
    CortPath =Path(__file__).parents[1] / '05_Homogenization/ROIs'
    Folders = [F for F in CortPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    Surface = np.zeros((len(Folders), 16))
    Volume = np.zeros((len(Folders), 16))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Read image
            Image = sitk.ReadImage(str(Folder / (ROI + '.mhd')))
        
            # Compute bone mask
            Array = sitk.GetArrayFromImage(Image)
            Bone = Array > 131

            # Compute surface
            Vertices, Faces, _, _ = measure.marching_cubes(Bone)
            Area = measure.mesh_surface_area(Vertices, Faces)
            
            # Store values
            Surface[f,r] = Area
            Volume[f,r] = Bone.sum()

    # Compute surface ratio and volume fraction
    BVTV = Volume / (77*77*77)
    BVBS = Volume / Surface
    Spec = Surface / (77*77*77)

    Figure, Axis = plt.subplots(1,1, dpi=192)
    Axis.plot(BVTV.ravel(), Spec.ravel(), linestyle='none', marker='o', color=(0,0,1), fillstyle='none', label='Cortical')
    Axis.plot(Rho, Sp, linestyle='none', marker='^', color=(1,0,0), fillstyle='none', label='Trabecular')
    Axis.set_xlabel(r'$\rho$ (-)')
    Axis.set_ylabel('Bone surface / total volume (mm$^{-1}$)')
    plt.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.1))
    plt.show(Figure)

    # Fit functions to data
    def OLS(X, a, b):
        return a + b*X
    Y_Cort = curve_fit(OLS, 1-BVTV.ravel(), 1/BVBS.ravel())
    X_Cort = np.linspace(0.0, max(1-BVTV.ravel()))

    Y_Trab = curve_fit(OLS, 1-Rho, 1/Phi)
    X_Trab = np.linspace(min(1-Rho), 1)

    Figure, Axis = plt.subplots(1,1, dpi=192)
    Axis.plot(1-BVTV.ravel(), 1/BVBS.ravel(), linestyle='none', marker='o', color=(0,0,1), fillstyle='none', label='Cortical')
    Axis.plot(1-Rho, 1/Phi, linestyle='none', marker='^', color=(1,0,0), fillstyle='none', label='Trabecular')
    Axis.plot(X_Cort, Y_Cort[0][0] + Y_Cort[0][1]*X_Cort, linestyle='--', color=(0,0,1))
    Axis.plot(X_Trab, Y_Trab[0][0] + Y_Trab[0][1]*X_Trab, linestyle='--', color=(1,0,0))
    Axis.plot([], linestyle='--', color=(0,0,0), label='Fit')
    Axis.set_xlabel(r'1-$\rho$ (-)')
    Axis.set_ylabel('Bone surface / volume (mm$^{-1}$)')
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))
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
