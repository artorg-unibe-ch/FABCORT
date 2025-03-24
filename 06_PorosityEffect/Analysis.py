#%% !/usr/bin/env python3

Description = """
Description
"""

__author__ = ['Mathieu Simon']
__date_created__ = '14-03-2025'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions
def PlotFabric(eValuesCort, eVectorsCort, eValuesTrab, eVectorsTrab):

    Figure, Axis = plt.subplots(1,3, figsize=(10,5), dpi=200, sharex=True, sharey=True)
    for i, Normal in enumerate(np.eye(3)):

        # Iterate for each plane
        e1s, e2s, e3s = [], [], []
        for eVal, eVec in zip(eValuesCort, eVectorsCort):

            # Build fabric tensor
            M = Tensor.Fabric(eVal, eVec)

            # Project to normal plane
            eVals, eVecs = Tensor.ProjectEllipsoid(M, Normal)

            # Generate points for the ellipse
            Theta = np.linspace(0, 2 * np.pi, 100)
            e1 = eVals[0] * np.cos(Theta)
            e2 = eVals[1] * np.sin(Theta)

            # Rotate the points
            R = np.column_stack([eVecs[0],eVecs[1]])
            e1e2e3 = R @ np.vstack([e1,e2])

            # Plot lines
            if i == 0:
                Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(0,0,1,0.2))
            elif i == 1:
                Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(0,0,1,0.2))
            else:
                Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(0,0,1,0.2))


            # Store values
            e1s.append(e1e2e3[0])
            e2s.append(e1e2e3[1])
            e3s.append(e1e2e3[2])

        Mean1 = np.mean(e1s, axis=0)
        Mean2 = np.mean(e2s, axis=0)
        Mean3 = np.mean(e3s, axis=0)
        Std1 = np.std(e1s, axis=0)
        Std2 = np.std(e2s, axis=0)
        Std3 = np.std(e3s, axis=0)

        for eVal, eVec in zip(eValuesTrab, eVectorsTrab):

            # Build fabric tensor
            M = Tensor.Fabric(eVal, eVec)

            # Project to normal plane
            eVals, eVecs = Tensor.ProjectEllipsoid(M, Normal)

            # Generate points for the ellipse
            Theta = np.linspace(0, 2 * np.pi, 100)
            e1 = eVals[0] * np.cos(Theta)
            e2 = eVals[1] * np.sin(Theta)

            # Rotate the points
            R = np.column_stack([eVecs[0],eVecs[1]])
            e1e2e3 = R @ np.vstack([e1,e2])

            # Plot lines
            if i == 0:
                Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(1,0,0,0.2))
            elif i == 1:
                Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(1,0,0,0.2))
            else:
                Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(1,0,0,0.2))


            # Store values
            e1s.append(e1e2e3[0])
            e2s.append(e1e2e3[1])
            e3s.append(e1e2e3[2])



        # Plot the ellipse
        if i == 0:
            # Axis[i].plot(Mean2, Mean3, color=(0,0,0))
            # Axis[i].fill_between(Mean2, Mean3-Std3, Mean3+Std3, color=(0.8,0.8,0.8))
            # Axis[i].fill_betweenx(Mean3, Mean2-Std2, Mean2+Std2, color=(0.8,0.8,0.8))
            Axis[i].set_xlabel('e$_2$')
            Axis[i].set_ylabel('e$_3$')

        elif i == 1:
            # Axis[i].plot(Mean1, Mean3, color=(0,0,0))
            # Axis[i].fill_between(Mean1, Mean3-Std3, Mean3+Std3, color=(0.8,0.8,0.8))
            # Axis[i].fill_betweenx(Mean3, Mean1-Std1, Mean1+Std1, color=(0.8,0.8,0.8))
            Axis[i].set_xlabel('e$_1$')
            Axis[i].set_ylabel('e$_3$')

        else:
            # Axis[i].plot(Mean1, Mean2, color=(0,0,0))
            # Axis[i].fill_between(Mean1, Mean2-Std2, Mean2+Std2, color=(0.8,0.8,0.8))
            # Axis[i].fill_betweenx(Mean2, Mean1-Std1, Mean1+Std1, color=(0.8,0.8,0.8))
            Axis[i].set_xlabel('e$_1$')
            Axis[i].set_ylabel('e$_2$')

        Axis[i].axhline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].axvline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].grid()
        Axis[i].set_aspect('equal')
    plt.show(Figure)

    return

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
    FabPath = Path(__file__).parent / 'Fabric'
    ElaPath = Path(__file__).parent / 'Elasticity'

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
    Stiffness = np.zeros((len(Files),6,6))
    for i, F in enumerate(Files):
        ComplianceMatrix = Read.ComplianceDat(ElaPath / (F+'.dat'))
        ComplianceMatrix = (ComplianceMatrix+ComplianceMatrix.T)/2
        StiffnessMatrix = np.linalg.inv(ComplianceMatrix)
        StiffnessMatrix = (StiffnessMatrix+StiffnessMatrix.T)/2
        Stiffness[i] = StiffnessMatrix


    # Read cortical bone values
    CortPath =Path(__file__).parents[1] / '05_Homogenization/Elasticity'
    Folders = [F for F in CortPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Cortical = np.zeros((len(Folders), 16, 6, 6))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get homogenization stress results
            File = open(Folder / (ROI + f'_Isotropic.out'), 'r').readlines()
            Stress = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stress[i,j] = float(File[i+4].split()[j+1])

            # Compute stiffness
            for i in range(6):
                for j in range(6):
                    Cortical[f,r,i,j] = Stress[i,j] / Strain[i]

            # Symetrize matrix
            Cortical[f,r] = 1/2 * (Cortical[f,r] + Cortical[f,r].T)

    CortFab =Path(__file__).parents[1] / '05_Homogenization/Fabric'
    Folders = [F for F in CortFab.iterdir() if F.is_dir()]
    FABCORT = np.zeros((len(Folders), 16, 3))
    CortVectors = []
    FABRho = np.zeros((len(Folders), 16))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get homogenization stress results
            Fabric = Read.Fabric(Folder / (ROI+'.fab'))
            FABCORT[f,r] = Fabric[0]
            CortVectors.append(Fabric[1])
            FABRho[f,r] = Fabric[2]

    # Plot results
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-FABRho, FABCORT[:,:,2] / FABCORT[:,:,0],
              color=(0,0,1), marker='o', linestyle='none',
              fillstyle='none')
    Axis.plot(1-BVTV, eValues[:,2] / eValues[:,0],
              color=(1,0,0), marker='o', linestyle='none',
              fillstyle='none')
    Axis.plot([], fillstyle='none', label='Cortical',
              color=(0,0,1), marker='o', linestyle='none')
    Axis.plot([], fillstyle='none', label='Trabecular',
              color=(1,0,0), marker='o', linestyle='none')
    plt.xlabel('Porosity (-)')
    plt.ylabel('Degree of Anisotropy (-)')
    plt.legend(loc='upper left')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1, dpi=192)
    Axis.plot(1-FABRho.ravel(), FABCORT[...,0].ravel(), linestyle='none', marker='o', color=(1,0,0), fillstyle='none')
    Axis.plot(1-FABRho.ravel(), FABCORT[...,1].ravel(), linestyle='none', marker='o', color=(0,0,1), fillstyle='none')
    Axis.plot(1-FABRho.ravel(), FABCORT[...,2].ravel(), linestyle='none', marker='o', color=(0,0,0), fillstyle='none')
    Axis.plot(1-BVTV, eValues[:,0], linestyle='none', marker='^', color=(1,0,0), fillstyle='none')
    Axis.plot(1-BVTV, eValues[:,1], linestyle='none', marker='^', color=(0,0,1), fillstyle='none')
    Axis.plot(1-BVTV, eValues[:,2], linestyle='none', marker='^', color=(0,0,0), fillstyle='none')
    Axis.plot([], linestyle='none', marker='o', color=(0.5,0.5,0.5), fillstyle='none', label='Cortical')
    Axis.plot([], linestyle='none', marker='^', color=(0.5,0.5,0.5), fillstyle='none', label='Trabecular')
    Axis.plot([], color=(1,0,0), label='m$_1$')
    Axis.plot([], color=(0,0,1), label='m$_2$')
    Axis.plot([], color=(0,0,0), label='m$_3$')
    Axis.set_xlabel(r'1-$\rho$ (-)')
    Axis.set_ylabel('Fabric eigenvalues (-)')
    plt.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.15))
    plt.show(Figure)

    # Plot fabric
    PlotFabric(eValues, eVectors, FABCORT, CortVectors)


    # Plot function
    X = np.linspace(0.0, 1.0)
    Y = 1 + np.exp(X*0.5)
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(X,Y)
    plt.show(Figure)

    # Fit functions to data
    def OLS(X, a, b):
        return a + b*X
    P_Cort = curve_fit(OLS, 1-FABRho.ravel(), Cortical[:,:,2,2].ravel() / Cortical[:,:,0,0].ravel())
    X_Cort = np.linspace(0.0, 0.5)
    P_Trab = curve_fit(OLS, 1-BVTV, Stiffness[:,2,2] / Stiffness[:,0,0])
    X_Trab = np.linspace(0.5, 1.0)

    Phi = np.concatenate([1-FABRho.ravel(), 1-BVTV])
    Ratios = np.concatenate([Cortical[:,:,2,2].ravel() / Cortical[:,:,0,0].ravel(), Stiffness[:,2,2] / Stiffness[:,0,0]])
    def Function(X, a, b):
        return a + np.exp(X*b)
    P = curve_fit(Function, Phi, Ratios)
    X = np.linspace(0.0, 1.0)

    Figure, Axis = plt.subplots(1,1, dpi=192)
    Axis.plot(1-FABRho, Cortical[:,:,2,2] / Cortical[:,:,0,0],
              color=(0,0,1), marker='o', linestyle='none',
              fillstyle='none')
    Axis.plot(1-BVTV, Stiffness[:,2,2] / Stiffness[:,0,0],
              color=(1,0,0), marker='o', linestyle='none',
              fillstyle='none')
    Axis.plot([], fillstyle='none', label='Cortical',
              color=(0,0,1), marker='o', linestyle='none')
    Axis.plot([], fillstyle='none', label='Trabecular',
              color=(1,0,0), marker='o', linestyle='none')
    Axis.plot(X_Cort, OLS(X_Cort, *P_Cort[0]),
              linestyle='--', color=(0,0,1))
    Axis.plot(X_Trab, OLS(X_Trab, *P_Trab[0]),
              linestyle='--', color=(1,0,0))
    # Axis.plot(X, Function(X, *P[0]),
    #           linestyle='--', color=(0,0,0), label='Fit')
    Axis.plot([], linestyle='--', color=(0,0,0), label='Fit')
    plt.xlabel('Porosity (-)')
    plt.ylabel('Stiffness Ratio (-)')
    plt.legend(loc='upper left')
    plt.show(Figure)
    
    Figure, Axis = plt.subplots(1,1, dpi=192)
    Axis.plot(1-FABRho, Cortical[:,:,0,0], color=(1,0,0), marker='o', linestyle='none', fillstyle='none')
    Axis.plot(1-FABRho, Cortical[:,:,1,1], color=(0,0,1), marker='o', linestyle='none', fillstyle='none')
    Axis.plot(1-FABRho, Cortical[:,:,2,2], color=(0,0,0), marker='o', linestyle='none', fillstyle='none')
    Axis.plot(1-BVTV, Stiffness[:,0,0], color=(1,0,0), marker='^', linestyle='none', fillstyle='none')
    Axis.plot(1-BVTV, Stiffness[:,1,1], color=(0,0,1), marker='^', linestyle='none', fillstyle='none')
    Axis.plot(1-BVTV, Stiffness[:,2,2], color=(0,0,0), marker='^', linestyle='none', fillstyle='none')
    Axis.plot([], linestyle='none', marker='o', color=(0.5,0.5,0.5), fillstyle='none', label='Cortical')
    Axis.plot([], linestyle='none', marker='^', color=(0.5,0.5,0.5), fillstyle='none', label='Trabecular')
    Axis.plot([], color=(1,0,0), label='S$_{11}$')
    Axis.plot([], color=(0,0,1), label='S$_{22}$')
    Axis.plot([], color=(0,0,0), label='S$_{33}$')
    plt.xlabel('Porosity (-)')
    plt.ylabel('Stiffness (MPa)')
    plt.legend(loc='upper right')
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
