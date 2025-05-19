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
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def PlotFabric(eValues, eVectors):

    Figure, Axis = plt.subplots(1,3, figsize=(10,5), dpi=200, sharex=True, sharey=True)
    for i, Normal in enumerate(np.eye(3)):

        # Iterate for each plane
        for eVal, eVec in zip(eValues, eVectors):

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
                Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(0,0,0,0.2))
            elif i == 1:
                Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(0,0,0,0.2))
            else:
                Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(0,0,0,0.2))


        # Plot the ellipse
        if i == 0:
            Axis[i].set_xlabel('e$_2$')
            Axis[i].set_ylabel('e$_3$')

        elif i == 1:
            Axis[i].set_xlabel('e$_1$')
            Axis[i].set_ylabel('e$_3$')

        else:
            Axis[i].set_xlabel('e$_1$')
            Axis[i].set_ylabel('e$_2$')

        Axis[i].axhline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].axvline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].grid()
        Axis[i].set_aspect('equal')
    
    # Plot average fabric tensor
    eValMean = np.mean(eValues, axis=0)
    eVecMean = np.mean(eVectors, axis=0)
    eVecMean = eVecMean / np.linalg.norm(eVecMean, axis=0)
    eVecMean[1] = np.cross(eVecMean[0], eVecMean[2])
    eVecMean[0] = np.cross(eVecMean[1], eVecMean[2])

    # Build fabric tensor
    M = Tensor.Fabric(eValMean, eVecMean)

    # Project to normal plane
    for i, Normal in enumerate(np.eye(3)):
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
            Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(0,0.75,1))
        elif i == 1:
            Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(0,0.75,1))
        else:
            Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(0,0.75,1))

    # Plot transverse isotropic fabric tensor
    eValMean = np.mean(eValues, axis=0)
    eValMean[0] = (eValMean[0] + eValMean[1]) / 2
    eValMean[1] = eValMean[0]

    # Build fabric tensor
    M = Tensor.Fabric(eValMean, eVecMean)

    # Project to normal plane
    for i, Normal in enumerate(np.eye(3)):
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
            Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(1,0,0))
        elif i == 1:
            Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(1,0,0))
        else:
            Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(1,0,0))

    Axis[1].plot([], color=(0,0,0), label='ROIs fabric')
    Axis[1].plot([], color=(0,0.75,1), label='Average fabric')
    Axis[1].plot([], color=(1,0,0), label='Transverse isotropic fabric')
    Axis[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
    plt.savefig(Path(__file__).parents[1] / '_Results/Fabric/Fabric.png')
    plt.show(Figure)

    return

#%% Main

def Main():

    # Read fabric values
    DataPath =Path(__file__).parents[1] / '_Results/Fabric'
    Folders = [F for F in DataPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    FabricDensity = np.zeros((len(Folders), len(ROIs)))
    FabricValues = np.zeros((len(Folders), len(ROIs), 3))
    FabricVectors = np.zeros((len(Folders), len(ROIs), 3, 3))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get MIL results
            Fabric = Read.Fabric(Folder / (ROI+'.fab'))
            FabricValues[f,r] = Fabric[0]
            FabricVectors[f,r] = Fabric[1]
            FabricDensity[f,r] = Fabric[2]

    # Reshape arrays
    FabricValues = np.reshape(FabricValues, (-1,3))
    FabricVectors = np.reshape(FabricVectors, (-1,3,3))
    FabricDensity = np.reshape(FabricDensity, (-1))

    # Plot fabric
    PlotFabric(FabricValues, FabricVectors)

    # Plot anisotropy as function of density
    Figure, Axis = plt.subplots(1,1, figsize=(5,5), dpi=200)
    Axis.plot(1-FabricDensity, FabricValues[:,2] / FabricValues[:,0],
              marker='o', fillstyle='none', color=(1,0,0), linestyle='none')
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Degree of anisotropy (-)')
    plt.savefig(Path(__file__).parents[1] / '_Results/Fabric/Anisotropy.png')
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
