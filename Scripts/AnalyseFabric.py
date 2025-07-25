#%% !/usr/bin/env python3

Description = """
Read and visualize fabric tensor projections from cortical and trabecular ROIs.

This script reads precomputed fabric tensors (stored in `.fab` files) for both cortical and trabecular
regions of interest (ROIs), computes 2D projections of the fabric ellipsoids, and plots them on their
principal anatomical planes. The visualization compares anisotropy patterns across both bone types.

Functionality:
- Reads MIL-based fabric tensor eigenvalues, eigenvectors, and density for each ROI:
  - Cortical data from: Results/Cortical/Fabric/ROI_xxyyzz.fab
  - Trabecular data from: Data_Trabecular/*.fab
- Projects each 3D fabric tensor onto the XY, YZ, and XZ planes using ellipsoid projection.
- Generates a figure with fabric ellipses:
  - Blue for trabecular ROIs
  - Red for cortical ROIs
- Saves the figure to `Results/Fabric.png`.

Dependencies:
    - '.fab' files are produced from medtool.
    - `Utils/Read.py`, `Utils/Tensor.py`, and `Utils/Time.py` must be available in the parent directory.

Outputs:
    - A comparison plot of fabric tensor projections: `Results/Fabric.png`

Example:
    python PlotFabric.py
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
from Utils import Read, Tensor, Time

#%% Functions

def PlotFabric(eValuesTrab, eVectorsTrab, eValuesCort, eVectorsCort):

    Figure, Axis = plt.subplots(1,3, figsize=(10,5), dpi=200, sharex=True, sharey=True)
    for i, Normal in enumerate(np.eye(3)):

        # Plot trabecular eigen values
        for eVal, eVec in zip(eValuesTrab, eVectorsTrab):

            # Build fabric tensor
            eVec = np.eye(3)
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

        # Plot cortical eigen values
        for eVal, eVec in zip(eValuesCort, eVectorsCort):

            # Build fabric tensor
            eVec = np.eye(3)
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

        # Set axis labels
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
    
    # # Plot average fabric tensor
    # eValMean = np.mean(eValues, axis=0)
    # eVecMean = np.mean(eVectors, axis=0)
    # eVecMean = eVecMean / np.linalg.norm(eVecMean, axis=0)
    # eVecMean[1] = np.cross(eVecMean[0], eVecMean[2])
    # eVecMean[0] = np.cross(eVecMean[1], eVecMean[2])

    # # Build fabric tensor
    # M = Tensor.Fabric(eValMean, eVecMean)

    # # Project to normal plane
    # for i, Normal in enumerate(np.eye(3)):
    #     eVals, eVecs = Tensor.ProjectEllipsoid(M, Normal)

    #     # Generate points for the ellipse
    #     Theta = np.linspace(0, 2 * np.pi, 100)
    #     e1 = eVals[0] * np.cos(Theta)
    #     e2 = eVals[1] * np.sin(Theta)

    #     # Rotate the points
    #     R = np.column_stack([eVecs[0],eVecs[1]])
    #     e1e2e3 = R @ np.vstack([e1,e2])

    #     # Plot lines
    #     if i == 0:
    #         Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(0,0.75,1))
    #     elif i == 1:
    #         Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(0,0.75,1))
    #     else:
    #         Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(0,0.75,1))

    # # Plot transverse isotropic fabric tensor
    # eValMean = np.mean(eValues, axis=0)
    # eValMean[0] = (eValMean[0] + eValMean[1]) / 2
    # eValMean[1] = eValMean[0]

    # # Build fabric tensor
    # M = Tensor.Fabric(eValMean, eVecMean)

    # # Project to normal plane
    # for i, Normal in enumerate(np.eye(3)):
    #     eVals, eVecs = Tensor.ProjectEllipsoid(M, Normal)

    #     # Generate points for the ellipse
    #     Theta = np.linspace(0, 2 * np.pi, 100)
    #     e1 = eVals[0] * np.cos(Theta)
    #     e2 = eVals[1] * np.sin(Theta)

    #     # Rotate the points
    #     R = np.column_stack([eVecs[0],eVecs[1]])
    #     e1e2e3 = R @ np.vstack([e1,e2])

    #     # Plot lines
    #     if i == 0:
    #         Axis[i].plot(e1e2e3[1], e1e2e3[2], color=(1,0,0))
    #     elif i == 1:
    #         Axis[i].plot(e1e2e3[0], e1e2e3[2], color=(1,0,0))
    #     else:
    #         Axis[i].plot(e1e2e3[0], e1e2e3[1], color=(1,0,0))

    Axis[1].plot([], color=(1,0,0), label='Cortical')
    Axis[1].plot([], color=(0,0,1), label='Trabecular')
    Axis[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    plt.savefig(Path(__file__).parents[1] / 'Results/Fabric.png')
    plt.show()

    return

#%% Main

def Main():

    # Start timer
    Time.Process(1,'Get cortical data')

    # Read cortical fabric data
    DataPath = Path(__file__).parents[1] / 'Results/Cortical/Fabric'
    Folders = [F for F in DataPath.iterdir() if F.is_dir()]

    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    FabricDensityCort = np.zeros((len(Folders), len(ROIs)))
    FabricValuesCort = np.zeros((len(Folders), len(ROIs), 3))
    FabricVectorsCort = np.zeros((len(Folders), len(ROIs), 3, 3))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get MIL results
            Fabric = Read.Fabric(Folder / (ROI+'.fab'))
            FabricValuesCort[f,r] = Fabric[0]
            FabricVectorsCort[f,r] = Fabric[1]
            FabricDensityCort[f,r] = Fabric[2]

    # Reshape arrays
    FabricValuesCort = np.reshape(FabricValuesCort, (-1,3))
    FabricVectorsCort = np.reshape(FabricVectorsCort, (-1,3,3))
    FabricDensityCort = np.reshape(FabricDensityCort, (-1))

    # Update timer
    Time.Update(0.5,'Get trabecular data')

    # Read trabecular fabric data
    DataPath = Path(__file__).parents[1] / 'Data_Trabecular'
    Files = [F for F in DataPath.iterdir() if F.name.endswith('.fab')]

    FabricDensityTrab = np.zeros((len(Files)))
    FabricValuesTrab = np.zeros((len(Files), 3))
    FabricVectorsTrab = np.zeros((len(Files), 3, 3))
    for f, File in enumerate(Files):

        # Get MIL results
        Fabric = Read.Fabric(File)
        FabricValuesTrab[f] = Fabric[0]
        FabricVectorsTrab[f] = Fabric[1]
        FabricDensityTrab[f] = Fabric[2]

    # # Plot anisotropy as function of density
    # Figure, Axis = plt.subplots(1,1, figsize=(5,5), dpi=200)
    # Axis.plot(FabricDensityTrab, FabricValuesTrab[:,2] / FabricValuesTrab[:,0],
    #           marker='o', fillstyle='none', color=(0,0,1), linestyle='none', label='Trabecular')
    # Axis.plot(FabricDensityCort, FabricValuesCort[:,2] / FabricValuesCort[:,0],
    #           marker='o', fillstyle='none', color=(1,0,0), linestyle='none', label='Cortical')
    # Axis.set_xlabel(r'$\rho$ (-)')
    # Axis.set_ylabel('Degree of anisotropy (-)')
    # plt.legend()
    # plt.show()

    # Plot fabric
    Time.Update(0.8,'Plot fabric')
    PlotFabric(FabricValuesTrab, FabricVectorsTrab, FabricValuesCort, FabricVectorsCort)
    
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
