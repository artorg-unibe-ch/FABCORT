#%% !/usr/bin/env python3

Description = """
Analyse impact of the different resolution on the fabric
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__date__ = '07-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Read, Tensor


#%% Main

def Main():

    # Print time
    Time.Process(1, 'Collect data')

    # Set data paths and list files
    FabricPath = Path(__file__).parent / 'Fabric'
    Files = [F.name[:-4] for F in FabricPath.iterdir() if F.name.endswith('.fab')]

    # Array to store results (file x resolution x threhold)
    Factors = [1, 2, 4]
    BVTVs = np.zeros((len(Files)//len(Factors), len(Factors)))
    eValues = np.zeros((len(Files)//len(Factors), len(Factors), 3))
    eVectors = np.zeros((len(Files)//len(Factors), len(Factors), 3, 3))

    for i, File in enumerate(Files):

        # Get factor index
        f = Factors.index(int(File[-1]))

        # Read fabric in original resolution
        FabricFile = str(File) + '.fab'
        Fabric = Read.Fabric(FabricPath / FabricFile)
        eValues[i//len(Factors),f] = Fabric[0]
        eVectors[i//len(Factors),f] = Fabric[1]

    # Plot eigenvalues
    FName = str(Path(__file__).parent / 'Results' / 'eValues.png')
    Args = np.argsort(eValues[:,0,0])
    Colors = [(1,0,0),(0,1,0),(0,0,1)]
    Figure, Axis = plt.subplots(1,3, figsize=(7,4), dpi=200, sharex=True, sharey=True)
    for Arg in Args:
        for i in range(3):
            Axis[i].plot(Factors, eValues[Arg,:,i], color=Colors[Arg])
            Axis[i].plot(Factors, eValues[Arg,:,i], color=Colors[Arg])
            Axis[i].plot(Factors, eValues[Arg,:,i], color=Colors[Arg])
            Axis[i].set_title(f'm$_{i+1}$')
    Axis[0].set_ylabel('Eigenvalues')
    Axis[1].set_xlabel('Downscale factor (-)')
    Axis[1].set_xticks(Factors)
    plt.savefig(FName)
    plt.close(Figure)
    
    # Iterate for each file
    Time.Update(0.5, 'Plot fabric')

    # Project every fabric on each plane
    NPoints = 100
    e1s = np.zeros((len(Files)//len(Factors), len(Factors), 3, NPoints))
    e2s = np.zeros((len(Files)//len(Factors), len(Factors), 3, NPoints))
    e3s = np.zeros((len(Files)//len(Factors), len(Factors), 3, NPoints))
    for e, (eVal, eVec) in enumerate(zip(eValues, eVectors)):
        
        for f in range(len(Factors)):
        
            # Build fabric tensor
            M = Tensor.Fabric(eVal[f], eVec[f])
            for n, Normal in enumerate(np.eye(3)):

                # Project to normal plane
                eVals, eVecs = Tensor.ProjectEllipsoid(M, Normal)

                # Generate points for the ellipse
                Theta = np.linspace(0, 2 * np.pi, NPoints)
                e1 = eVals[0] * np.cos(Theta)
                e2 = eVals[1] * np.sin(Theta)

                # Rotate the points
                R = np.column_stack([eVecs[0],eVecs[1]])
                e1e2e3 = R @ np.vstack([e1,e2])

                # Store values
                e1s[e,f,n] = e1e2e3[0]
                e2s[e,f,n] = e1e2e3[1]
                e3s[e,f,n] = e1e2e3[2]



        # Plot the ellipse
    
    # Plot fabric
    FName = str(Path(__file__).parent / 'Results' / 'Fabric.png')
    Figure, Axis = plt.subplots(1,3, figsize=(10,5), dpi=200, sharex=True, sharey=True)
    for i in range(3):
        for j in range(3):
            for r, c in enumerate(Colors):
                if i == 0:
                    Axis[i].plot(e2s[j,r,0], e3s[j,r,0], color=c)
                    Axis[i].set_xlabel('e$_2$')
                    Axis[i].set_ylabel('e$_3$')
                elif i == 1:
                    Axis[i].plot(e1s[j,r,1], e3s[j,r,1], color=c)
                    Axis[i].set_xlabel('e$_1$')
                    Axis[i].set_ylabel('e$_3$')
                else:
                    Axis[i].plot(e1s[j,r,2], e2s[j,r,2], color=c)
                    Axis[i].set_xlabel('e$_1$')
                    Axis[i].set_ylabel('e$_2$')

        Axis[i].axhline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].axvline(0, color=(0,0,0), linewidth=0.5, linestyle='--')
        Axis[i].grid()
        Axis[i].set_aspect('equal')
    for c, f in zip(Colors, Factors):
        Axis[2].plot([], color=c, label=f'Downsampling factor: {f}')
    plt.legend(loc='upper center', bbox_to_anchor=(-0.75, 1.15), ncol=3)
    plt.savefig(FName)
    plt.close(Figure)
    Time.Process(0)

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
