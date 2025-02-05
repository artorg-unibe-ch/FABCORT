#%% !/usr/bin/env python3

Description = """
Plot the relationship between the bone volume fraction
obtained using the mean Otsu's threshold and the one
of the original publication with CT Analyzer. Additionally,
Plot the distribution of the fabric in the image coordinate
system.
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__date__ = '06-12-2024'
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
from Utils import Time, Plot, Read, Tensor

#%% Main

def Main():

    Time.Process(1, 'Compare BV/TV')

    # List data files
    FilePath = Path(__file__).parent / 'Results'
    Files = sorted([F for F in Path.iterdir(FilePath) if F.name.endswith('.fab')])

    # Read data of CT Analyzer
    DataPath = Path(__file__).parents[1] / '00_data'
    Data = pd.read_excel(DataPath / 'Cij.xls')

    # Collect BV/TV and fabric
    X, Y = [], []
    eValues, eVectors = [], []
    for File in Files:

        # Read file
        Fabric = Read.Fabric(File)

        # Find sample in data frame
        FileParts = File.name.split('_')
        SampleName = FileParts[3] + '_' + FileParts[2]
        Filter1 = Data['sample'] == SampleName
        Filter2 = Data['Octant'] == FileParts[4][1]
        Porosity = Data[Filter1 & Filter2]['porosity'].values[0]

        # Store values
        X.append(1 - Porosity)
        Y.append(Fabric[2])
        eValues.append(Fabric[0])
        eVectors.append(Fabric[1])

    # Plot relation to published data        
    FName = str(Path(__file__).parent / 'Results' / 'BVTV_OLS.png')
    OLS = Plot.OLS(X, Y, XLabel='CTAnalyser', YLabel='Otsu Threshold', FileName=FName, Show=False)
    Time.Process(0)

    # Iterate for each file
    FName = str(Path(__file__).parent / 'Results' / 'Fabric.png')
    Time.Process(1, 'Plot fabric')
    Figure, Axis = plt.subplots(1,3, figsize=(10,5), dpi=200, sharex=True, sharey=True)
    for i, Normal in enumerate(np.eye(3)):

        # Iterate for each plane
        e1s, e2s, e3s = [], [], []
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
