#%% !/usr/bin/env python3

Description = """
Investigates the relationships between resolution
and bone volume fraction using constant segmentation
threshold
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__date__ = '05-12-2024'
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
from matplotlib.cm import winter

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time


#%% Main

def Main():

    # List data files
    DataPath = Path(__file__).parent / 'ROIs'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Read mean Otsu threshold
    ThresholdFile = Path(__file__).parents[1] / '01_Fabric' / 'Parameters.csv'
    Threshold = pd.read_csv(ThresholdFile, sep=',')
    Threshold = Threshold.loc[0,'$Threshold']

    # Investigate different resolution factors
    Factors = [1,2,4]
    BVTVs = np.zeros((len(Files)//len(Factors),len(Factors)))
    print('\nSegment scans and compute bone volume fraction')
    Time.Process(1,'Scan n°1')
    for i, File in enumerate(Files):

        # Read file
        Scan = sitk.ReadImage(str(File))

        # Compute Otsu's threshold
        Filter = sitk.BinaryThresholdImageFilter()
        Filter.SetInsideValue(0)
        Filter.SetOutsideValue(1)
        Filter.SetUpperThreshold(int(Threshold))
        Seg = Filter.Execute(Scan)

        # Compute bone volume fraction
        Array = sitk.GetArrayFromImage(Seg)
        VolumeFraction = Array.sum() / Array.size

        # Store values
        f = Factors.index(int(File.name[-5]))
        BVTVs[i//len(Factors),f] = VolumeFraction

        # Update time
        Time.Update(i/(len(Files)-1),f'Scan n°{i+1} / {len(Files)}')

    # Sort by bone volume fraction for plotting
    Args = np.argsort(BVTVs[:,0])

    # Plot bone volume fraction as function of downscaling
    Time.Update(1,'Plot results')
    Figure, Axis = plt.subplots(1,1, figsize=(5.5,4.5), dpi=200)
    for Arg, Color in zip(Args, [(1,0,0),(0,1,0),(0,0,1)]):
        Axis.plot(Factors, BVTVs[Arg], color=Color)
    Axis.set_xlabel('Downscale factor (-)')
    Axis.set_xticks(Factors)
    Axis.set_ylabel('Bone volume fraction (-)')
    plt.savefig(Path(__file__).parent / 'Results' / 'BVTV_Downscaling.png')
    plt.close(Figure)

    # Write parameter file for medtool
    with open(Path(__file__).parent / 'Parameters.csv','w') as F:
        F.write('$N,$Res,$Threshold,\n')
        for i in range(3):
            for f in Factors:
                Text = f'{i+1},{f},{Threshold},\n'
                F.write(Text)
    
    # Print elapsed time
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
