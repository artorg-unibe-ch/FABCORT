#%% !/usr/bin/env python3

Description = """
Select and write 3 regions of interest (ROIs) 1x1x1 mm.
These ROIs are selected in the sample with the lower
bone volume fraction (BV/TV) and downsampled to investigate
the effect of resolution.
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '07-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
from pathlib import Path
import SimpleITK as sitk

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Plot, Read, Image


#%% Main

def Main():

    # Print time
    Time.Process(1, 'Select ROIs')

    # List bone volume fraction of samples
    DataPath = Path(__file__).parents[1] / '01_Fabric/Results'
    Files = [F for F in DataPath.iterdir() if F.name.endswith('.fab')]
    BVTVs = []
    for F in Files:
        _, _, BVTV = Read.Fabric(F)
        BVTVs.append(BVTV)

    # Load sample with lower BV/TV
    ScanPath = Path(__file__).parents[1] / '00_Data'
    File = Files[BVTVs.index(min(BVTVs))].name[:-4]
    Scan = sitk.ReadImage(str(ScanPath / File) + '.mhd')

    # Select ROIs, downsample, plot and write
    Size = round(1 / Scan.GetSpacing()[0])
    Centers = np.zeros((3,3))
    Centers[0] = np.array(Scan.GetSize()) / 2   # 1st ROI around the center
    Centers[1] = np.array(Scan.GetSize()) / 3   # 2nd ROI in top third
    Centers[2] = np.array(Scan.GetSize()) / 3*2 # 3rd ROI in last third

    for c, Center in enumerate(Centers):
        X,Y,Z = np.round(Center - Size / 2)
        ROI = Scan[int(X):int(X+Size), int(Y):int(Y+Size), int(Z):int(Z+Size)]

        # Resample ROI for different resolutions
        for F in [1,2,4]:
            if F > 1:
                Resampled = Image.Resample(ROI, Factor=F)
            else:
                Resampled = Scan
            sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_3_F{F}.mhd'))

        # Print progress
        Time.Update((c+1)/3)

    # Print elapsed time
    Time.Process(0)

if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Run script
    Main()
