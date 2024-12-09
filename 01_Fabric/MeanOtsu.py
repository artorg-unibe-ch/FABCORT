#%% !/usr/bin/env python3

Description = """
Investigates the relationships between resolution,
Otsu threshold for segmentation, and resulting BVTV
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__date__ = '04-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import SimpleITK as sitk
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time


#%% Main

def Main():

    # List data files
    DataPath = Path(__file__).parents[1] / '00_Data'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Compute Otsu's threshold for each scan
    Otsus = []
    print('Compute mean Otsu')
    Time.Process(1,f'Scan n°{1} / {len(Files)}')
    for f, File in enumerate(Files):

        # Read file
        Image = sitk.ReadImage(str(File))

        # Compute Otsu's threshold
        Otsu = sitk.OtsuThresholdImageFilter()
        Otsu.SetInsideValue(0)
        Otsu.SetOutsideValue(1)
        Seg = Otsu.Execute(Image)
        Thresold = Otsu.GetThreshold()

        # Store values
        Otsus.append(Thresold)

        # Print time
        Time.Update(f/len(Files), f'Scan n°{f+1} / {len(Files)}')
    Time.Process(0)

    # Write parameter file for Medtool
    with open(Path(__file__).parent / 'Parameters.csv','w') as F:
        F.write('$Sample,$Threshold,\n')
        for File in Files:
            Text = f'{File.name[:-4]},{int(round(np.mean(Otsus)))},\n'
            F.write(Text)
    
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
