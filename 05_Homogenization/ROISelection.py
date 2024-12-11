#%% !/usr/bin/env python3

Description = """
Read ISQ files and plot them in 3D using pyvista
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '10-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
import pyvista as pv
from pathlib import Path
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from matplotlib.patches import Rectangle

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Read, Image

#%% Main

def Main():

    # Print time
    Time.Process(1, 'Select scan')

    # List samples
    DataPath = Path(__file__).parents[1] / '00_Data'
    Files = sorted([F for F in DataPath.iterdir() if F.name.endswith('.mhd')])

    # Medtool parameter file
    Text = '$Folder,$ROI,$Threshold,\n'

    # Get mean Otsu Threshold
    ThresholdFile = Path(__file__).parents[1] / '01_Fabric' / 'Parameters.csv'
    Threshold = pd.read_csv(ThresholdFile, sep=',')
    Threshold = Threshold.loc[0,'$Threshold']

    # Define ROI size and grid
    Size = round(1/0.0065+0.5)
    Nx, Ny, Nz = 2, 2, 4
    for File in Files:

        # Create folder for ROIs
        FileParts = File.name.split('_')
        SampleName = FileParts[3] + '_' + FileParts[2] + '_' + FileParts[4][1]
        FPath = Path(__file__).parent / 'ROIs' / SampleName
        Path.mkdir(FPath, exist_ok=True)

        # Load sample
        Scan = sitk.ReadImage(str(DataPath / File))
        Array = sitk.GetArrayFromImage(Scan).T
        Shape = np.array(Array.shape)

        # Define ROI Map
        XSpacing = (Shape[0] - Nx*Size) / (Nx+1)
        YSpacing = (Shape[1] - Ny*Size) / (Ny+1)
        ZSpacing = (Shape[2] - Nz*Size) / (Nz+1)
        XPos = np.linspace(XSpacing+Size/2,Shape[0]-XSpacing-Size/2,Nx)
        YPos = np.linspace(YSpacing+Size/2,Shape[1]-YSpacing-Size/2,Ny)
        ZPos = np.linspace(ZSpacing+Size/2,Shape[2]-ZSpacing-Size/2,Nz)

        # Select, downscale and write ROIs
        for x, X in enumerate(XPos):
            for y, Y in enumerate(YPos):
                for z, Z in enumerate(ZPos):
                    X1, X2 = int(X - Size//2), int(X + Size//2)
                    Y1, Y2 = int(Y - Size//2), int(Y + Size//2)
                    Z1, Z2 = int(Z - Size//2), int(Z + Size//2)
                    ROI = Scan[X1:X2,Y1:Y2,Z1:Z2]

                    # Resample ROI and write it
                    Resampled = Image.Resample(ROI, Factor=2)
                    sitk.WriteImage(Resampled,str(FPath / f'ROI_{x+1}{y+1}{z+1}.mhd'))

                    # Medtool parameter file
                    Text += f'{SampleName},{x+1}{y+1}{z+1},{Threshold},\n'

    # Write parameter file
    with open(Path(__file__).parent / 'Parameters.csv','w') as F:
        F.write(Text)

    # Plot ROIs position with pyvista
    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle='',
                ytitle='',
                ztitle='',
                use_3d_text=False
                )
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(Array, opacity=0, cmap='bone', show_scalar_bar=False)
    for X in XPos:
        for Y in YPos:
            for Z in ZPos:
                Mesh = pv.Cube((X,Y,Z), Size, Size, Size)
                pl.add_mesh(Mesh, color=(0.0,0,0), opacity=1, style='wireframe', line_width=3)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = 120
    pl.camera.zoom(0.9)
    pl.show_bounds(**Args)
    pl.screenshot(Path(__file__).parent / 'Plots/ROIs.png', return_img=False)

    # Print elapsed time
    Time.Process(0, 'ROI selected')

if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Run script
    Main()
