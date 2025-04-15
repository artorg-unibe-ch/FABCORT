#%% !/usr/bin/env python3

Description = """
Compute the mean Otsu threshold for all
scans, select regions of interest (ROIs)
and write a parameter file for Medtool
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pyvista as pv
import SimpleITK as sitk
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Image


#%% Main

def Main():

    # List data files
    DataPath = Path(__file__).parents[1] / '_Data'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Compute Otsu's threshold for each scan
    Otsus = []
    print('Compute mean Otsu')
    Time.Process(1,f'Scan n°{1} / {len(Files)}')
    for f, File in enumerate(Files):

        # Read file
        Scan = sitk.ReadImage(str(File))

        if f == 0:
            Size = round(1/Scan.GetSpacing()[0]+0.5)

        # Compute Otsu's threshold
        Otsu = sitk.OtsuThresholdImageFilter()
        Otsu.SetInsideValue(0)
        Otsu.SetOutsideValue(1)
        Seg = Otsu.Execute(Scan)
        Thresold = Otsu.GetThreshold()

        # Store values
        Otsus.append(Thresold)

        # Print time
        Time.Update((f+1)/len(Files), f'Scan n°{f+1:02d} / {len(Files)}')
    Otsu = round(np.mean(Otsus))
    Time.Process(0)

    # Medtool parameter file
    Text = '$Folder,$ROI,$Threshold,\n'

    # Define ROI size and grid
    Nx, Ny, Nz = 2, 2, 4
    Time.Process(1, 'Select ROIs')
    for f, File in enumerate(Files):

        # Create folder for ROIs
        SampleName = File.name[:-4]
        FPath = Path(__file__).parents[1] / '_Results' / 'ROIs' / SampleName
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
                    Text += f'{SampleName},{x+1}{y+1}{z+1},{Otsu},\n'
        
        # Print time
        Time.Update((f+1)/len(Files), f'Scan n°{f+1:02d} / {len(Files)}')
    Time.Process(0, 'ROIs selected')

    # Write parameter file for Medtool
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
    pl.screenshot(Path(__file__).parents[1] / '_Results' / 'ROIs' / 'ROIs.png', return_img=False)

    
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
