#%% !/usr/bin/env python3

Description = """
Compute the mean Otsu threshold across all input scans, extract regions of interest (ROIs),
and generate a parameter file for Medtool.

This script processes a folder of 3D CT scan files (in `.mhd` format). It performs the following steps:
1. Calculates the Otsu threshold for each scan and computes the mean threshold.
2. Selects a 3D grid of ROIs (2x2x4) per scan, each with a fixed voxel size based on resolution.
3. Downsamples each ROI at given factor and saves them to the disk for further medtool analysis.
4. Writes a `Parameters.csv` file containing folder names, ROI identifiers, and the mean threshold.
5. Optionally generates a 3D preview image showing the spatial layout of ROIs using PyVista.

Arguments:
    DataPath (str): Path to the directory containing input `.mhd` scan files.
    Plot (bool): Whether to generate a 3D visualization of ROI positions (default: False).

Example:
    python ComputeROIs.py ../Data_Cortical 2 True

Outputs:
    - ROIs saved in: ../Results/Cortical/ROIs/<SampleName>/
    - Parameter file: ./Medtool/Parameters.csv
    - ROI preview image (if Plot=True): ../Results/Cortical/ROIs/ROIs.png
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

def Main(DataPath, Factor, Plot):

    # List data files
    DataPath = Path(DataPath)
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Compute Otsu's threshold for all scans
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
        FPath = Path(__file__).parents[1] / 'Results' / 'Cortical' / 'ROIs' / SampleName
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
                    Resampled = Image.Resample(ROI, Factor=Factor)
                    sitk.WriteImage(Resampled,str(FPath / f'ROI_{x+1}{y+1}{z+1}.mhd'))

                    # Medtool parameter file
                    Text += f'{SampleName},{x+1}{y+1}{z+1},{Otsu},\n'
        
        # Print time
        Time.Update((f+1)/len(Files), f'Scan n°{f+1:02d} / {len(Files)}')
    Time.Process(0, 'ROIs selected')

    # Write parameter file for Medtool
    with open(Path(__file__).parent / 'Medtool' / 'Parameters.csv','w') as F:
        F.write(Text)

    if Plot:

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
        pl.screenshot(Path(__file__).parents[1] / 'Results' / 'Cortical' / 'ROIs' / 'ROIs.png', return_img=False)
        pl.close()

    return

    
if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Define default paths relative to this script's location
    ScriptDir = Path(__file__).parent
    DefaultPath = ScriptDir.parent / 'Data_Cortical'

    # Add positional arguments with defaults
    Parser.add_argument('DataPath', nargs='?', default=str(DefaultPath), help='Path to scan files')
    Parser.add_argument('Factor', nargs='?', default=int(2), help='Factor for sample coarsening')
    Parser.add_argument('Plot', nargs='?', default=False, help='Plot ROIs positions within sample')

    # Parse arguments
    Arguments = Parser.parse_args()

    # Call main function with parsed or default arguments
    Main(Arguments.DataPath, Arguments.Factor, Arguments.Plot)

#%%
