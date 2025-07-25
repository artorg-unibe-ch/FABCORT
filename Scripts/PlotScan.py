#%% !/usr/bin/env python3

Description = """
Plot a 3D CT scan and save a volumetric rendering as an image.

This script loads a 3D image file (e.g., NIfTI, DICOM) using SimpleITK,
converts it to a NumPy array, and visualizes it using PyVista. The output
is a 3D rendering of the scan saved as an image file. Useful for quickly
previewing scan volumes without manual inspection. Note that the image is
inverted before plotting (255-Image).

Usage:
    python script_name.py <ScanFile> <OutputFile>

Arguments:
    ScanFile (str): Path to the input 3D scan file.
    OutputFile (str): Path to save the rendered screenshot.

Example:
    python PlotScan.py scan.mhd output.png
"""

__author__ = ['Mathieu Simon']
__date_created__ = '20-11-2024'
__license__ = 'MIT'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pyvista as pv
from pathlib import Path
import SimpleITK as sitk

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time


#%% Main

def Main(ScanFile, OutputFile):

    # Start timer
    Time.Process(1, 'Plot scan')

    Image = sitk.ReadImage(ScanFile)
    Array = sitk.GetArrayFromImage(Image).T
    Spacing = np.array(Image.GetSpacing())[::-1]
    Shape = np.round(np.array(Array.shape) * Spacing,1)

    Args = dict(font_family='times', 
                font_size=30,
                location='outer',
                show_xlabels=False,
                show_ylabels=False,
                show_zlabels=False,
                all_edges=True,
                fmt='%i',
                xtitle=f'{Shape[0]} mm',
                ytitle=f'{Shape[1]} mm',
                ztitle=f'{Shape[2]} mm',
                use_3d_text=False
                )

    # Plot using pyvista
    pl = pv.Plotter(off_screen=True)
    actors = pl.add_volume(255 - Array,cmap='bone',show_scalar_bar=False)
    actors.prop.interpolation_type = 'linear'
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = -60
    pl.camera.zoom(1)
    pl.show_bounds(**Args)
    pl.screenshot(OutputFile, return_img=False)
    pl.close()

    # Stop timer
    Time.Process(0,'Scan plotted')

    return


if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Define default paths relative to this script's location
    ScriptDir = Path(__file__).parent
    DefaultScan = ScriptDir.parent / 'Data_Cortical' / '2009_213_L.mhd'
    DefaultOutput = ScriptDir.parent / 'Results' / 'Tests' / '2009_213_L.png'

    # Add positional arguments with defaults
    Parser.add_argument('ScanFile', nargs='?', default=str(DefaultScan), help='Path to input scan file')
    Parser.add_argument('OutputFile', nargs='?', default=str(DefaultOutput), help='Path to output image file')

    # Parse arguments
    Arguments = Parser.parse_args()

    # Call main function with parsed or default arguments
    Main(Arguments.ScanFile, Arguments.OutputFile)
