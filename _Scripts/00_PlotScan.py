#%% !/usr/bin/env python3

Description = """
Script used to plot scans
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

def Main():

    # List data
    DataPath = Path(__file__).parents[1] / '_Data'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Select a file
    for i, File in enumerate(Files):
        Image = sitk.ReadImage(str(File))
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
        Time.Process(1,f'\nPlot scan {i+1}/{len(Files)}')
        pl = pv.Plotter(off_screen=True)
        actors = pl.add_volume(255 - Array,
                    cmap='bone',
                    show_scalar_bar=False)
        actors.prop.interpolation_type = 'linear'
        pl.camera_position = 'xz'
        pl.camera.roll = 0
        pl.camera.elevation = 30
        pl.camera.azimuth = -60
        pl.camera.zoom(1)
        pl.show_bounds(**Args)
        pl.screenshot(Path(__file__).parents[1] / f'_Results/Scans/{File.name[:-4]}.png', return_img=False)
        Time.Process(0,'Plot scan')

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
