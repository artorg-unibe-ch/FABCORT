#%% !/usr/bin/env python3

Description = """
Script used to plot a scan and its fabric
"""

__author__ = ['Mathieu Simon']
__date_created__ = '20-11-2024'
__date__ = '05-12-2024'
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
from Utils import Time, Read

#%% Functions

def PlotFabricROI(ROI:np.array, eValues:np.array, eVectors:np.array, FileName:Path) -> None:

    """
    Plots a 3D ellipsoid representing a region of interest (ROI) with scaling based on the
    eigenvalues and eigenvectors provided. The ellipsoid is overlaid on a binary structure mesh,
    and the plot is generated with the ability to visualize the MIL (Mean Intercept Length) values.

    Parameters:
    -----------
    ROI (3D array): A 3D binary array representing the region of interest (ROI).
        
    eValues (1D array): A 1D array containing the eigenvalues of the fabric.
        
    eVectors (3D array) : A 2D array (shape: 3x3) containing the eigenvectors of the fabric.
        
    Returns:
    --------
    None
    """

    # Create a unit sphere and transform it to an ellipsoid
    Sphere = pv.Sphere(radius=ROI.shape[0]/2, theta_resolution=50, phi_resolution=50)

    # Scale the sphere by the square roots of the eigenvalues
    ScaleMatrix = np.diag(np.sqrt(eValues))
    TransformMatrix = np.matmul(eVectors, ScaleMatrix)

    # Transform the sphere points to ellipsoid points
    Points = np.matmul(Sphere.points, TransformMatrix.T)

    # Center the ellipsoid at the structure's midpoint
    Offset = np.array(ROI.shape) / 2
    EllispoidPoints = Points + Offset
    Ellispoid = pv.PolyData(EllispoidPoints, Sphere.faces)

    # Calculate the radius for each ellipsoid point to color by radius
    Radii = np.linalg.norm(Ellispoid.points - Offset, axis=1)
    Radii = (Radii - min(Radii)) / (max(Radii) - min(Radii))
    Radii = Radii * (max(eValues) - min(eValues)) + min(eValues)
    Ellispoid['MIL'] = Radii

    # Plotting
    sargs = dict(font_family='times', 
                    width=0.05,
                    height=0.75,
                    vertical=True,
                    position_x=0.9,
                    position_y=0.125,
                    title_font_size=30,
                    label_font_size=20
                    )
    
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(ROI, cmap='bone', show_scalar_bar=False, opacity=0.5)
    pl.add_mesh(Ellispoid, scalars='MIL', cmap='jet', scalar_bar_args=sargs)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = -60
    pl.camera.zoom(1.0)
    pl.add_bounding_box(color=(0,0,0), line_width=1)
    pl.screenshot(FileName)

    return


#%% Main

def Main():

    # List data
    DataPath = Path(__file__).parents[1] / '00_Data'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])

    # Select a file
    File = Files[0]
    Image = sitk.ReadImage(str(File))
    Array = sitk.GetArrayFromImage(Image).T

    Shape = np.round(np.array(Array.shape) * 0.0065,1)

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
    Time.Process(1,'\nPlot scan')
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
    pl.screenshot(Path(__file__).parent / 'Plots/Scan.png', return_img=False)
    Time.Process(0,'Plot scan')

    # Plot fabric
    Time.Process(1,'\nPlot fabric')
    FabricFile = Path(__file__).parent / 'Results' / (File.name[:-4] + '.fab')
    eValues, eVectors, BVTV = Read.Fabric(FabricFile)
    FileName = Path(__file__).parent / 'Plots/Fabric.png'
    PlotFabricROI(Array, eValues, eVectors, FileName)
    Time.Process(0,'Plot fabric')

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
