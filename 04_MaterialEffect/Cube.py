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
import pyvista as pv
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Read, Tensor

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
    actors = pl.add_volume(ROI,
                           cmap='bone',
                           show_scalar_bar=False, opacity=[0.005,0])
    actors.prop.interpolation_type = 'linear'
    pl.add_mesh(Ellispoid, scalars='MIL', cmap='jet', scalar_bar_args=sargs)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = 30
    pl.camera.zoom(1.0)
    pl.add_bounding_box(color=(0,0,0), line_width=1)
    # pl.add_axes(viewport=(0,0,0.25,0.25))
    pl.screenshot(FileName)
    # pl.show()

    return

def PlotStiffnessROI(ROI:np.array, StiffnessTensor:np.array, FileName:Path) -> None:

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
    Sphere = pv.Sphere(radius=1, theta_resolution=50, phi_resolution=50)

    # Compute elongation and bulk modulus
    I = np.eye(3)
    ElongationModulus = np.zeros(Sphere.points.shape)
    BulkModulus = np.zeros(len(Sphere.points))
    for p, Point in enumerate(Sphere.points):
        N = Tensor.DyadicProduct(Point, Point)
        SN = Tensor.Transform(StiffnessTensor, N)
        ElongationModulus[p] = Tensor.FrobeniusProduct(N, SN) * Point
        BulkModulus[p] = Tensor.FrobeniusProduct(I, SN)

    # # Scale the sphere by the square roots of the eigenvalues
    # Scale = ROI.shape[0]/2 / max(np.linalg.norm(ElongationModulus, axis=1))
    # Points = ElongationModulus * Scale

    # # Center the ellipsoid at the structure's midpoint
    # Offset = np.array(ROI.shape) / 2
    # EllispoidPoints = Points + Offset
    Ellispoid = pv.PolyData(ElongationModulus, Sphere.faces)
    Ellispoid['Bulk Modulus'] = BulkModulus

    # Plotting
    SArgs = dict(font_family='times', 
                 width=0.05,
                 height=0.75,
                 vertical=True,
                 position_x=0.85,
                 position_y=0.125,
                 title_font_size=30,
                 label_font_size=20
                 )
    
    BArgs = dict(font_family='times', 
                 font_size=30,
                 location='back',
                 n_xlabels=3,
                 n_ylabels=3,
                 n_zlabels=3,
                 all_edges=True,
                 fmt='%i',
                 xtitle='',
                 ytitle='',
                 ztitle='',
                 use_3d_text=False
                 )
    
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(Ellispoid, scalars='Bulk Modulus', cmap='jet', scalar_bar_args=SArgs)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = 30
    pl.camera.zoom(0.8)
    pl.add_axes(viewport=(0,0,0.25,0.25),
                label_size=(0.065, 0.065),
                xlabel='e1',
                ylabel='e2',
                zlabel='e3')
    pl.show_bounds(**BArgs)
    pl.screenshot(FileName, return_img=False, scale=2)
    # pl.show()

    return

#%% Main

def Main():

    DataPath = Path(__file__).parent
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])

    # Get fabric info
    Fabric = Read.Fabric(DataPath / 'Fabric.fab')
    eValues = Fabric[0]
    eVectors = Fabric[1]

    # Get homogenization stress results
    Isotropic = open(DataPath / 'Cube_Isotropic.out', 'r').readlines()
    Transverse = open(DataPath / 'Cube_Transverse.out', 'r').readlines()
    Fabric = open(DataPath / 'Cube_Fabric.out', 'r').readlines()

    IsoStress = np.zeros((6,6))
    TransStress = np.zeros((6,6))
    FabStress = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            IsoStress[i,j] = float(Isotropic[i+4].split()[j+1])
            TransStress[i,j] = float(Transverse[i+4].split()[j+1])
            FabStress[i,j] = float(Fabric[i+4].split()[j+1])

    IsoStiffness = np.zeros((6,6))
    TransStiffness = np.zeros((6,6))
    FabStiffness = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            IsoStiffness[i,j] = IsoStress[i,j] / Strain[i]
            TransStiffness[i,j] = TransStress[i,j] / Strain[i]
            FabStiffness[i,j] = FabStress[i,j] / Strain[i]

    # Symetrize matrices
    IsoStiffness = 1/2 * (IsoStiffness + IsoStiffness.T)
    TransStiffness = 1/2 * (TransStiffness + TransStiffness.T)
    FabStiffness = 1/2 * (FabStiffness + FabStiffness.T)

    # Compute material parameters
    Material = open(DataPath / 'Material_Transverse.inp', 'r').readlines()
    Material = Material[-4][:-1] + ', ' + Material[-3][:-1]
    Material = [float(M) for M in Material.split(', ')]
    E1, E2, E3 = Material[0], Material[1], Material[2]
    Nu23, Nu31, Nu12 = Material[3], Material[4], Material[5]
    Mu23, Mu31, Mu12 = Material[8], Material[7], Material[6]
    S = np.array([[1/E1, -Nu23/E1, -Nu31/E1, 0, 0, 0],
                  [-Nu23/E1, 1/E2, -Nu12/E2, 0, 0, 0],
                  [-Nu31/E1, -Nu12/E2, 1/E3, 0, 0, 0],
                  [0, 0, 0, 1/Mu23, 0, 0],
                  [0, 0, 0, 0, 1/Mu31, 0],
                  [0, 0, 0, 0, 0, 1/Mu12]])
    C = np.linalg.inv(S)

    # Build 3x3x3x3 tensor
    IsoStiffness = Tensor.IsoMorphism66_3333(IsoStiffness)
    TransStiffness = Tensor.IsoMorphism66_3333(TransStiffness)
    FabStiffness = Tensor.IsoMorphism66_3333(FabStiffness)

    # Plot
    FName = str(Path(__file__).parent / 'FI_Stiffness.png')
    Tensor.PlotTensor(IsoStiffness/1E3, FileName=FName)

    FName = str(Path(__file__).parent / 'FT_Stiffness.png')
    Tensor.PlotTensor(TransStiffness/1E3, FileName=FName)


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
