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

    DataPath = Path(__file__).parents[1] / 'Tests'
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

    # Transform fabric stiffness into fabric coordinate system
    Stiffness = 1/2 * (FabStiffness + FabStiffness.T)

    # Write tensor into mandel notation
    Mandel = Tensor.Engineering2MandelNotation(Stiffness)

    # Step 3: Transform tensor into fabric coordinate system
    I = np.eye(3)
    Q = np.array(eVectors)
    Transformed = Tensor.Transform(Mandel, I, Q)

    # Project onto orthotropy
    Orthotropic = np.zeros(Transformed.shape)
    for i in range(Orthotropic.shape[0]):
        for j in range(Orthotropic.shape[1]):
            if i < 3 and j < 3:
                Orthotropic[i, j] = Transformed[i, j]
            elif i == j:
                Orthotropic[i, j] = Transformed[i, j]

    # Get tensor back to engineering notation
    Stiffness = Tensor.Mandel2EngineeringNotation(Orthotropic)

    # Print results
    print('\nStiffness matrix of isotropic material')
    print(IsoStiffness)

    print('\nStiffness matrix of transverse isotropic material')
    print(TransStiffness)

    print('\nStiffness matrix of transverse isotropic material in fabric coordinate system')
    print(Stiffness)

    E = 0.001
    Nu = 0.3
    np.matmul(IsoStiffness, np.array([E, -E*Nu, -E*Nu, 0, 0, 0,]))

    E = 0.001
    np.matmul(TransStiffness, np.array([E, -E*0.2, -E*0.3, 0, 0, 0,]))


if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--Sample', help='Sample main file name', type=str)
    Parser.add_argument('--OutputPath', help='Output path for the ROI and png image of the plot', default=Path(__file__).parents[1] / '02_Results/Scans')
    Parser.add_argument('--NROIs', help='Number of region of interests to extract', type=int, default=3)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)
