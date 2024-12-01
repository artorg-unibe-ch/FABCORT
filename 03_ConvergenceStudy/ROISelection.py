#%% !/usr/bin/env python3

Description = """
Read ISQ files and plot them in 3D using pyvista
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '15-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
import numpy as np
from pathlib import Path
import SimpleITK as sitk
from scipy.spatial import Voronoi

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

#%% Functions

def GetFabric(FileName):

    Text = open(FileName,'r').readlines()
    BVTV = float(Text[12].split('=')[1])
    eValues = np.array(Text[18].split(':')[1].split(),float)
    eVectors = np.zeros((3,3))
    for i in range(3):
        eVectors[i] = Text[19+i].split(':')[1].split()

    return eValues, eVectors, BVTV

def ROISelection(Shape, NROIs, ROISize):
    """
    Implements Lloyd's algorithm (Voronoi Relaxation) to distribute ROIs evenly in a 3D rectangular array
    using NumPy vectorization.
    
    Parameters:
        Shape (tuple): Dimensions of the 3D array (m, n, o).
        NROIs (int): Number of ROIs to generate.
        num_iterations (int): Number of iterations for Lloyd's algorithm.

    Returns:
        np.ndarray: Final ROI positions [(x, y, z), ...].
    """
    # Initialize random points within the bounds of the array
    Xmin, Xmax = ROISize, Shape[0]-ROISize
    Ymin, Ymax = ROISize, Shape[1]-ROISize
    Zmin, Zmax = ROISize, Shape[2]-ROISize
    ROIPoints = np.random.uniform([Xmin, Ymin, Zmin],[Xmax, Ymax, Zmax],(NROIs,3))

    # Add ghost points along boundaries to constrain regions
    Boundaries = np.array([[ Xmin, Ymin, Zmin],
                           [ Xmax, Ymin, Zmin],
                           [ Xmin, Ymax, Zmin],
                           [ Xmin, Ymin, Zmax],
                           [ Xmax, Ymax, Zmin],
                           [ Xmax, Ymin, Zmax],
                           [ Xmin, Ymax, Zmax],
                           [ Xmax, Ymax, Zmax]])

    Points = np.vstack([Boundaries, ROIPoints])
    for Iteration in range(100):

        # Compute the Voronoi diagram
        VDiagram = Voronoi(Points)
        
        # Initialize an array to store the new points
        NewPoints = np.copy(Points)

        # Iterate over Voronoi regions
        for iPoint, PointRegion in enumerate(VDiagram.point_region):

            # Select region
            Region = VDiagram.regions[PointRegion]
                        
            # Skip open regions or degenerate regions
            if -1 in Region or len(Region) == 0:
                continue
            
            # Extract vertices of the region
            Vertices = VDiagram.vertices[Region]
            
            # Clip vertices to stay within the bounds
            Vertices = np.clip(Vertices, [0, 0, 0], Shape)
            
            # Compute the centroid of the region
            Centroid = Vertices.mean(axis=0)
            
            # Add the centroid to the new position
            NewPoints[iPoint] = Centroid

        # Compute distance between actual and new points
        Distances = np.linalg.norm(Points - NewPoints, axis=1)
        if max(Distances) < 0.5:
            break
        
        # Update points
        Points = np.vstack([Boundaries, NewPoints[8:]])

    # Print iterations to reach convergence
    print(f'Converged in {Iteration} iterations')

    return np.round(Points[8:]).astype(int)

def Resample(Image, Factor=None, Size=[None], Spacing=[None], Order=1):

    """
    Resample a SimpleITK image by either a given factor, a new size, or
    a new voxel spacing. Order stands for interpolation order e.g.
    Order = 1: Linear interpolation 
    """

    Dimension = Image.GetDimension()
    OriginalSpacing = np.array(Image.GetSpacing())
    OriginalSize = np.array(Image.GetSize())
    PhysicalSize = OriginalSize * OriginalSpacing

    Origin = Image.GetOrigin()
    Direction = Image.GetDirection()

    if Factor:
        NewSize = [round(Size/Factor) for Size in Image.GetSize()] 
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Size[0]:
        NewSize = Size
        NewSpacing = [PSize/Size for Size, PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing[0]:
        NewSpacing = Spacing
        NewSize = [np.floor(Size/Spacing).astype('int') + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]

    NewArray = np.zeros(NewSize[::-1],'int')
    NewImage = sitk.GetImageFromArray(NewArray)
    NewImage.SetOrigin(Origin - OriginalSpacing/2)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    Resampled = sitk.Resample(Image, NewImage, Transform, Order)
    
    return Resampled

#%% Main

def Main():

    # List bone volume fraction of samples
    DataPath = Path(__file__).parents[1] / '01_Fabric/Results'
    Files = [F for F in DataPath.iterdir()]
    BVTVs = []
    for F in Files:
        _, _, BVTV = GetFabric(F)
        BVTVs.append(BVTV)

    # Load sample with lower BV/TV
    ScanPath = Path(__file__).parents[1] / '00_Data'
    File = Files[BVTVs.index(min(BVTVs))].name[:-4]
    Scan = sitk.ReadImage(str(ScanPath / File) + '.mhd')

    # Select ROIs
    Shape = Scan.GetSize()
    Size = 1/0.0065
    Text = '$Sampling,$ROI,\n'
    for NROIs in range(1,11):
        Coords = ROISelection(Shape, NROIs, Size)

        # Resample ROI resolutions
        for i, Coord in enumerate(Coords):
            X1, X2 = int(Coord[0] - Size // 2), int(Coord[0] + Size // 2)
            Y1, Y2 = int(Coord[1] - Size // 2), int(Coord[1] + Size // 2)
            Z1, Z2 = int(Coord[2] - Size // 2), int(Coord[2] + Size // 2)
            ROI = Scan[X1:X2,Y1:Y2,Z1:Z2]
            Resampled = Resample(ROI, Factor=2)
            Resampled.SetOrigin((0,0,0))
            sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_S{NROIs}_N{i+1}.mhd'))
            Text += f'{NROIs},{i+1},\n'

    # Write parameter file
    with open('Parameters.csv','w') as F:
        F.write(Text)

if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Run script
    Main()
