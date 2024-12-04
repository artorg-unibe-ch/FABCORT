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

import pickle
import argparse
import numpy as np
from pathlib import Path
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from matplotlib.patches import Rectangle


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

def HRASampling(Boundaries, NROIs, Factor=10):
    """
    Implements High-Resolution Antialiasing-inspired sampling for even ROI distribution.
    
    Parameters:
        Shape (tuple): Dimensions of the 3D array (m, n, o).
        NROIs (int): Number of ROIs to generate.
        Factor (int): Factor for oversampling the space.

    Returns:
        np.ndarray: ROI positions [(x, y, z), ...].
    """
    
    # Step 1: Oversample the space with a dense grid of candidate points
    Candidates = Factor * NROIs
    Min = Boundaries.min(axis=0)
    Max = Boundaries.max(axis=0)
    Candidates = np.linspace(Min, Max, NROIs*Factor)

    # Step 2: Initialize first ROI point
    if NROIs % 2 == 0:
        ROIs = [Min]
    else:
        ROIs = [Boundaries.mean(axis=0) / 2]

    # Step 3: Iteratively select points to maximize minimum distance
    for i in range(NROIs - 1):

        # Compute distances from all candidates to existing ROIs
        Differences = Candidates[:, np.newaxis, :] - np.array(ROIs)[np.newaxis, :, :]
        Distances = np.linalg.norm(Differences, axis=-1)  # Shape: (num_candidates, len(rois))

        # Select the candidate with the maximum minimum distance
        BestCandidate = np.argmax(Distances.min(axis=1))
        ROIs.append(Candidates[BestCandidate])
        
        # Remove the selected candidate to avoid duplicates
        Candidates = np.delete(Candidates, BestCandidate, axis=0)

    return np.array(ROIs)

def LloydsAlgorithm(Shape, NROIs, ROISize):
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

def FibonacciSphere(N):

    """
    Generate a sphere with (almost) equidistant points
    https://arxiv.org/pdf/0912.4540.pdf
    """

    if N > 1:

        i = np.arange(N)
        Phi = np.pi * (np.sqrt(5) - 1)  # golden angle in radians

        Z = 1 - (i / float(N - 1)) * 2  # z goes from 1 to -1
        Radius = np.sqrt(1 - Z*Z)     # radius at z

        Theta = Phi * i                 # golden angle increment

        X = np.cos(Theta) * Radius
        Y = np.sin(Theta) * Radius

        return np.vstack([X,Y,Z]).T
    
    else:
        return np.array([[0,0,0]])

def PlotROIs(Shape, Boundaries, ROIs, ROISize):
    """
    Plot the ROI results in 2D projections (xy, xz, yz) with contours, boundaries, and ROI squares.
    
    Parameters:
    - Shape (tuple): Dimensions of the 3D space, specified as (X, Y, Z).
    - Boundaries (numpy.ndarray): Array of boundary points, shape (8, 3).
    - ROIs (list): List of ROI coordinates, where each is a tuple (x, y, z).
    - ROISize (float): Size of each ROI square.
    """
    # Define projection planes
    projections = ['xy', 'xz', 'yz']
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, proj in enumerate(projections):
        ax = axes[i]
        ax.set_title(f"Projection: {proj}")
        
        # Set limits and aspect ratio
        if proj == 'xy':
            ax.plot([0, Shape[0]], [0, 0], color=(0,0,1), label='Contour', linewidth=2)
            ax.plot([Shape[0], Shape[0]], [0, Shape[1]], color=(0,0,1), linewidth=2)
            ax.plot([Shape[0], 0], [Shape[1], Shape[1]], color=(0,0,1), linewidth=2)
            ax.plot([0, 0], [Shape[1],0], color=(0,0,1), linewidth=2)
            x1 = Boundaries.min(axis=0)[0] - ROISize / 2
            x2 = Boundaries.max(axis=0)[0] + ROISize / 2
            y1 = Boundaries.min(axis=0)[1] - ROISize / 2
            y2 = Boundaries.max(axis=0)[1] + ROISize / 2
            ax.fill_between([x1,x2], [y1,y1], [y2,y2], color=(0,0,0,0.2), label='Boundaries')
            idx1, idx2 = 0, 1
        elif proj == 'xz':
            ax.plot([0, Shape[0]], [0, 0], color=(0,0,1), label='Contour', linewidth=2)
            ax.plot([Shape[0], Shape[0]], [0, Shape[2]], color=(0,0,1), linewidth=2)
            ax.plot([Shape[0], 0], [Shape[2], Shape[2]], color=(0,0,1), linewidth=2)
            ax.plot([0, 0], [Shape[2],0], color=(0,0,1), linewidth=2)
            x1 = Boundaries.min(axis=0)[0] - ROISize / 2
            x2 = Boundaries.max(axis=0)[0] + ROISize / 2
            y1 = Boundaries.min(axis=0)[2] - ROISize / 2
            y2 = Boundaries.max(axis=0)[2] + ROISize / 2
            ax.fill_between([x1,x2], [y1,y1], [y2,y2], color=(0,0,0,0.2), label='Boundaries')
            idx1, idx2 = 0, 2
        elif proj == 'yz':
            ax.plot([0, Shape[1]], [0, 0], color=(0,0,1), label='Contour', linewidth=2)
            ax.plot([Shape[1], Shape[1]], [0, Shape[2]], color=(0,0,1), linewidth=2)
            ax.plot([Shape[1], 0], [Shape[2], Shape[2]], color=(0,0,1), linewidth=2)
            ax.plot([0, 0], [Shape[2],0], color=(0,0,1), linewidth=2)
            x1 = Boundaries.min(axis=0)[1] - ROISize / 2
            x2 = Boundaries.max(axis=0)[1] + ROISize / 2
            y1 = Boundaries.min(axis=0)[2] - ROISize / 2
            y2 = Boundaries.max(axis=0)[2] + ROISize / 2
            ax.fill_between([x1,x2], [y1,y1], [y2,y2], color=(0,0,0,0.2), label='Boundaries')
            idx1, idx2 = 1, 2
        ax.set_aspect('equal')

        # Step 2: Draw the boundaries in red
        x = Boundaries[:, idx1]
        y = Boundaries[:, idx2]

        # Step 3: Plot ROI squares in blue transparent
        for roi in ROIs:
            x_roi = roi[idx1] - ROISize / 2
            y_roi = roi[idx2] - ROISize / 2
            ax.add_patch(Rectangle(
                            (x_roi, y_roi),
                            ROISize, ROISize,
                            edgecolor='red',
                            facecolor='none',
                            linewidth=2.5,
                            label="ROI Square"
                        ))

        # Add legend
        # ax.legend()

    plt.tight_layout()
    plt.savefig(str(Path(__file__).parent / 'ROIs' / f'N{len(ROIs)}.png'))
    plt.close(fig)

def ROISelectionFibonacci(Shape, NROIs, ROISize):
    """
    Generate a specified number of Regions of Interest (ROIs) within a 3D space, ensuring uniform distribution 
    using Fibonacci sphere sampling and adapting the ROIs to the given constraints.

    Parameters:
    - Shape (tuple): Dimensions of the 3D space, specified as (X, Y, Z).
    - NROIs (int): Number of ROIs to generate.
    - ROISize (float): The size of each ROI, which determines the offset from the edges of the space.

    Returns:
    - list: A list of ROI coordinates, where each coordinate is a tuple (x, y, z).
    """

    # Set limits
    Xmin, Xmax = ROISize/4*3, Shape[0]-ROISize/4*3
    Ymin, Ymax = ROISize/4*3, Shape[1]-ROISize/4*3
    Zmin, Zmax = ROISize/4*3, Shape[2]-ROISize/4*3

    # Set boundary points
    Boundaries = np.array([[ Xmin, Ymin, Zmin],
                           [ Xmax, Ymin, Zmin],
                           [ Xmin, Ymax, Zmin],
                           [ Xmin, Ymin, Zmax],
                           [ Xmax, Ymax, Zmin],
                           [ Xmax, Ymin, Zmax],
                           [ Xmin, Ymax, Zmax],
                           [ Xmax, Ymax, Zmax]])
        
    # Compute center and range
    Center = Boundaries.mean(axis=0)
    Range = (Boundaries.max(axis=0) - Boundaries.min(axis=0))/2

    # Generate directions
    Directions = FibonacciSphere(NROIs)

    # Scale each point to the boundaries
    ROIs = []
    for dx, dy, dz in Directions:

        # Determine scaling factor based on boundaries in the direction
        MaxRadius = min(
            Range[0] / abs(dx) if dx != 0 else float('inf'),
            Range[1] / abs(dy) if dy != 0 else float('inf'),
            Range[2] / abs(dz) if dz != 0 else float('inf')
        )
        if MaxRadius == np.inf:
            MaxRadius = 0
        # Scale and center
        xROI = Center[0] + dx * MaxRadius
        yROI = Center[1] + dy * MaxRadius
        zROI = Center[2] + dz * MaxRadius
        ROIs.append((xROI, yROI, zROI))

    PlotROIs(Shape, Boundaries, ROIs, ROISize)
    
    return ROIs

def ROISelection(Shape, NROIs, ROISize):

    # Set limits
    Xmin, Xmax = ROISize/4*3, Shape[0]-ROISize/4*3
    Ymin, Ymax = ROISize/4*3, Shape[1]-ROISize/4*3
    Zmin, Zmax = ROISize/4*3, Shape[2]-ROISize/4*3

    if NROIs == 1:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)[1:-1]
        YPos = np.linspace(Ymin, Ymax, 3)[1:-1]
        XPos = np.linspace(Xmin, Xmax, 3)[1:-1]

    elif NROIs == 2:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)[1:-1:2]
        YPos = np.linspace(Ymin, Ymax, 3)[1:-1]
        XPos = np.linspace(Xmin, Xmax, 3)[1:-1]

        # Reshape and concatenate
        YPos = np.repeat(YPos, 2)
        XPos = np.repeat(XPos, 2)

    elif NROIs == 3:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 3)[1:-1]
        XPos = np.linspace(Xmin, Xmax, 3)[1:-1]

        # Reshape and concatenate
        YPos = np.repeat(YPos, 3)
        XPos = np.repeat(XPos, 3)

    elif NROIs == 4:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)[1:-1:2]
        YPos = np.linspace(Ymin, Ymax, 2)
        XPos = np.linspace(Xmin, Xmax, 2)

        # Repeat to create grid
        ZPos = np.repeat(ZPos, 2)
        YPos = np.tile(YPos, 2)
        XPos = np.concatenate([XPos, XPos[::-1]])

    elif NROIs == 5:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[1], ZPos[2], ZPos[2]])
        YPos = np.array([YPos[0], YPos[2], YPos[1], YPos[0], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[1], XPos[2], XPos[0]])

    elif NROIs == 6:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 2)
        XPos = np.linspace(Xmin, Xmax, 3)[1:-1]

        # Repeat to create grid
        ZPos = np.repeat(ZPos, 2)
        YPos = np.tile(YPos,3)
        XPos = np.tile(XPos,6)

    elif NROIs == 7:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[4], ZPos[4], ZPos[1], ZPos[2], ZPos[3]])
        YPos = np.array([YPos[0], YPos[2], YPos[0], YPos[2], YPos[1], YPos[1], YPos[1]])
        XPos = np.array([XPos[0], XPos[2], XPos[2], XPos[0], XPos[1], XPos[1], XPos[1]])

    elif NROIs == 8:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)[1:-1:2]
        YPos = np.linspace(Ymin, Ymax, 2)
        XPos = np.linspace(Xmin, Xmax, 2)

        # Repeat to create grid
        ZPos = np.repeat(ZPos, 4)
        YPos = np.repeat(np.tile(YPos,2),2)
        XPos = np.tile(XPos,4)

    elif NROIs == 9:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0],
                         ZPos[0], ZPos[1], ZPos[2],
                         ZPos[2], ZPos[2], ZPos[2]])
        YPos = np.array([YPos[0], YPos[0], YPos[2],
                         YPos[2], YPos[1], YPos[0],
                         YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0],
                         XPos[2], XPos[1], XPos[0],
                         XPos[2], XPos[0], XPos[2]])

    elif NROIs == 10:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 4)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0],
                         ZPos[0], ZPos[1], ZPos[2], ZPos[3],
                         ZPos[3], ZPos[3], ZPos[3]])
        YPos = np.array([YPos[0], YPos[0], YPos[2],
                         YPos[2], YPos[1], YPos[1], YPos[0],
                         YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0],
                         XPos[2], XPos[1], XPos[1], XPos[0],
                         XPos[2], XPos[0], XPos[2]])

    elif NROIs == 11:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0], ZPos[0],
                         ZPos[1], ZPos[2], ZPos[3],
                         ZPos[4], ZPos[4], ZPos[4], ZPos[4]])
        YPos = np.array([YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[1], YPos[1], YPos[1],
                         YPos[0], YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[1], XPos[1], XPos[1],
                         XPos[0], XPos[2], XPos[0], XPos[2]])

    elif NROIs == 12:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 2)
        XPos = np.linspace(Xmin, Xmax, 2)

        # Repeat to create grid
        ZPos = np.repeat(ZPos, 4)
        YPos = np.repeat(np.tile(YPos,3),2)
        XPos = np.tile(XPos,6)

    elif NROIs == 13:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 3)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0], ZPos[0],
                         ZPos[1], ZPos[1], ZPos[1], ZPos[1], ZPos[1],
                         ZPos[2], ZPos[2], ZPos[2], ZPos[2]])
        YPos = np.array([YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[0], YPos[0], YPos[2], YPos[2], YPos[1],
                         YPos[0], YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[0], XPos[2], XPos[0], XPos[2], XPos[1],
                         XPos[0], XPos[2], XPos[0], XPos[2]])

    elif NROIs == 14:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0], ZPos[0],
                         ZPos[1],
                         ZPos[2], ZPos[2], ZPos[2], ZPos[2],
                         ZPos[3],
                         ZPos[4], ZPos[4], ZPos[4], ZPos[4]])
        YPos = np.array([YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[1],
                         YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[1],
                         YPos[0], YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[1],
                         XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[1],
                         XPos[0], XPos[2], XPos[0], XPos[2]])

    elif NROIs == 15:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 5)
        YPos = np.linspace(Ymin, Ymax, 3)
        XPos = np.linspace(Xmin, Xmax, 3)

        # Repeat to create grid
        ZPos = np.array([ZPos[0], ZPos[0], ZPos[0], ZPos[0],
                         ZPos[1],
                         ZPos[2], ZPos[2], ZPos[2], ZPos[2], ZPos[2],
                         ZPos[3],
                         ZPos[4], ZPos[4], ZPos[4], ZPos[4]])
        YPos = np.array([YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[1],
                         YPos[1], YPos[0], YPos[0], YPos[2], YPos[2],
                         YPos[1],
                         YPos[0], YPos[0], YPos[2], YPos[2]])
        XPos = np.array([XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[1],
                         XPos[1], XPos[0], XPos[2], XPos[0], XPos[2],
                         XPos[1],
                         XPos[0], XPos[2], XPos[0], XPos[2]])

    elif NROIs == 16:
        # Compute positions
        ZPos = np.linspace(Zmin, Zmax, 4)
        YPos = np.linspace(Ymin, Ymax, 2)
        XPos = np.linspace(Xmin, Xmax, 2)

        # Repeat to create grid
        ZPos = np.repeat(ZPos, 4)
        YPos = np.repeat(np.tile(YPos,4),2)
        XPos = np.tile(XPos,8)

    # Reshape and concatenate
    ROIs = np.vstack([XPos, YPos, ZPos]).T

    Boundaries = np.array([[ Xmin, Ymin, Zmin],
                           [ Xmax, Ymin, Zmin],
                           [ Xmin, Ymax, Zmin],
                           [ Xmin, Ymin, Zmax],
                           [ Xmax, Ymax, Zmin],
                           [ Xmax, Ymin, Zmax],
                           [ Xmin, Ymax, Zmax],
                           [ Xmax, Ymax, Zmax]])
        
    PlotROIs(Shape, Boundaries, ROIs, ROISize)

    return ROIs

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
    Text = '$ROI,\n'
    ROICoords = []
    ROINumbers = []
    for NROIs in range(1,17):

        # Select ROIs positions
        Coords = ROISelection(Shape, NROIs, Size)

        # Iterate over every positions
        ROINumber = []
        for Coord in Coords:
            
            # If already in list, skip
            if any([(Coord == C).all() for C in ROICoords]):
                for i, C in enumerate(ROICoords):
                    if (Coord == C).all():
                        ROINumber.append(i)
                continue

            # If not, store coords, resample and write ROI
            else:
                ROICoords.append(Coord)
                ROINumber.append(len(ROICoords))
                # X1, X2 = int(Coord[0] - Size // 2), int(Coord[0] + Size // 2)
                # Y1, Y2 = int(Coord[1] - Size // 2), int(Coord[1] + Size // 2)
                # Z1, Z2 = int(Coord[2] - Size // 2), int(Coord[2] + Size // 2)
                # ROI = Scan[X1:X2,Y1:Y2,Z1:Z2]
                # Resampled = Resample(ROI, Factor=2)
                # Resampled.SetOrigin((0,0,0))
                # sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_{len(ROICoords):03d}.mhd'))
                # Text += f'{len(ROICoords):03d},\n'
        ROINumbers.append(ROINumber)

    # # Write parameter file
    # with open(Path(__file__).parent / 'Parameters.csv','w') as F:
    #     F.write(Text)

    # Save map
    with open(Path(__file__).parent / 'ROIMap.', 'wb') as F:
        pickle.dump(ROINumbers, F)

if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Run script
    Main()
