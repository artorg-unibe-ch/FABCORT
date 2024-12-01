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

    # Select first ROI of about 1 mm3 around the center
    Center = np.array(Scan.GetSize()) / 2
    Size = round(1 / 0.0065)
    X,Y,Z = np.round(Center - Size / 2)
    ROI = Scan[int(X):int(X+Size), int(Y):int(Y+Size), int(Z):int(Z+Size)]

    # Resample ROI for different resolutions
    for F in [1,2,4]:
        Resampled = Resample(ROI, Factor=F)
        sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_1_F{F}.mhd'))

    # Select second ROI in top third
    Center = np.array(Scan.GetSize()) / 3
    Size = round(1 / 0.0065)
    X,Y,Z = np.round(Center - Size / 2)
    ROI = Scan[int(X):int(X+Size), int(Y):int(Y+Size), int(Z):int(Z+Size)]

    # Resample ROI for different resolutions
    for F in [1,2,4]:
        Resampled = Resample(ROI, Factor=F)
        sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_2_F{F}.mhd'))

    # Select third ROI in last third
    Center = np.array(Scan.GetSize()) / 3 * 2
    Size = round(1 / 0.0065)
    X,Y,Z = np.round(Center - Size / 2)
    ROI = Scan[int(X):int(X+Size), int(Y):int(Y+Size), int(Z):int(Z+Size)]

    # Resample ROI for different resolutions
    for F in [1,2,4]:
        Resampled = Resample(ROI, Factor=F)
        sitk.WriteImage(Resampled,str(Path(__file__).parent / 'ROIs' / f'ROI_3_F{F}.mhd'))


if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Run script
    Main()
