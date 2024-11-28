#%% !/usr/bin/env python3

Description = """
To write
"""

__author__ = ['Mathieu Simon']
__date_created__ = '21-11-2024'
__date__ = '21-11-2024'
__license__ = 'MIT'
__version__ = '1.0'


#%% Imports

import time
import argparse
import numpy as np
import pandas as pd
import pyvista as pv
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# %% Time class
class Time():

    def __init__(self):
        self.Width = 15
        self.Length = 16
        self.Text = 'Process'
        self.Tic = time.time()
        pass
    
    def Set(self, Tic=None):
        
        if Tic == None:
            self.Tic = time.time()
        else:
            self.Tic = Tic

    def Print(self, Tic=None,  Toc=None):

        """
        Print elapsed time in seconds to time in HH:MM:SS format
        :param Tic: Actual time at the beginning of the process
        :param Toc: Actual time at the end of the process
        """

        if Tic == None:
            Tic = self.Tic
            
        if Toc == None:
            Toc = time.time()


        Delta = Toc - Tic

        Hours = np.floor(Delta / 60 / 60)
        Minutes = np.floor(Delta / 60) - 60 * Hours
        Seconds = Delta - 60 * Minutes - 60 * 60 * Hours

        print('\nProcess executed in %02i:%02i:%02i (HH:MM:SS)' % (Hours, Minutes, Seconds))

        return

    def Update(self, Progress, Text=''):

        Percent = int(round(Progress * 100))
        Np = self.Width * Percent // 100
        Nb = self.Width - Np

        if len(Text) == 0:
            Text = self.Text
        else:
            self.Text = Text

        Ns = self.Length - len(Text)
        if Ns >= 0:
            Text += Ns*' '
        else:
            Text = Text[:self.Length]
        
        Line = '\r' + Text + ' [' + Np*'=' + Nb*' ' + ']' + f' {Percent:.0f}%'
        print(Line, sep='', end='', flush=True)

    def Process(self, StartStop:bool, Text=''):

        if len(Text) == 0:
            Text = self.Text
        else:
            self.Text = Text

        if StartStop*1 == 1:
            self.Tic = time.time()
            self.Update(0, Text)

        elif StartStop*1 == 0:
            self.Update(1, Text)
            self.Print()

Time = Time()

#%% Functions

def GetFabric(FileName):

    Text = open(FileName,'r').readlines()
    BVTV = float(Text[12].split('=')[1])
    eValues = np.array(Text[18].split(':')[1].split(),float)
    eVectors = np.zeros((3,3))
    for i in range(3):
        eVectors[i] = Text[19+i].split(':')[1].split()

    return eValues, eVectors, BVTV

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

def ComputeFabric(S, BVTV=1):

    """
    Compute the parameters of the Zysset-Curnier model from a stiffness matrix.

    This function extracts the material parameters of the Zysset-Curnier model based on 
    a given stiffness matrix `S` and a specified bone volume fraction (BV/TV). It uses 
    the diagonal terms of the stiffness matrix to compute the fabric eigenvalues and 
    subsequently the anisotropic material constants.

    Parameters
    ----------
    S : numpy.ndarray
        A 6x6 symmetric stiffness matrix representing the elastic properties of trabecular bone.
        The matrix structure corresponds to the Zysset-Curnier model.
    BVTV : float, optional
        Bone volume fraction (BV/TV), representing the ratio of bone volume to total volume.
        Default is 0.95.

    Returns
    -------
    Lambda0 : float
        The bulk modulus parameter scaled by BV/TV, representing the isotropic elasticity.
    Mu0 : float
        The shear modulus parameter scaled by BV/TV, representing the anisotropic shear response.
    k : float
        The exponent describing the relationship between stiffness and BV/TV.
    l : float
        The anisotropy exponent describing the scaling of stiffness with principal directions.
    m1, m2, m3 : float
        Fabric eigenvalues (normalized), representing the bone's anisotropy.

    Notes
    -----
    The function implements the relationships described in the Zysset-Curnier model:
        - Fabric eigenvalues (`m1`, `m2`, `m3`) are computed from the ratios of the diagonal 
          components of the stiffness matrix.
        - Material constants scaled by BV/TV (`Lambda0_BVTV`, `Mu0_BVTV`) are calculated 
          from the diagonal terms and shear components.
        - The BV/TV scaling exponent `k` is determined based on the relationship between 
          bulk modulus and BV/TV.
        - The Lamé constants (`Lambda0`, `Mu0`) are then computed by normalizing with respect 
          to BV/TV scaling.

    References
    ----------
    Zysset, P. K. "A review of morphology–elasticity relationships in human trabecular bone: 
    theories and experiments." Journal of Biomechanics.

    Example
    -------
    Given a stiffness matrix `S`, the function can be called as follows:

    >>> S = np.array([[1.2, 0.3, 0.3, 0, 0, 0],
                      [0.3, 1.1, 0.3, 0, 0, 0],
                      [0.3, 0.3, 1.0, 0, 0, 0],
                      [0, 0, 0, 0.8, 0, 0],
                      [0, 0, 0, 0, 0.8, 0],
                      [0, 0, 0, 0, 0, 0.8]])
    >>> Lambda0, Mu0, k, l, m1, m2, m3 = ComputeFabric(S, BVTV=0.95)
    >>> print(f"Lambda0={Lambda0}, Mu0={Mu0}, k={k}, l={l}, m1={m1}, m2={m2}, m3={m3}")
    """

    # Extract diagonal terms
    S11, S22, S33 = S[0, 0], S[1, 1], S[2, 2]
    S44, S55, S66 = S[3, 3], S[4, 4], S[5, 5]

    # Ratios to compute fabric eigenvalues (m1, m2, m3)
    m1_m2_l_ratio = (S11 / S22)**0.5
    m1_m3_l_ratio = (S11 / S33)**0.5

    # Normalize fabric eigenvalues
    m1_l = 1 / (1 + 1/m1_m2_l_ratio + 1/m1_m3_l_ratio)
    m2_l = m1_l / m1_m2_l_ratio
    m3_l = m1_l / m1_m3_l_ratio

    # Compute LambdaS_BVTV, Mu0_BVTV, La
    Lambda0P_BVTV = S11 / (m1_l**2)
    Mu0_BVTV = S44 / (2 * m2_l * m3_l)




    # Compute l (anisotropy exponent)
    l = 0.5 * np.log(S22 / S33) / np.log(m2_l / m3_l)



    # Compute Lambda0_BVTV, Mu0_BVTV
    Lambda0P_BVTV = S11 / (m1_l**2)
    Mu0_BVTV = S44 / (2 * m2_l * m3_l)
    Lambda0_BVTV = Lambda0P_BVTV - 2 * Mu0_BVTV

    # Compute BVTV^k
    BVTV_k = (Lambda0_BVTV + 2 * Mu0_BVTV) / Lambda0P_BVTV

    # Compute k (BVTV scaling exponent)
    k = np.log(BVTV_k) / np.log(BVTV)

    # Compute Lambda0 and Mu0 (Lamé constants)
    Lambda0 = Lambda0_BVTV / BVTV_k
    Mu0 = Mu0_BVTV / BVTV_k

    return Lambda0, Mu0, k, l, m1, m2, m3

def Stiffness(Lambda0, Lambda0p, Mu0, m1, l=1):
    
    m3 = 3 - 2*m1

    S = np.array([[(Lambda0+2*Mu0)*m1**(2*l), Lambda0p*m1**l*m1**l, Lambda0p*m1**l*m3**l, 0, 0, 0],
                  [Lambda0p*m1**l*m1**l, (Lambda0+2*Mu0)*m1**(2*l), Lambda0p*m1**l*m3**l, 0, 0, 0],
                  [Lambda0p*m3**l*m1**l, Lambda0p*m3**l*m1**l, (Lambda0+2*Mu0)*m3**(2*l), 0, 0, 0],
                  [0, 0, 0, 2*Mu0*m1**l*m3**l, 0, 0],
                  [0, 0, 0, 0, 2*Mu0*m3**l*m1**l, 0],
                  [0, 0, 0, 0, 0, 2*Mu0*m1**l*m1**l]])

    return S

def Cost(params, C):
    Lambda0, Lambda0p, Mu0, m1 = params
    S = Stiffness(Lambda0, Lambda0p, Mu0, m1)
    return np.linalg.norm(S - C)

def FitFabric(S):
    
    # Initial guess
    Guess = [10, 10, 5, 1]

    # Bounds for parameters (optional, adjust as necessary)
    Bounds = [(0, None), (0, None), (0, None), (0, 3)]

    # Minimize the objective function
    Results = minimize(Cost, Guess, args=(S), bounds=Bounds)
    Lambda0, Lambda0p, Mu0, m1 = Results.x

    return Lambda0, Lambda0p, Mu0, m1


#%% Main

def Main(Arguments):

    DataPath = Path(__file__).parents[1] / '00_Data'
    Data = pd.read_excel(DataPath / 'Cij.xls')

    # Build stiffness tensor
    i = 0
    Lambda0s, Lambda0ps, Mu0s, m1s = [], [], [], []
    for i in Data.index:
        C11 = Data.loc[i,'RUS:C11']
        C33 = Data.loc[i,'RUS:C33']
        C44 = Data.loc[i,'RUS:C44']
        C66 = Data.loc[i,'RUS:C66']
        C13 = Data.loc[i,'RUS:C13']
        C12 = C11 - 2*C66

        C = np.array([[C11, C12, C13,   0,   0,   0],
                      [C12, C11, C13,   0,   0,   0],
                      [C13, C13, C33,   0,   0,   0],
                      [  0,   0,   0, C44,   0,   0],
                      [  0,   0,   0,   0, C44,   0],
                      [  0,   0,   0,   0,   0, C66]])
    
    """ Constants according to
    Mirzaali et al.
    Mechanical properties of cortical bone and their relationships with age,
    gender, composition and microindentation properties in the elderly
    Bone, 2016
    *ELASTIC, TYPE=ENGINEERING CONSTANTS
    11200, 11200, 19900, 0.3, 0.4, 0.225, 4300, 5700
    5700
    """
    

    # Get structural fabric
    Label = Data.loc[i,'Measurement label.1']
    Label = f'{Label[5:8]}_{Label[:4]}_{Label[9:11]}'
    Measurement = [F.name[:-4] for F in DataPath.iterdir() if Label in F.name and F.name.endswith('.mhd')]
    Name = Path(__file__).parents[1] / '02_Results/Fabric' / (Measurement[0] + '.fab')
    eValues, eVectors, BVTV = GetFabric(Name)




    #
    FileName = str(Name)[:-4] + '.png'
    PlotFabricROI(Array, eValues, eVectors, FileName)




    return



if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--OutputPath', help='Output path for the plots', default=Path(__file__).parents[1] / '02_Results')

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)

#%%
