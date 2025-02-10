#%% !/usr/bin/env python3

Description = """
Set material parameters for single element tests
"""

__author__ = ['Mathieu Simon']
__date_created__ = '13-12-2024'
__date__ = '13-12-2024'
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



#%% Main

def Main(E=10000, Nu=0.3, DA=2.0):

    # Build isotropic stiffness tensor
    Si = Tensor.Isotropic(E, Nu)

    # Compute transverse isotropic constants for given DA
    E1, E2, E3 = E/DA, E/DA, E
    Nu23, Nu31, Nu12 = Nu, Nu, Nu
    Mu23 = E1 / (2*(1+Nu23))
    Mu31 = E2 / (2*(1+Nu31))
    Mu12 = E3 / (2*(1+Nu12))

    # Build transverse isotropic tensor
    St = Tensor.Othotropic(E1, E2, E3, Mu23, Mu31, Mu12, Nu23, Nu31, Nu12)

    # Fabric
    Fabric = Read.Fabric(Path(__file__).parent / 'Fabric.fab')
    Basis = Fabric[1].T
    Sf = Tensor.Othotropic(E1, E2, E3, Mu23, Mu31, Mu12, Nu23, Nu31, Nu12, Basis)
    
    # Plots
    FName = str(Path(__file__).parent / 'Fabric.png')
    Tensor.PlotFabric(Fabric[0], Fabric[1], FileName=FName)

    FName = str(Path(__file__).parent / 'I_Stiffness.png')
    Tensor.PlotTensor(Si/1E3, FileName=FName)

    FName = str(Path(__file__).parent / 'T_Stiffness.png')
    Tensor.PlotTensor(St/1E3, FileName=FName)

    FName = str(Path(__file__).parent / 'F_Stiffness.png')
    Tensor.PlotTensor(Sf/1E3, FileName=FName)

    # Material parameters form Philippe
    E1 = 12.4425
    E2 = 12.4425
    E3 = 23.7708
    Nu23 = 0.245987
    Nu31 = 0.469944
    Nu12 = 0.34
    Mu23 = 6.41714
    Mu31 = 6.41714
    Mu12 = 4.64274

    # Parameters for Abaqus
    print('Material parameters for Abaqus')
    print('E1 =', E1)
    print('E2 =', E2)
    print('E3 =', E3)
    print('Nu12 =', Nu12)
    print('Nu13 =', Nu31 / E3 * E1)
    print('Nu23 =', Nu23)
    print('G12 =', Mu12)
    print('G13 =', Mu31)
    print('G23 =', Mu23)

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

#%%
