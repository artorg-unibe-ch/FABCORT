#%% !/usr/bin/env python3

Description = """
Read ISQ files and plot them in 3D using pyvista
"""

__author__ = ['Mathieu Simon']
__date_created__ = '11-11-2024'
__date__ = '11-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
from pathlib import Path
import SimpleITK as sitk

#%% Functions

def WriteBash(Sample, Size, Symmetry):

    Shift = 9
    Main = Sample[:-Shift] + '_Main.inp'
    Job = Sample[:-Shift]
    ODBFile = Sample[:-Shift] + '.odb'
    OutFile = Sample[:-Shift] + '.out'

    Text = f"""abaqus interactive job={Job} inp="/home/ms20s284/FABCORT/Homogenization_{Symmetry}/{Main}" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/{ODBFile}"  out="/home/ms20s284/FABCORT/Homogenization_{Symmetry}/{OutFile}"  size="{Size[0]};{Size[1]};{Size[2]}" spec="Stress" 
rm * 
"""

    with open('RunAbaqus.bash','a') as File:
        File.write(Text)

    return

#%% Main

def Main(Arguments):

    # Read Arguments
    if Arguments.AbaqusInp:
        InputISQs = [Arguments.InputISQ]
    else:
        DataPath = Path(__file__).parents[1] / '02_Results/Mesh/'
        AbaqusInps = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('Mesh.inp')])

    for Input in AbaqusInps[10:12]:
        Image = sitk.ReadImage(str(Input)[:-9] + '.mhd')
        Size = [int(Si * Sp) for Si,Sp in zip(Image.GetSize(), Image.GetSpacing())]
        for Symmetry in ['Isotropic']:
            WriteBash(Input.name, Size, Symmetry)



if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--AbaqusInp', help='Abaqus input file', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)
