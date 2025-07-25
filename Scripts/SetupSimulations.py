#%% !/usr/bin/env python3

Description = """
Prepare Abaqus input and batch files for homogenization simulations of bone ROIs.

This script automates the modification of Abaqus `.inp` files and the generation of
corresponding bash scripts to run simulations on multiple cubic regions of interest (ROIs),
typically extracted from 3D bone imaging data.

Functionality:
1. Iterates over all ROI folders in the specified directory.
2. For each ROI and material type (e.g., Isotropic, Transverse), it:
   - Computes the physical ROI size from image metadata.
   - Modify a medtool generated Abaqus input file defining mesh, material, boundary conditions,
     and six simulation steps (3 tension, 3 pure shear).
   - Appends a command to a batch-specific bash script that runs the simulation and extracts results.
3. Batch scripts are split evenly across a fixed number of jobs (e.g., 6), making them easy to submit to HPC systems.

Inputs:
    DataPath (str): Path to the folder containing ROI directories with Abaqus mesh files and `.mhd` images.

Outputs:
    - Abaqus input files: Written into each ROI directory.
    - Bash scripts: Written into `Scripts/Medtool/RunAbaqus_<n>.bash` for batch processing.

Requirements:
    - ROI naming format: ROI_<x><y><z>_<Type>
    - Image metadata in accompanying `.mhd` files
    - Pre-existing mesh (`Mesh.inp`) and boundary condition (`KUBCs.inp`) files

Example:
    python GenerateAbaqusInputs.py ../Results/Cortical/Homogenisation

"""

__author__ = ['Mathieu Simon']
__date_created__ = '10-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
from pathlib import Path
import SimpleITK as sitk

#%% Functions

def WriteMain(Sample, FileName, Size, Type='Isotropic'):

    L1, L2, L3 = Size
    Shift = 8
    Mesh = FileName.name[:Shift] + 'Mesh.inp'
    BCs = FileName.name[:Shift] + 'KUBCs.inp'
    Main = FileName.name[:Shift] + f'Main_{Type}.inp'

    Text = f"""**********************************************
**
**         Main Abaqus input File
**
**       Homogenization of cubic ROI
**     
**    Mathieu Simon, ARTORG Center, 2024 
**
**********************************************
** Paramter Definition for Steps
**  (Unit Cell Dimensions l1,l2,l3=h)
**********************************************
*PARAMETER
u1  = {L1/1000}
u2  = {L2/1000}
u3  = {L3/1000}
**
** Node, Element, and Material Definitons 
**********************************************
*INCLUDE, INPUT=/home/ms20s284/FABCORT/Homogenization/{Sample}/{Mesh}
*INCLUDE, INPUT=/home/ms20s284/FABCORT/Homogenization/Material_{Type}.inp
**
** Interactions (*Equation and *Nset)
**********************************************
*INCLUDE, INPUT=/home/ms20s284/FABCORT/Homogenization/{Sample}/{BCs}
**
** Steps Definitions
***************** Tensile 1 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 1, <u1>
SEB, 2, 2, 0
SEB, 3, 3, 0
NWB, 1, 3, 0
NEB, 1, 1, <u1>
NEB, 2, 2, 0
NEB, 3, 3, 0
SWT, 1, 3, 0
SET, 1, 1, <u1>
SET, 2, 2, 0
SET, 3, 3, 0
NWT, 1, 3, 0
NET, 1, 1, <u1>
NET, 2, 2, 0
NET, 3, 3, 0
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
***************** Tensile 2 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 3, 0
NWB, 1, 1, 0
NWB, 2, 2, <u2>
NWB, 3, 3, 0
NEB, 1, 1, 0
NEB, 2, 2, <u2>
NEB, 3, 3, 0
SWT, 1, 3, 0
SET, 1, 3, 0
NWT, 1, 1, 0
NWT, 2, 2, <u2>
NWT, 3, 3, 0
NET, 1, 1, 0
NET, 2, 2, <u2>
NET, 3, 3, 0
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
***************** Tensile 3 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 3, 0
NWB, 1, 3, 0
NEB, 1, 3, 0
SWT, 1, 1, 0
SWT, 2, 2, 0
SWT, 3, 3, <u3>
SET, 1, 1, 0
SET, 2, 2, 0
SET, 3, 3, <u3>
NWT, 1, 1, 0
NWT, 2, 2, 0
NWT, 3, 3, <u3>
NET, 1, 1, 0
NET, 2, 2, 0
NET, 3, 3, <u3>
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
****************** Shear 23 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 3, 0
NWB, 1, 1, 0
NWB, 2, 2, 0
NWB, 3, 3, <u2>
NEB, 1, 1, 0
NEB, 2, 2, 0
NEB, 3, 3, <u2>
SWT, 1, 1, 0
SWT, 2, 2, <u3>
SWT, 3, 3, 0
SET, 1, 1, 0
SET, 2, 2, <u3>
SET, 3, 3, 0
NWT, 1, 1, 0
NWT, 2, 2, <u3>
NWT, 3, 3, <u2>
NET, 1, 1, 0
NET, 2, 2, <u3>
NET, 3, 3, <u2>
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
****************** Shear 13 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 1, 0
SEB, 2, 2, 0
SEB, 3, 3, <u1>
NWB, 1, 3, 0
NEB, 1, 1, 0
NEB, 2, 2, 0
NEB, 3, 3, <u1>
SWT, 1, 1, <u3>
SWT, 2, 2, 0
SWT, 3, 3, 0
SET, 1, 1, <u3>
SET, 2, 2, 0
SET, 3, 3, <u1>
NWT, 1, 1, <u3>
NWT, 2, 2, 0
NWT, 3, 3, 0
NET, 1, 1, <u3>
NET, 2, 2, 0
NET, 3, 3, <u1>
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
****************** Shear 21 ******************
*STEP
*STATIC
*BOUNDARY, TYPE=DISPLACEMENT
SWB, 1, 3, 0
SEB, 1, 1, 0
SEB, 2, 2, <u1>
SEB, 3, 3, 0
NWB, 1, 1, <u2>
NWB, 2, 2, 0
NWB, 3, 3, 0
NEB, 1, 1, <u2>
NEB, 2, 2, <u1>
NEB, 3, 3, 0
SWT, 1, 3, 0
SET, 1, 1, 0
SET, 2, 2, <u1>
SET, 3, 3, 0
NWT, 1, 1, <u2>
NWT, 2, 2, 0
NWT, 3, 3, 0
NET, 1, 1, <u2>
NET, 2, 2, <u1>
NET, 3, 3, 0
** Element Output 
*OUTPUT, FIELD
*ELEMENT OUTPUT
IVOL, S, E
*END STEP
"""


    FilePath = FileName.parent / Main
    with open(FilePath,'w') as File:
        File.write(Text)


    # Check that neither step or material is written in the mesh file
    with open(FileName.parent / Mesh, 'r') as File:
        Text = File.read()
    StepStart = Text.find('*STEP')
    MaterialStart = Text.find('** Material')

    # If yes, remove it
    if StepStart > 0 or MaterialStart > 0:
        if StepStart < 0:
            Stop = MaterialStart
        elif MaterialStart < 0:
            Stop = StepStart
        else:
            Stop = min(StepStart, MaterialStart)
        with open(FileName.parent / Mesh, 'w') as File:
            File.write(Text[:Stop])

    return

def WriteBash(BatchNumber, Sample, FileName, Size, Type='Isotropic'):

    Shift = 8
    Main = FileName[:Shift] + f'Main_{Type}.inp'
    Job = FileName[:Shift] + f'{Type}'
    ODBFile = FileName[:Shift] + f'{Type}.odb'
    OutFile = FileName[:Shift] + f'{Type}.out'

    Text = f"""abaqus interactive job={Job} inp="/home/ms20s284/FABCORT/Homogenization/{Sample}/{Main}" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations{BatchNumber}/{ODBFile}"  out="/home/ms20s284/FABCORT/Homogenization/{Sample}/{OutFile}"  size="{Size[0]};{Size[1]};{Size[2]}" spec="Stress" 
rm * 
"""

    with open(Path(__file__).parent / f'RunAbaqus_{BatchNumber}.bash','a') as File:
        File.write(Text)

    return

#%% Main

def Main(DataPath):

    # List folders
    DataPath = Path(DataPath)
    Folders = sorted(F for F in DataPath.iterdir() if F.is_dir())

    # Abaqus input files
    AbaqusInps = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                for Type in ['Isotropic', 'Transverse']:
                    AbaqusInps.append(f'ROI_{x+1}{y+1}{z+1}_{Type}')
        
    # Full file path
    AllInps = []
    for F in Folders:
        for I in AbaqusInps:
            AllInps.append(F / I)

    # Write bash script
    NBatch = 6
    BatchSize = len(AllInps) // NBatch

    for N in range(NBatch):

        with open(Path(__file__).parent / 'Medtool' / f'RunAbaqus_{N+1}.bash','w') as File:
            File.write('# Bash script to run abaqus simulations and get homogenized stress\n')

        for Input in AllInps[N*BatchSize:(N+1)*BatchSize]:

            # Define simulation type
            Type = Input.name.split('_')[-1]

            # If temporary file, remove it
            TempFile = str(Input)[:-len(Type)] + '_temp.inp'
            if Path(TempFile).exists():
                Path(TempFile).unlink()

            # Get image size
            Image = sitk.ReadImage(str(Input)[:-len(Type)-1] + '.mhd')
            Size = [int(Si * Sp) for Si,Sp in zip(Image.GetSize(), Image.GetSpacing())]

            # Write files
            WriteMain(Input.parent.name, Input, Size, Type)
            WriteBash(N+1, Input.parent.name, Input.name, Size, Type)


if __name__ == '__main__':
    
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Define default paths relative to this script's location
    ScriptDir = Path(__file__).parent
    DefaultPath = ScriptDir.parent / 'Results' / 'Cortical' / 'Homogenisation'

    # Add positional arguments with defaults
    Parser.add_argument('DataPath', nargs='?', default=str(DefaultPath), help='Path to scan files')

    # Parse arguments
    Arguments = Parser.parse_args()

    # Call main function with parsed or default arguments
    Main(Arguments.DataPath)

#%%
