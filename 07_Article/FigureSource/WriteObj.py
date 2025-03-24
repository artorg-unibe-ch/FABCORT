#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used to write .obj files from uCT scans
    These files can then be imported in blender for rendering

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: November 2023
    """

#%% Imports
# Modules import

import time
import argparse
import numpy as np
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt
from skimage import measure

#%% Classes
# Define classes
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
            print('')
            self.Tic = time.time()
            self.Update(0, Text)

        elif StartStop*1 == 0:
            self.Update(1, Text)
            self.Print()

Time = Time()

#%% Functions
# Define functions

def WriteOBJ(Array, FileName):

    Time.Process(1, 'Get surfaces')

    # Write obj file
    Vertices, Faces, Normals, Values = measure.marching_cubes(Array)
    Faces = Faces + 1

    Time.Update(2/5, 'Write vertices')

    OBJ = open(FileName, 'w')
    for V in Vertices:
        OBJ.write('v {0} {1} {2}\n'.format(V[0],V[1],V[2]))

    Time.Update(3/5, 'Write normals')

    for N in Normals:
        OBJ.write('vn {0} {1} {2}\n'.format(N[0],N[1],N[2]))

    Time.Update(4/5, 'Write faces')

    for F in Faces:
        OBJ.write('f {0}//{0} {1}//{1} {2}//{2}\n'.format(F[0],F[1],F[2]))

    OBJ.close()

    Time.Process(0, 'File written')

    return

#%% Main
# Main code

def Main():

    CWD = Path(__file__).parent
    File = CWD / 'CorticalROI.mhd'

    Scan = sitk.ReadImage(str(File))
    Array = sitk.GetArrayFromImage(Scan)

    # Threshold implant and bone (only slice for faster computation)
    Otsu = sitk.OtsuMultipleThresholdsImageFilter()
    Otsu.SetNumberOfThresholds(1)
    Otsu.Execute(Scan)
    Thresholds = Otsu.GetThresholds()

    # Get binary images
    Bone = Array > Thresholds[0]

    # Write .obj files
    FileName = File.parent / (File.name[:-4] + '.obj')
    WriteOBJ(np.pad(Bone,1), FileName)

    return

#%% Execution part
# Execution as main
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main()