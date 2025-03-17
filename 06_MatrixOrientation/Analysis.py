#%% !/usr/bin/env python3

Description = """
Investigate the effect of different matrix orientations
(fixed vs around pores)
"""

__author__ = ['Mathieu Simon']
__date_created__ = '03-03-2024'
__date__ = '03-03-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, morphology, measure


#%% Functions


#%% Main

def Main():

    # Read images
    Stained = io.imread('Stained.png')
    Segmented = io.imread('Segmented.png')

    # Crop images
    X0, X1 = 70, 390
    Y0, Y1 = 220, 460
    Stained = Stained[Y0:Y1,X0:X1]
    Segmented = Segmented[Y0:Y1,X0:X1]

    # Smooth cement line mask
    Disk = morphology.disk(40)
    Mask = Segmented[:,:,2] == 255
    Mask = morphology.binary_dilation(Mask, Disk)
    Mask = morphology.binary_erosion(Mask, Disk)

    # Compute isodistances from canal to cement line
    Label = Mask*1 + (Segmented[:,:,1] == 255)*1
    MedialAxis, Distances = morphology.medial_axis(1-Label, return_distance=True)
    Contours = measure.find_contours(Distances, level=1) # Equals to Distance == 1

    Figure, Axis = plt.subplots(1,1)
    Axis.imshow(Stained)
    for C in Contours:
        Axis.plot(C[:,1],C[:,0],color=(1,0,0))
    # Axis.imshow(Segmented)
    plt.show(Figure)

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
