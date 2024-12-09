#%% !/usr/bin/env python3

Description = """
Class to read different type of files
"""

__author__ = ['Mathieu Simon']
__date_created__ = '05-12-2024'
__date__ = '05-12-2024'
__license__ = 'MIT'
__version__ = '1.0'


import numpy as np

class Read():

    def __init__(self):
        pass

    def Fabric(self, FileName):

        Text = open(FileName,'r').readlines()
        BVTV = float(Text[12].split('=')[1])
        eValues = np.array(Text[18].split(':')[1].split(),float)
        eVectors = np.zeros((3,3))
        for i in range(3):
            eVectors[i] = Text[19+i].split(':')[1].split()

        # Sort eigen values and eigen vectors
        Args = np.argsort(eValues)
        eValues = eValues[Args]
        eVectors = eVectors[Args]

        return eValues, eVectors, BVTV
