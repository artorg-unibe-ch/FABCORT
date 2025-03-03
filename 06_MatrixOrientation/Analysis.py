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

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


#%% Functions


#%% Main

def Main():

    # Define grid
    nx, ny = 100, 100  # Grid size
    x = np.linspace(-2, 2, nx)
    y = np.linspace(-2, 2, ny)
    X, Y = np.meshgrid(x, y)

    # Define parameters
    U_inf = 1  # Free-stream velocity
    R = 0.5  # Radius of the disk
    Xc, Yc = 0, 0  # Center of the disk

    # Convert to polar coordinates
    r = np.sqrt((X - Xc)**2 + (Y - Yc)**2)
    theta = np.arctan2(Y - Yc, X - Xc)

    # Define velocity components using potential flow theory
    u = U_inf * (1 - (R**2 / r**2)) * np.cos(theta)
    v = U_inf * (1 + (R**2 / r**2)) * np.sin(theta)

    # Mask values inside the disk
    u[r < R] = np.nan
    v[r < R] = np.nan

    # Plot streamlines
    plt.figure(figsize=(6, 6))
    plt.streamplot(X, Y, u, v, density=2, color='b')
    circle = plt.Circle((Xc, Yc), R, color='k', fill=True)
    plt.gca().add_patch(circle)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Streamlines Around a Circular Obstacle')
    plt.show()


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
