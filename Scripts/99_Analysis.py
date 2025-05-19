#%% !/usr/bin/env python3

Description = """
Read results from homogenisation and fit
them to modified Zysset-Curnier model
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats.distributions import t

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def ProposedModel(Lambda0, Lambda0p, Mu0, k1, k2, l, Rho, eValues, eVectors=np.eye(3)):

    # k = [k1, k2, k3]
    # m1, m2, m3 = eVectors
    # M1, M2, M3 = np.outer(m1,m1), np.outer(m2,m2), np.outer(m3,m3)
    # M = [M1, M2, M3]
    # m1, m2, m3 = eValues
    # m = [m1,m2,m3]

    # S4 = np.zeros((3,3,3,3))
    # for i in range(3):
    #     S4 += (Lambda0 + 2*Mu0) * Rho**k[i] * m[i]**l * m[i]**l   * Tensor.Dyadic(M[i],M[i])
    #     for j in range(3):
    #         if i != j:
    #             S4 += Lambda0p * Rho**((k[i]+k[j])/2) * m[i]**l * m[j]**l * Tensor.Dyadic(M[i],M[j])
    #             S4 +=    2*Mu0 * Rho**((k[i]+k[j])/2) * m[i]**l * m[j]**l * Tensor.Symmetric(M[i],M[j])


    m1, m2, m3 = eVectors
    M1, M2, M3 = np.outer(m1,m1), np.outer(m2,m2), np.outer(m3,m3)
    M = [M1, M2, M3]
    m1, m2, m3 = eValues
    m = [m1,m2,m3]

    S4 = np.zeros((3,3,3,3))
    for i in range(3):
        S4 += Lambda0  * Rho**(k1*(m[i]+m[i])/(2*m[i]*m[i]))     * m[i]**l * m[i]**l   * Tensor.Dyadic(M[i],M[i])
        S4 += Lambda0p * Rho**(k1*(m[i]+m[i-1])/(2*m[i]*m[i-1])) * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i],M[i-1])
        S4 += Lambda0p * Rho**(k1*(m[i]+m[i-1])/(2*m[i]*m[i-1])) * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i-1],M[i])
        S4 +=   2*Mu0  * Rho**(k2*(m[i]+m[i])/(2*m[i]*m[i]))     * m[i]**l * m[i]**l   * Tensor.Symmetric(M[i],M[i])
        S4 +=   2*Mu0  * Rho**(k2*(m[i]+m[i-1])/(2*m[i]*m[i-1])) * m[i]**l * m[i-1]**l * Tensor.Symmetric(M[i-1],M[i])
    
    return S4

def Ref(Lambda0, Lambda0p, Mu0, k1, k2, l, Rho, eValues):

    m1, m2, m3 = eValues

    S = np.zeros((6,6))
    S[0,0] = m1**(2*l) * Rho**(k1/m1) * Lambda0 + 2*Mu0 * m1**(2*l) * Rho**(k2/m1)
    S[1,1] = m2**(2*l) * Rho**(k1/m2) * Lambda0 + 2*Mu0 * m2**(2*l) * Rho**(k2/m2)
    S[2,2] = m3**(2*l) * Rho**(k1/m3) * Lambda0 + 2*Mu0 * m3**(2*l) * Rho**(k2/m3)
    S[0,1] = m1**l * m2**l * Rho**(k1*(m1+m2)/(2*m1*m2)) * Lambda0p
    S[1,2] = m2**l * m3**l * Rho**(k1*(m2+m3)/(2*m2*m3)) * Lambda0p
    S[2,0] = m3**l * m1**l * Rho**(k1*(m3+m1)/(2*m3*m1)) * Lambda0p
    S[1,0] = m2**l * m1**l * Rho**(k1*(m2+m1)/(2*m2*m1)) * Lambda0p
    S[2,1] = m3**l * m2**l * Rho**(k1*(m3+m2)/(2*m3*m2)) * Lambda0p
    S[0,2] = m1**l * m3**l * Rho**(k1*(m1+m3)/(2*m1*m3)) * Lambda0p
    S[3,3] = 2*Mu0 * m2**l * m3**l * Rho**(k2*(m2+m3)/(2*m2*m3))
    S[4,4] = 2*Mu0 * m3**l * m1**l * Rho**(k2*(m3+m1)/(2*m3*m1))
    S[5,5] = 2*Mu0 * m1**l * m2**l * Rho**(k2*(m1+m2)/(2*m1*m2))

    return S

def PlotFitK(Y, X, B, FName=''):

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.exp(Y)
    Y_Fit = np.exp(X*B)
    Residuals = Y_Obs - Y_Fit

    # Compute R2adj and NE
    RSS = np.sum([R**2 for R in Residuals])
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-12-1)

    NE = np.zeros(len(Y_Obs))
    for i, s in enumerate(Y_Obs):
        Numerator = np.sum([T**2 for T in (Y_Obs[i]-Y_Fit[i])])
        Denominator = np.sum([T**2 for T in Y_Obs[i]])
        NE[i] = np.sqrt(Numerator/Denominator)

    # Prepare data for plot
    Line = np.linspace(min(Y_Obs[Y_Obs>0].min(), Y_Fit[Y_Fit>0].min()),
                       max(Y_Obs.max(), Y_Fit.max()), len(Y_Obs))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    for i in range(3):
        Axes.plot(Y_Obs[i*4::12,0], Y_Fit[i*4::12,0], color=Colors[0], linestyle='none', marker='s')
    for i in range(3):
        Axes.plot(Y_Obs[i+1::12,0], Y_Fit[i+1::12,0], color=Colors[1], linestyle='none', marker='o')
        Axes.plot(Y_Obs[i+5::12,0], Y_Fit[i+5::12,0], color=Colors[1], linestyle='none', marker='o')
    for i in range(3):
        Axes.plot(Y_Obs[i+9::12,0], Y_Fit[i+9::12,0], color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(Y_Obs)//12), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y_Obs)), xy=(0.3, 0.025), xycoords='axes fraction')
    Axes.annotate(r'$R^2_{ajd}$: ' + format(round(R2adj, 3),'.3f'), xy=(0.65, 0.1), xycoords='axes fraction')
    Axes.annotate(r'NE : ' + format(round(NE.mean(), 2), '.2f') + '$\pm$' + format(round(NE.std(), 2), '.2f'), xy=(0.65, 0.025), xycoords='axes fraction')
    Axes.set_xlabel('Observed $\mathrm{\mathbb{S}}$ (MPa)')
    Axes.set_ylabel('Fitted $\mathrm{\mathbb{S}}$ (MPa)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return (R2adj, NE)

#%% Main

def Main():

    # Define paths
    ResultsPath =Path(__file__).parents[1] / '_Results' / 'Cortical' / 'Homogenisation'
    FabricPath = Path(__file__).parents[1] / '_Results' / 'Cortical' / 'Fabric'
    Folders = [Folder.name for Folder in FabricPath.iterdir() if Folder.is_dir()]

    # List simulations output files
    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    # Collect data
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Isotropic = np.zeros((len(Folders), len(ROIs), 6, 6))
    Transverse = np.zeros((len(Folders), len(ROIs), 6, 6))
    Compliance = np.zeros((len(Folders), len(ROIs), 6, 6))
    Fabric = np.zeros((len(Folders), len(ROIs), 3))
    Rho = np.zeros((len(Folders), len(ROIs)))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get MIL results
            FabricData = Read.Fabric(FabricPath / Folder / (ROI+'.fab'))
            Fabric[f,r] = FabricData[0]
            Rho[f,r] = FabricData[2]

            # Get simulation results
            for S, sType in zip([Isotropic, Transverse],['Isotropic','Transverse']):

                # Get homogenization stress results
                File = open(ResultsPath / Folder / (ROI + f'_{sType}.out'), 'r').readlines()
                Stress = np.zeros((6,6))
                for i in range(6):
                    for j in range(6):
                        Stress[i,j] = float(File[i+4].split()[j+1])

                # Compute stiffness
                for i in range(6):
                    for j in range(6):
                        S[f,r,i,j] = Stress[i,j] / Strain[i]

                # Symetrize matrix
                S[f,r] = 1/2 * (S[f,r] + S[f,r].T)

                # Project onto orthotropy
                for i in range(6):
                    for j in range(6):
                        if i>2 or j>2:
                            if i != j:
                                S[f,r,i,j] = 0

                # Compute corresponding compliance
                Compliance[f,r] = np.linalg.inv(S[f,r])

    # Reshape arrays
    Isotropic = np.reshape(Isotropic, (-1, 6, 6))
    Transverse = np.reshape(Transverse, (-1, 6, 6))
    Compliance = np.reshape(Compliance, (-1, 6, 6))
    Fabric = np.reshape(Fabric, (-1, 3))
    Rho = np.reshape(Rho, (-1))

    # Build linear system
    Step = 12
    X = np.matrix(np.zeros((len(Rho)*Step, 13)))
    Y = np.matrix(np.zeros((len(Rho)*Step, 1)))
    m1, m2, m3 = np.mean(Fabric, axis=0)
    m1 = (m1 + m2) / 2
    m2 = m1
    for f in range(len(Rho)):
        
        Start, Stop = Step*f, Step*(f+1)
        
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, 0, np.log(m1 * m1)],
                                  [0, 1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, np.log(m1 * m2)],
                                  [0, 1, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, np.log(m1 * m3)],
                                  [0, 1, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, np.log(m2 * m1)],
                                  [1, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, np.log(m2 * m2)],
                                  [0, 1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, np.log(m2 * m3)],
                                  [0, 1, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, np.log(m3 * m1)],
                                  [0, 1, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, np.log(m3 * m2)],
                                  [1, 0, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, np.log(m3 * m3)],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, np.log(m2 * m3)],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, np.log(m3 * m1)],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), np.log(m1 * m2)]])
        
        Y[Start:Stop] = np.log([[Transverse[f][0,0]],
                                [Transverse[f][0,1]],
                                [Transverse[f][0,2]],
                                [Transverse[f][1,0]],
                                [Transverse[f][1,1]],
                                [Transverse[f][1,2]],
                                [Transverse[f][2,0]],
                                [Transverse[f][2,1]],
                                [Transverse[f][2,2]],
                                [Transverse[f][3,3]],
                                [Transverse[f][4,4]],
                                [Transverse[f][5,5]]])
    
    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    Lambda0 = np.exp(B[0,0]) - np.exp(B[2,0])
    Lambda0p = np.exp(B[1,0])
    Mu0 = np.exp(B[2,0]) / 2
    k11 = B[3,0]
    k12 = B[4,0]
    k13 = B[5,0]
    k21 = B[6,0]
    k22 = B[7,0]
    k23 = B[4,0]
    k31 = B[5,0]
    k32 = B[6,0]
    k33 = B[8,0]
    k44 = B[9,0]
    k55 = B[10,0]
    k66 = B[11,0]
    l = B[12,0]


    FName = Path(__file__).parents[1] / '_Results/ProposedModel.png'
    R2adj, NE = PlotFitK(Y, X, B)# FName=str(FName))

    print('Lamba0 = ', Lambda0)
    print('Lamba0p = ', Lambda0p)
    print('Mu0 = ', Mu0)
    print('k11 = ', k11)
    print('k12 = ', k12)
    print('k13 = ', k13)
    print('k21 = ', k21)
    print('k22 = ', k22)
    print('k23 = ', k23)
    print('k31 = ', k31)
    print('k32 = ', k32)
    print('k33 = ', k33)
    print('k44 = ', k44)
    print('k55 = ', k55)
    print('k66 = ', k66)
    

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
