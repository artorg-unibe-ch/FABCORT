#%% !/usr/bin/env python3

Description = """
Read results from homogenisation and fit
them to model proposed by Cowin and Yang (1999)
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
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def Q66(Q):

    Q11, Q12, Q13 = Q[0,0]**2, Q[0,1]**2, Q[0,2]**2
    Q21, Q22, Q23 = Q[1,0]**2, Q[1,1]**2, Q[1,2]**2
    Q31, Q32, Q33 = Q[2,0]**2, Q[2,1]**2, Q[2,2]**2

    Q14, Q15, Q16 = Q[0,1]*Q[0,2], Q[0,0]*Q[0,2], Q[0,0]*Q[0,1]
    Q24, Q25, Q26 = Q[1,1]*Q[1,2], Q[1,0]*Q[1,2], Q[1,1]*Q[1,0]
    Q34, Q35, Q36 = Q[2,2]*Q[2,1], Q[2,2]*Q[2,0], Q[2,0]*Q[2,1]

    Q41, Q51, Q61 = Q[1,0]*Q[2,0], Q[1,1]*Q[2,1], Q[1,2]*Q[2,2]
    Q42, Q52, Q62 = Q[0,0]*Q[2,0], Q[0,1]*Q[2,1], Q[0,2]*Q[2,2]
    Q43, Q53, Q63 = Q[0,0]*Q[2,0], Q[0,1]*Q[1,1], Q[0,2]*Q[1,2]

    Q44 = Q[1,1]*Q[2,2] + Q[1,2]*Q[2,1]
    Q45 = Q[1,0]*Q[2,2] + Q[2,0]*Q[1,2]
    Q46 = Q[1,0]*Q[2,1] + Q[2,0]*Q[1,1]

    Q54 = Q[0,1]*Q[2,2] + Q[2,1]*Q[0,2]
    Q55 = Q[0,0]*Q[2,2] + Q[0,2]*Q[2,1]
    Q56 = Q[0,0]*Q[2,1] + Q[2,0]*Q[0,1]

    Q64 = Q[0,1]*Q[1,2] + Q[1,1]*Q[0,2]
    Q65 = Q[0,0]*Q[1,2] + Q[1,0]*Q[0,2]
    Q66 = Q[0,0]*Q[1,1] + Q[1,0]*Q[0,1]

    Q = np.array([[Q11, Q12, Q13, Q14*2**0.5, Q15*2**0.5, Q16*2**0.5],
                  [Q21, Q22, Q23, Q24*2**0.5, Q25*2**0.5, Q26*2**0.5],
                  [Q31, Q32, Q33, Q34*2**0.5, Q35*2**0.5, Q36*2**0.5],
                  [Q41*2**0.5, Q42*2**0.5, Q43*2**0.5, Q44, Q45, Q46],
                  [Q51*2**0.5, Q52*2**0.5, Q53*2**0.5, Q54, Q55, Q56],
                  [Q61*2**0.5, Q62*2**0.5, Q63*2**0.5, Q64, Q65, Q66]])
    
    return Q

def FitEigenValues(Rho, eValues):

    X = np.ones((len(Rho), 2))
    X[:,1] = np.log(Rho)
    X = np.matrix(X)

    Y = np.matrix(np.log(eValues)).T

    # Solve least-squares in batch: X = (A^T A)^{-1} A^T y
    XTXi = np.linalg.inv(X.T @ X)
    B = XTXi @ X.T @ Y

    return (np.exp(B[0,0]), B[1,0])

def PlotFit(S, St, FName=''):

    # Compute residuals, variance, and covariance matrix
    Y_Obs = S
    Y_Fit = St
    Residuals = Y_Obs - Y_Fit

    # Compute R2adj and NE
    RSS = np.sum([R**2 for R in Residuals])
    TSS = np.sum([R**2 for R in (S[S > 0] - S[S > 0].mean())])
    R2adj = 1 - RSS/TSS * (len(S)-1)/(len(S)-6-1)

    NE = np.zeros(len(S))
    for i, s in enumerate(S):
        Numerator = np.sum([T**2 for T in (Y_Obs[i]-Y_Fit[i])])
        Denominator = np.sum([T**2 for T in Y_Obs[i]])
        NE[i] = np.sqrt(Numerator/Denominator)

    # Prepare data for plot
    Line = np.linspace(min(Y_Obs[Y_Obs>0].min(), Y_Fit[Y_Fit>0].min()),
                       max(Y_Obs.max(), Y_Fit.max()), len(S))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    for i in range(3):
        Axes.plot(Y_Obs[:,i,i], Y_Fit[:,i,i],
                  color=Colors[0], linestyle='none', marker='s')
    for i in range(3):
        for j in range(3):
            if i != j:
                Axes.plot(Y_Obs[:,i,j], Y_Fit[:,i,j],
                          color=Colors[1], linestyle='none', marker='o')
    for i in range(3):
        Axes.plot(Y_Obs[:,i+3,i+3], Y_Fit[:,i+3,i+3],
                  color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(S)), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(S)*12), xy=(0.3, 0.025), xycoords='axes fraction')
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

    # Read CV data
    CVData = pd.read_csv(Path(__file__).parents[1] / '_Results' / 'Cortical' / 'CV.csv', index_col=[0,1])

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
    CV = np.zeros((len(Folders), len(ROIs)))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get sample CV
            CV[f,r] = CVData.loc[(Folder, ROI),'CV']

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

                # Compute corresponding compliance
                Compliance[f,r] = np.linalg.inv(S[f,r])

    # Reshape arrays
    Isotropic = np.reshape(Isotropic, (-1, 6, 6))
    Transverse = np.reshape(Transverse, (-1, 6, 6))
    Compliance = np.reshape(Compliance, (-1, 6, 6))
    Fabric = np.reshape(Fabric, (-1, 3))
    Rho = np.reshape(Rho, (-1))
    CV = np.reshape(CV, (-1))

    # Remove ROIs with too high heterogeneity
    Filter = CV < 0.263
    Rho = Rho[Filter]
    Fabric = Fabric[Filter]
    Transverse = Transverse[Filter]
    Isotropic = Isotropic[Filter]
    Compliance = Compliance[Filter]

    # Permute tensor coordinates to always have the same structure
    Idx = np.where(Transverse[:,0,0] > Transverse[:,1,1])[0]
    Q = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    Q = Q66(Q)
    for i in Idx:
        Transverse[i] = Q @ Transverse[i] @ Q.T
        Isotropic[i] = Q @ Isotropic[i] @ Q.T
        Compliance[i] = Q @ Compliance[i] @ Q.T

    # Project tensors onto transverse isotropy
    Projected, NE = Tensor.TransProj(Transverse)

    # Compute eigenvalues
    Alpha = np.arctan(2**0.5 / (4*Projected[:,0,2]) * (Projected[:,0,0] + Projected[:,0,1] - Projected[:,2,2]))
    L1 = Projected[:,2,2] + 2**0.5 * Projected[:,0,2] * (np.tan(Alpha) + 1/np.cos(Alpha))
    L2 = Projected[:,2,2] + 2**0.5 * Projected[:,0,2] * (np.tan(Alpha) - 1/np.cos(Alpha))
    L34 = Projected[:,0,0] - Projected[:,0,1]
    L56 = Projected[:,3,3]
    Lambdas = [L1, L2, L34, L56]

    # Compute eigentensors
    E11 = 1/4 * (1 + np.sin(Alpha)) * (Projected[:,0,0] + Projected[:,1,1])
    E11 += 2**0.5/4 * np.cos(Alpha) * Projected[:,2,2]
    E33 = 2**0.5/4 * np.cos(Alpha) * (Projected[:,0,0] + Projected[:,1,1])
    E33 += 1/2 * (1 - np.sin(Alpha)) * Projected[:,2,2]
    E1 = np.zeros((len(L1), 3, 3))
    E1[:,0,0] = E11
    E1[:,1,1] = E11
    E1[:,2,2] = E33
    E1 = (E1.T / np.linalg.norm(E1, axis=(1,2))).T

    E11 = 1/4 * (1 - np.sin(Alpha)) * (Projected[:,0,0] + Projected[:,1,1])
    E11 -= 2**0.5/4 * np.cos(Alpha) * Projected[:,2,2]
    E33 = -2**0.5/4 * np.cos(Alpha) * (Projected[:,0,0] + Projected[:,1,1])
    E33 += 1/2 * (1 + np.sin(Alpha)) * Projected[:,2,2]
    E2 = np.zeros((len(L2), 3, 3))
    E2[:,0,0] = E11
    E2[:,1,1] = E11
    E2[:,2,2] = E33
    E2 = (E2.T / np.linalg.norm(E2, axis=(1,2))).T

    E3 = np.zeros((len(L34), 3, 3))
    E3[:,0,0] = Projected[:,0,0]
    E3[:,1,1] = -Projected[:,1,1]
    E3 = (E3.T / np.linalg.norm(E3, axis=(1,2))).T

    E4 = np.zeros((len(L34), 3, 3))
    E4[:,0,1] = Projected[:,0,1]
    E4[:,1,0] = -Projected[:,1,0]
    E4 = (E4.T / np.linalg.norm(E4, axis=(1,2))).T

    E5 = np.zeros((len(L56), 3, 3))
    E5[:,0,2] = Projected[:,0,2]
    E5[:,2,0] = Projected[:,2,0]
    E5 = (E5.T / np.linalg.norm(E5, axis=(1,2))).T

    E6 = np.zeros((len(L56), 3, 3))
    E6[:,1,2] = Projected[:,1,2]
    E6[:,2,1] = Projected[:,2,1]
    E6 = (E6.T / np.linalg.norm(E6, axis=(1,2))).T

    # Plot eigen values as function of BVTV
    X = np.linspace(Rho.min(), Rho.max(), 100)
    Figure, Axis = plt.subplots(2,2, sharex=True, figsize=(10,7))
    for i in range(2):
        for j in range(2):
            Y = Lambdas[2*i+j]/1E3
            P = FitEigenValues(Rho, Y)
            YFit = P[0] * X ** P[1]
            Residuals = Y - (P[0] * Rho ** P[1])
            RSS = np.sum([R**2 for R in Residuals])
            TSS = np.sum([R**2 for R in (Y - Y.mean())])
            RegSS = TSS - RSS
            R2 = RegSS / TSS
            Axis[i,j].plot(Rho, Y, marker='o', linestyle='none', color=(1,0,0))
            Axis[i,j].plot(X, YFit, color=(0,0,1))
            Axis[i,j].set_xlabel(r'$\rho$ (-)')
            Axis[i,j].set_ylabel(f'$\Lambda_{2*i+j+1}$ (GPa)')
            Axis[i,j].annotate(f'$\Lambda_{2*i+j+1}$ = {round(P[0])}' + r'$\rho$' + '$^{%.2f}$'%(round(P[1],2)), xy=(0.1, 0.9), xycoords='axes fraction')
            Axis[i,j].annotate(r'$R^2$: ' + format(round(R2, 3),'.3f'), xy=(0.1, 0.8), xycoords='axes fraction')
    plt.tight_layout()
    plt.show(Figure)

    # Corresponding bases in 6-dimensional space
    N1 = np.array([E1[:,0,0],
                   E1[:,1,1],
                   E1[:,2,2],
                   E1[:,1,2],
                   E1[:,0,2],
                   E1[:,0,1]]).T
    
    N2 = np.array([E2[:,0,0],
                   E2[:,1,1],
                   E2[:,2,2],
                   E2[:,1,2],
                   E2[:,0,2],
                   E2[:,0,1]]).T
    
    N3 = np.array([E3[:,0,0],
                   E3[:,1,1],
                   E3[:,2,2],
                   E3[:,1,2],
                   E3[:,0,2],
                   E3[:,0,1]]).T

    N4 = np.array([E4[:,0,0],
                   E4[:,1,1],
                   E4[:,2,2],
                   E4[:,1,2],
                   E4[:,0,2],
                   E4[:,0,1]]).T  * 2**0.5
    
    N5 = np.array([E5[:,0,0],
                   E5[:,1,1],
                   E5[:,2,2],
                   E5[:,1,2],
                   E5[:,0,2],
                   E5[:,0,1]]).T * 2**0.5

    N6 = np.array([E6[:,0,0],
                   E6[:,1,1],
                   E6[:,2,2],
                   E6[:,1,2],
                   E6[:,0,2],
                   E6[:,0,1]]).T * 2**0.5

    Bases = np.array([-N1.T,N2.T,-N3.T,N4.T,N5.T,N6.T]).T

    # Arithmetic average
    N_NA = np.sum(Bases, axis=0) / len(Bases)

    # J
    J_NA = np.zeros((6,6))
    for i in range(6):
        J_NA += np.outer(N_NA[i], N_NA[i])
    L, D = np.linalg.eig(J_NA)
    J_NA = D @ np.diag(L**(-1/2)) @ D.T

    # Average basis
    N_AVG = np.zeros((6,6))
    for i in range(6):
        N_AVG[i] = J_NA @ N_NA[i]

    # Build again tensor
    L = [L1, L2, L34, L34, L56, L56]
    Yang = np.zeros(Projected.shape)
    for i in range(6):

        # Basis
        NN = np.ones(Yang.shape) * np.outer(N_AVG[i], N_AVG[i])

        # Modulus
        P = FitEigenValues(Rho, L[i])
        M = np.expand_dims(P[0] * Rho**P[1],-1)

        Yang += np.expand_dims(M,-1) * NN

    FName = Path(__file__).parents[1] / '_Results/YangCowinModel.png'
    R2, NE = PlotFit(Projected, Yang, FName=str(FName))

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
