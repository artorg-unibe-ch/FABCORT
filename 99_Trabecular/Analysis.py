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

def AverageBasis(Stiffness):

    # Compute eigenvalues
    Alpha = np.arctan(2**0.5 / (4*Stiffness[:,0,2]) * (Stiffness[:,0,0] + Stiffness[:,0,1] - Stiffness[:,2,2]))
    L1 = Stiffness[:,2,2] + 2**0.5 * Stiffness[:,0,2] * (np.tan(Alpha) + 1/np.cos(Alpha))
    L2 = Stiffness[:,2,2] + 2**0.5 * Stiffness[:,0,2] * (np.tan(Alpha) - 1/np.cos(Alpha))
    L34 = Stiffness[:,0,0] - Stiffness[:,0,1]
    L56 = Stiffness[:,3,3]
    Lambdas = [L1, L2, L34, L34, L56, L56]

    # Compute eigentensors
    E11 = 1/4 * (1 + np.sin(Alpha)) * (Stiffness[:,0,0] + Stiffness[:,1,1])
    E11 += 2**0.5/4 * np.cos(Alpha) * Stiffness[:,2,2]
    E33 = 2**0.5/4 * np.cos(Alpha) * (Stiffness[:,0,0] + Stiffness[:,1,1])
    E33 += 1/2 * (1 - np.sin(Alpha)) * Stiffness[:,2,2]
    E1 = np.zeros((len(L1), 3, 3))
    E1[:,0,0] = E11
    E1[:,1,1] = E11
    E1[:,2,2] = E33
    E1 = (E1.T / np.linalg.norm(E1, axis=(1,2))).T

    E11 = 1/4 * (1 - np.sin(Alpha)) * (Stiffness[:,0,0] + Stiffness[:,1,1])
    E11 -= 2**0.5/4 * np.cos(Alpha) * Stiffness[:,2,2]
    E33 = -2**0.5/4 * np.cos(Alpha) * (Stiffness[:,0,0] + Stiffness[:,1,1])
    E33 += 1/2 * (1 + np.sin(Alpha)) * Stiffness[:,2,2]
    E2 = np.zeros((len(L2), 3, 3))
    E2[:,0,0] = E11
    E2[:,1,1] = E11
    E2[:,2,2] = E33
    E2 = (E2.T / np.linalg.norm(E2, axis=(1,2))).T

    E3 = np.zeros((len(L34), 3, 3))
    E3[:,0,0] = Stiffness[:,0,0]
    E3[:,1,1] = -Stiffness[:,1,1]
    E3 = (E3.T / np.linalg.norm(E3, axis=(1,2))).T

    E4 = np.zeros((len(L34), 3, 3))
    E4[:,0,1] = Stiffness[:,0,1]
    E4[:,1,0] = -Stiffness[:,1,0]
    E4 = (E4.T / np.linalg.norm(E4, axis=(1,2))).T

    E5 = np.zeros((len(L56), 3, 3))
    E5[:,0,2] = Stiffness[:,0,2]
    E5[:,2,0] = Stiffness[:,2,0]
    E5 = (E5.T / np.linalg.norm(E5, axis=(1,2))).T

    E6 = np.zeros((len(L56), 3, 3))
    E6[:,1,2] = Stiffness[:,1,2]
    E6[:,2,1] = Stiffness[:,2,1]
    E6 = (E6.T / np.linalg.norm(E6, axis=(1,2))).T

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

    return Lambdas, N_AVG

def FitEigenValues(Rho, eValues):

    X = np.ones((len(Rho), 2))
    X[:,1] = np.log(Rho)
    X = np.matrix(X)

    Y = np.matrix(np.log(eValues)).T

    # Solve least-squares in batch: X = (A^T A)^{-1} A^T y
    XTXi = np.linalg.inv(X.T @ X)
    B = XTXi @ X.T @ Y

    return (np.exp(B[0,0]), B[1,0])

def PlotFit(S, St, L, FName=''):

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

    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # for i in range(3):
    #     Axes.plot(Y_Obs[:,i,i], Y_Fit[:,i,i],
    #               color=Colors[0], linestyle='none', marker='s')
    # for i in range(3):
    #     for j in range(3):
    #         if i != j:
    #             Axes.plot(Y_Obs[:,i,j], Y_Fit[:,i,j],
    #                       color=Colors[1], linestyle='none', marker='o')
    # for i in range(3):
    #     Axes.plot(Y_Obs[:,i+3,i+3], Y_Fit[:,i+3,i+3],
    #               color=Colors[2], linestyle='none', marker='^')
    # Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    # Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    # Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    # Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
    # Axes.annotate(r'N ROIs   : ' + str(len(S)), xy=(0.3, 0.1), xycoords='axes fraction')
    # Axes.annotate(r'N Points : ' + str(len(S)*12), xy=(0.3, 0.025), xycoords='axes fraction')
    # Axes.annotate(r'$R^2_{ajd}$: ' + format(round(R2adj, 3),'.3f'), xy=(0.65, 0.1), xycoords='axes fraction')
    # Axes.annotate(r'NE : ' + format(round(NE.mean(), 2), '.2f') + '$\pm$' + format(round(NE.std(), 2), '.2f'), xy=(0.65, 0.025), xycoords='axes fraction')
    # Axes.set_xlabel('Observed $\mathrm{\mathbb{S}}$ (MPa)')
    # Axes.set_ylabel('Fitted $\mathrm{\mathbb{S}}$ (MPa)')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend(loc='upper left')
    # plt.subplots_adjust(left=0.15, bottom=0.15)
    # if len(FName) > 0:
    #     plt.savefig(FName)
    # plt.show()

    Figure, Axes = plt.subplots(1, 3, figsize=(15,4.5), dpi=DPI, sharex=True, sharey=True)
    for i in range(3):
        Axes[0].plot(Y_Obs[:L,i,i], Y_Fit[:L,i,i], fillstyle='none',
                  color=(0,0,1), linestyle='none', marker='o')
        Axes[0].plot(Y_Obs[L:,i,i], Y_Fit[L:,i,i], fillstyle='none',
                  color=(1,0,0), linestyle='none', marker='o')
        Axes[0].set_title('$\lambda_{ii}$')
    for i in range(3):
        for j in range(3):
            if i != j:
                Axes[1].plot(Y_Obs[:L,i,j], Y_Fit[:L,i,j], fillstyle='none',
                          color=(0,0,1), linestyle='none', marker='o')
                Axes[1].plot(Y_Obs[L:,i,j], Y_Fit[L:,i,j], fillstyle='none',
                          color=(1,0,0), linestyle='none', marker='o')
                Axes[1].set_title('$\lambda_{ij}$')
    for i in range(3):
        Axes[2].plot(Y_Obs[:L,i+3,i+3], Y_Fit[:L,i+3,i+3], fillstyle='none',
                  color=(0,0,1), linestyle='none', marker='o')
        Axes[2].plot(Y_Obs[L:,i+3,i+3], Y_Fit[L:,i+3,i+3], fillstyle='none',
                  color=(1,0,0), linestyle='none', marker='o')
        Axes[2].set_title('$\mu_{ij}$')
    for i in range(3):
        Axes[i].plot([],fillstyle='none', color=(0,0,1), linestyle='none', marker='o', label='Trabecular')
        Axes[i].plot([],fillstyle='none', color=(1,0,0), linestyle='none', marker='o', label='Cortical')
        Axes[i].legend()
        Axes[i].plot(Line, Line, color=(0, 0, 0), linestyle='--')
        Axes[i].set_xlabel('Observed $\mathrm{\mathbb{S}}$ (MPa)')
    Axes[0].set_ylabel('Fitted $\mathrm{\mathbb{S}}$ (MPa)')
    plt.xscale('log')
    plt.yscale('log')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()



    return (R2adj, NE)

#%% Main

def Main():

    # Read trabecular data
    Data = pd.read_csv(Path(__file__).parent / 'Data.csv')

    # Filter ROIs with too high heterogeneity
    Filter = Data['Variation Coefficient'] < 0.263
    Data = Data[Filter]

    # Get volume fraction
    TrabRho = Data['BV/TV']

    # Get stiffness
    TrabStiffness = np.zeros((len(Data),6,6))
    TrabStiffness[:,0,0] = Data['S11']
    TrabStiffness[:,0,1] = Data['S12']
    TrabStiffness[:,0,2] = Data['S13']
    TrabStiffness[:,1,0] = Data['S21']
    TrabStiffness[:,1,1] = Data['S22']
    TrabStiffness[:,1,2] = Data['S23']
    TrabStiffness[:,2,0] = Data['S31']
    TrabStiffness[:,2,1] = Data['S32']
    TrabStiffness[:,2,2] = Data['S33']
    TrabStiffness[:,3,3] = Data['S44']
    TrabStiffness[:,4,4] = Data['S55']
    TrabStiffness[:,5,5] = Data['S66']

    # Read cortical CV data
    CVData = pd.read_csv(Path(__file__).parents[1] / 'Results' / 'CV.csv', index_col=[0,1])

    # Define paths
    ResultsPath =Path(__file__).parents[1] / 'Results' / 'Homogenisation'
    FabricPath = Path(__file__).parents[1] / 'Results' / 'Fabric'
    Folders = [Folder.name for Folder in FabricPath.iterdir() if Folder.is_dir()]

    # List simulations output files
    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    # Collect data
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    CortStiffness = np.zeros((len(Folders), len(ROIs), 6, 6))
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

            # Get homogenization stress results
            File = open(ResultsPath / Folder / (ROI + f'_Isotropic.out'), 'r').readlines()
            Stress = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stress[i,j] = float(File[i+4].split()[j+1])

            # Compute stiffness
            for i in range(6):
                for j in range(6):
                    CortStiffness[f,r,i,j] = Stress[i,j] / Strain[i]

            # Symetrize matrix
            CortStiffness[f,r] = 1/2 * (CortStiffness[f,r] + CortStiffness[f,r].T)

    # Reshape arrays
    CV = np.reshape(CV, (-1))
    Rho = np.reshape(Rho, (-1))
    CortStiffness = np.reshape(CortStiffness, (-1, 6, 6))

    # Remove ROIs with too high heterogeneity
    Filter = CV < 0.263
    CortRho = Rho[Filter]
    CortStiffness = CortStiffness[Filter]

    # Transform 4th rank into 2nd rank
    C = np.zeros(CortStiffness.shape)
    for i in range(3):
        for j in range(3):
            C[:,i,j] = CortStiffness[:,i,j]

    for i in range(3):
        for j in range(3,6):
            C[:,i,j] = 2**0.5 * CortStiffness[:,i,j]

    for i in range(3,6):
        for j in range(3):
            C[:,i,j] = 2**0.5 * CortStiffness[:,i,j]

    for i in range(3):
        for j in range(3):
            C[:,i+3,j+3] = 2 * CortStiffness[:,i+3,j+3]

    CortStiffness = C

    # Project into transverse isotropy
    TrabStiffness, NE = Tensor.TransProj(TrabStiffness)
    CortStiffness, NE = Tensor.TransProj(CortStiffness)

    # Group data
    Tensors = np.concatenate([TrabStiffness, CortStiffness], axis=0)
    XCort = np.linspace(CortRho.min(), CortRho.max(), 100)
    XTrab = np.linspace(TrabRho.min(), TrabRho.max(), 100)
    Rhos = np.hstack([TrabRho, CortRho])

    # Compute average basis
    V, B = AverageBasis(Tensors)

    # Build again tensor
    Yang = np.zeros(Tensors.shape)
    for i in range(6):

        # Basis
        NN = np.ones(Tensors.shape) * np.outer(B[i], B[i])

        # Modulus
        P = FitEigenValues(Rhos, V[i])
        M = np.expand_dims(P[0] * Rhos**P[1],-1)

        Yang += np.expand_dims(M,-1) * NN

    FName = str(Path(__file__).parent / 'SpectralModel.png')
    R2, NE = PlotFit(Tensors, Yang, len(TrabStiffness), FName)

    # Plot eigen values as function of BVTV
    X = np.hstack([XTrab, XCort])
    Rhos = np.hstack([TrabRho, CortRho])
    Figure, Axis = plt.subplots(2,3, sharex=True, figsize=(15,7))
    for i in range(2):
        for j in range(3):                
            Y = V[2*i+j] / 1E3
            P = FitEigenValues(Rhos, Y)
            YFit = P[0] * X ** P[1]
            Residuals = Y - (P[0] * Rhos ** P[1])
            RSS = np.sum([R**2 for R in Residuals])
            TSS = np.sum([R**2 for R in (Y - Y.mean())])
            RegSS = TSS - RSS
            R2 = RegSS / TSS
            Axis[i,j].plot(Rhos[:len(TrabRho)], Y[:len(TrabRho)], marker='o', linestyle='none', color=(0,0,1), fillstyle='none')
            Axis[i,j].plot(Rhos[len(TrabRho):], Y[len(TrabRho):], marker='o', linestyle='none', color=(1,0,0), fillstyle='none')
            Axis[i,j].plot(X, YFit, color=(0,0,0))
            Axis[i,j].set_xlabel(r'$\rho$ (-)')
            Axis[i,j].set_ylabel(f'$\Lambda_{2*i+j+1}$ (GPa)')
            Axis[i,j].annotate(f'$\Lambda_{2*i+j+1}$ = {round(P[0])}' + r'$\rho$' + '$^{%.2f}$'%(round(P[1],2)), xy=(0.1, 0.9), xycoords='axes fraction')
            Axis[i,j].annotate(r'$R^2$: ' + format(round(R2, 3),'.3f'), xy=(0.1, 0.8), xycoords='axes fraction')
            # Axis[i,j].set_xscale('log')
            # Axis[i,j].set_yscale('log')
    plt.tight_layout()
    plt.savefig(Path(__file__).parent / 'LambdaVSRho.png')
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
