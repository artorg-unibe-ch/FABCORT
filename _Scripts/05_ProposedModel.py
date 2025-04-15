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

def ProposedModel(Lambda0, Lambda0p, Mu0, k1, k2, k3, l, Rho, eValues, eVectors):

    k = [k1, k2, k3]
    m1, m2, m3 = eVectors
    M1, M2, M3 = np.outer(m1,m1), np.outer(m2,m2), np.outer(m3,m3)
    M = [M1, M2, M3]
    m1, m2, m3 = eValues
    m = [m1,m2,m3]

    S4 = np.zeros((3,3,3,3))
    for i in range(3):
        S4 += (Lambda0 + 2*Mu0) * Rho**k[i] * m[i]**l * m[i]**l   * Tensor.Dyadic(M[i],M[i])
        for j in range(3):
            if i != j:
                S4 += Lambda0p * Rho**((k[i]+k[j])/2) * m[i]**l * m[j]**l * Tensor.Dyadic(M[i],M[j])
                S4 +=    2*Mu0 * Rho**((k[i]+k[j])/2) * m[i]**l * m[j]**l * Tensor.Symmetric(M[i],M[j])
    
    return S4

def Homogenization_OLS(X, Y, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.exp(Y)
    Y_Fit = np.exp(X * B)
    Residuals = Y - X*B
    DOFs = len(Y) - X.shape[1]
    Sigma = Residuals.T * Residuals / DOFs
    Cov = Sigma[0,0] * XTXi

    # Compute B confidence interval
    t_Alpha = t.interval(Alpha, DOFs)
    B_CI_Low = B.T + t_Alpha[0] * np.sqrt(np.diag(Cov))
    B_CI_Top = B.T + t_Alpha[1] * np.sqrt(np.diag(Cov))

    # Store parameters in data frame
    Parameters = pd.DataFrame(columns=['Lambda0','Lambda0p','Mu0','k1','k2','k3','l'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0], B[5,0], B[6,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Top[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4], B_CI_Low[0,5], B_CI_Low[0,6]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Low[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4], B_CI_Top[0,5], B_CI_Top[0,6]]

    # Round values
    Columns = Parameters.columns
    for C in Columns[:3]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x))
    for C in Columns[3:]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x,2))

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y - Y.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y)-1)/(len(Y)-X.shape[1]-1)

    NE = []
    for i in range(0,len(Y),12):
        T_Obs = Y_Obs[i:i+12]
        T_Fit = Y_Fit[i:i+12]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y.min(), (X*B).min()),
                       max(Y.max(), (X*B).max()), len(Y))

    # Plots
    DPI = 500
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(Y_Obs[X[:, 0] == 1], Y_Fit[X[:, 0] == 1],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[X[:, 1] == 1], Y_Fit[X[:, 1] == 1],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[X[:, 2] == 1], Y_Fit[X[:, 2] == 1],
              color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(np.exp(Line), np.exp(Line), color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(Y)//12), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y)), xy=(0.3, 0.025), xycoords='axes fraction')
    Axes.annotate(r'$R^2_{ajd}$: ' + format(round(R2adj, 3),'.3f'), xy=(0.65, 0.1), xycoords='axes fraction')
    Axes.annotate(r'NE : ' + format(round(NE.mean(), 2), '.2f') + '$\pm$' + format(round(NE.std(), 2), '.2f'), xy=(0.65, 0.025), xycoords='axes fraction')
    Axes.set_xlabel('Observed $\mathrm{\mathbb{S}}$ (MPa)')
    Axes.set_ylabel('Fitted $\mathrm{\mathbb{S}}$ (MPa)')
    # Axes.set_xlim([SMin, SMax])
    # Axes.set_ylim([SMin, SMax])
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return Parameters, R2adj, NE

def Homogenization_OLS2(X, Y, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.exp(Y)
    Y_Fit = np.exp(X * B)
    Residuals = Y - X*B
    DOFs = len(Y) - X.shape[1]
    Sigma = Residuals.T * Residuals / DOFs
    Cov = Sigma[0,0] * XTXi

    # Compute B confidence interval
    t_Alpha = t.interval(Alpha, DOFs)
    B_CI_Low = B.T + t_Alpha[0] * np.sqrt(np.diag(Cov))
    B_CI_Top = B.T + t_Alpha[1] * np.sqrt(np.diag(Cov))

    # Store parameters in data frame
    Parameters = pd.DataFrame(columns=['Lambda0','Lambda0p','Mu0','k11','k12','k13','k23','k22','k33','k44','k55','k66','l'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0], B[5,0], B[6,0], B[7,0], B[8,0], B[9,0], B[10,0], B[11,0], B[12,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Top[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4], B_CI_Low[0,5], B_CI_Low[0,6], B_CI_Low[0,7], B_CI_Low[0,8], B_CI_Low[0,9], B_CI_Low[0,10], B_CI_Low[0,11], B_CI_Low[0,12]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Low[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4], B_CI_Top[0,5], B_CI_Top[0,6], B_CI_Top[0,7], B_CI_Top[0,8], B_CI_Top[0,9], B_CI_Top[0,10], B_CI_Top[0,11], B_CI_Top[0,12]]

    # Round values
    Columns = Parameters.columns
    for C in Columns[:3]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x))
    for C in Columns[3:]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x,2))

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y - Y.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y)-1)/(len(Y)-X.shape[1]-1)

    NE = []
    for i in range(0,len(Y),12):
        T_Obs = Y_Obs[i:i+12]
        T_Fit = Y_Fit[i:i+12]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y.min(), (X*B).min()),
                       max(Y.max(), (X*B).max()), len(Y))

    # Plots
    DPI = 500
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(Y_Obs[X[:, 0] == 1], Y_Fit[X[:, 0] == 1],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[X[:, 1] == 1], Y_Fit[X[:, 1] == 1],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[X[:, 2] == 1], Y_Fit[X[:, 2] == 1],
              color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(np.exp(Line), np.exp(Line), color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(Y)//12), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y)), xy=(0.3, 0.025), xycoords='axes fraction')
    Axes.annotate(r'$R^2_{ajd}$: ' + format(round(R2adj, 3),'.3f'), xy=(0.65, 0.1), xycoords='axes fraction')
    Axes.annotate(r'NE : ' + format(round(NE.mean(), 2), '.2f') + '$\pm$' + format(round(NE.std(), 2), '.2f'), xy=(0.65, 0.025), xycoords='axes fraction')
    Axes.set_xlabel('Observed $\mathrm{\mathbb{S}}$ (MPa)')
    Axes.set_ylabel('Fitted $\mathrm{\mathbb{S}}$ (MPa)')
    # Axes.set_xlim([SMin, SMax])
    # Axes.set_ylim([SMin, SMax])
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return Parameters, R2adj, NE

#%% Main

def Main():

    # Define paths
    ResultsPath =Path(__file__).parents[1] / '_Results' / 'Homogenisation'
    FabricPath = Path(__file__).parents[1] / '_Results' / 'Fabric'
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

    # Test 
    S4 = ProposedModel(3, 2, 1.5, 2, 2, 1, 0.5, 0.8, [0.9,0.9,1.2], np.eye(3))
    S = Tensor.IsoMorphism3333_66(S4)
    print(S)

    # Fit homogenization with theorical model
    X = np.matrix(np.zeros((len(Rho)*12, 7)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    for f in range(len(Rho)):
        
        Start, Stop = 12*f, 12*(f+1)
        m1, m2, m3 = Fabric[f]

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0,   np.log(Rho[f]),                0,                0, np.log(m1 * m1)],
                                  [0, 1, 0, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)],
                                  [1, 0, 0,                0,   np.log(Rho[f]),                0, np.log(m2 * m2)],
                                  [0, 1, 0,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 1, 0,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [1, 0, 0,                0,                0,   np.log(Rho[f]), np.log(m3 * m3)],
                                  [0, 0, 1,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [0, 0, 1, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 0, 1, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)]])
        
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
    
    FName = Path(__file__).parents[1] / '_Results/ProposedModel.png'
    Parameters, R2adj, NE = Homogenization_OLS(X, Y, FName=str(FName))
    print(Parameters)

    # Fit homogenization with theorical model and average fabric
    X = np.matrix(np.zeros((len(Rho)*12, 7)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    m1, m2, m3 = np.mean(Fabric, axis=0)
    for f in range(len(Rho)):
        
        Start, Stop = 12*f, 12*(f+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0,   np.log(Rho[f]),                0,                0, np.log(m1 * m1)],
                                  [0, 1, 0, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)],
                                  [1, 0, 0,                0,   np.log(Rho[f]),                0, np.log(m2 * m2)],
                                  [0, 1, 0,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 1, 0,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [1, 0, 0,                0,                0,   np.log(Rho[f]), np.log(m3 * m3)],
                                  [0, 0, 1,                0, np.log(Rho[f])/2, np.log(Rho[f])/2, np.log(m2 * m3)],
                                  [0, 0, 1, np.log(Rho[f])/2,                0, np.log(Rho[f])/2, np.log(m1 * m3)],
                                  [0, 0, 1, np.log(Rho[f])/2, np.log(Rho[f])/2,                0, np.log(m1 * m2)]])
        
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
    
    Parameters, R2adj, NE = Homogenization_OLS(X, Y)
    print(Parameters)

    # Fit homogenization with other theorical model
    X = np.matrix(np.zeros((len(Rho)*12, 13)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    m1, m2, m3 = np.mean(Fabric, axis=0)
    for f in range(len(Rho)):
        
        Start, Stop = 12*f, 12*(f+1)
        # m1, m2, m3 = Fabric[f]

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, 0, np.log(m1 * m1)],
                                  [0, 1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, np.log(m1 * m2)],
                                  [0, 1, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, np.log(m1 * m3)],
                                  [0, 1, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, np.log(m1 * m2)],
                                  [1, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, np.log(m2 * m2)],
                                  [0, 1, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, 0, np.log(m2 * m3)],
                                  [0, 1, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, 0, np.log(m1 * m3)],
                                  [0, 1, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, 0, 0, np.log(m2 * m3)],
                                  [1, 0, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, 0, np.log(m3 * m3)],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, 0, np.log(m2 * m3)],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, np.log(Rho[f]), 0, np.log(m1 * m3)],
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
    
    FName = Path(__file__).parents[1] / '_Results/FullKModel_AverageFabric.png'
    Parameters, R2adj, NE = Homogenization_OLS2(X, Y, FName=str(FName))
    print(Parameters)
    
    # Fit compliance with 2 constants model

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
