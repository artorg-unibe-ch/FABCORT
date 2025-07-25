#%% !/usr/bin/env python3

Description = """
Read homogenization results and fit them to a modified Zysset-Curnier model for cortical and trabecular bone.

This script performs multiple linear regression to fit a modified Zysset-Curnier model to
homogenized stiffness tensors derived from finite element simulations. It evaluates cortical ROIs
using both isotropic and transverse isotropic assumptions, and compares them with trabecular ROIs.
The modification of the model consists of a split between bulk and shear components.

Functionality:
- Reads homogenized stiffness data (`Isotropic.csv`, `Transverse.csv`, and `ROIsData.csv`)
- Filters ROIs based on bone volume variation threshold (CV < 0.263)
- Constructs the design matrix and log-transformed stiffness response for each ROI
- Fits the modified Zysset-Curnier model via linear regression in log-log space
- Computes fit statistics: adjusted R², normalized error (NE), confidence intervals
- Generates scatter plots comparing observed vs fitted stiffness
- Exports a LaTeX-formatted table summarizing the fitted model parameters and metrics

Fitted parameters:
- λ₀, λ₀′, μ₀ (scaling coefficients)
- k (density exponent)
- l (fabric anisotropy exponent)

Inputs:
    - Cortical data:
        - `Results/Cortical/Isotropic.csv`
        - `Results/Cortical/Transverse.csv`
        - `Results/Cortical/CV.csv`
    - Trabecular data:
        - `Data_Trabecular/ROIsData.csv`

Outputs:
    - Parameter plots:
        - `Results/Panyasantisuk_Transverse.png`
        - `Results/Panyasantisuk_Isotropic.png`
        - `Results/Panyasantisukr_Trabecular.png`
    - LaTeX table of results:
        - `Results/Panyasantisuk.tex`

Example:
    python PanyasantisukModel.py
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

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def SplitHomogenization_OLS(X1, Y1, X2, Y2, Alpha=0.95, FName=''):

    # Solve linear systems
    XTXi = np.linalg.inv(X1.T * X1)
    B1 = XTXi * X1.T * Y1

    XTXi = np.linalg.inv(X2.T * X2)
    B2 = XTXi * X2.T * Y2

    # Compute predictions
    Y_Obs1 = np.exp(Y1)
    Y_Obs2 = np.exp(Y2)
    Y_Obs = np.vstack((Y_Obs1, Y_Obs2))
    Y_Fit1 = np.exp(X1 * B1)
    Y_Fit2 = np.exp(X2 * B2)
    Y_Fit = np.vstack((Y_Fit1, Y_Fit2))

    # Compute residuals, variance, and covariance matrix
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y_Obs) - 7

    # Store parameters in data frame
    Parameters = pd.DataFrame(columns=['Lambda0','Lambda0p','Mu0','kb','lb','ks','ls'])
    Parameters.loc['Value'] = [np.exp(B1[0,0]) - np.exp(B2[0,0]), np.exp(B1[1,0]), np.exp(B2[0,0])/2, B1[2,0], B1[3,0], B2[1,0], B2[2,0]]

    # Round values
    Columns = Parameters.columns
    for C in Columns[:3]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x))
    for C in Columns[3:]:
        Parameters[C] = Parameters[C].apply(lambda x: round(x,2))

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-7-1)

    NE = []
    for i in range(0,len(Y_Obs),12):
        T_Obs = Y_Obs[i:i+12]
        T_Fit = Y_Fit[i:i+12]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), Y_Fit.min()),
                       max(Y_Obs.max(), Y_Fit.max()), len(Y_Obs))

    # Plots
    DPI = 500
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(Y_Obs1[X1[:, 0] == 1], Y_Fit1[X1[:, 0] == 1],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs1[X1[:, 1] == 1], Y_Fit1[X1[:, 1] == 1],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs2[X2[:, 0] == 1], Y_Fit2[X2[:, 0] == 1],
              color=Colors[2], linestyle='none', marker='^')
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
    # Axes.set_xlim([SMin, SMax])
    # Axes.set_ylim([SMin, SMax])
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.close()

    return Parameters, R2adj, NE

#%% Main

def Main():

    # Read data
    CVData = pd.read_csv(Path(__file__).parents[1] / 'Results/Cortical/CV.csv', index_col=[0,1])
    IsoData = pd.read_csv(Path(__file__).parents[1] / 'Results/Cortical/Isotropic.csv')
    TransData = pd.read_csv(Path(__file__).parents[1] / 'Results/Cortical/Transverse.csv')
    TrabData = pd.read_csv(Path(__file__).parents[1] / 'Data_Trabecular/ROIsData.csv')

    # Store data into numpy arrays
    Isotropic = np.zeros((len(IsoData),6,6))
    Transverse = np.zeros((len(TransData),6,6))
    Trabecular = np.zeros((len(TrabData),6,6))
    for i in range(6):
        for j in range(6):
            if i == j:
                Isotropic[:,i,j] = IsoData[f'S{i+1}{j+1}']
                Transverse[:,i,j] = TransData[f'S{i+1}{j+1}']
                Trabecular[:,i,j] = TrabData[f'S{i+1}{j+1}']
            elif i < 3 and j < 3:
                Isotropic[:,i,j] = IsoData[f'S{i+1}{j+1}']
                Transverse[:,i,j] = TransData[f'S{i+1}{j+1}']
                Trabecular[:,i,j] = TrabData[f'S{i+1}{j+1}']

    # Remove ROIs with too high heterogeneity
    Threshold = 0.263
    FilterCort = CVData['CV'] < Threshold
    FilterTrab = TrabData['Variation Coefficient'] < Threshold
    RhoCort = IsoData['BV/TV'].values[FilterCort]
    FabricCort = IsoData[['m1','m2','m3']].values[FilterCort]
    RhoTrab = TrabData['BV/TV'].values[FilterTrab]
    FabricTrab = TrabData[['m1','m2','m3']].values[FilterTrab]
    Transverse = Transverse[FilterCort]
    Isotropic = Isotropic[FilterCort]
    Trabecular = Trabecular[FilterTrab]

    # Splitted fitting for shear and bulk components
    X1 = np.matrix(np.zeros((len(RhoCort)*9, 4)))
    Y1 = np.matrix(np.zeros((len(RhoCort)*9, 1)))
    X2 = np.matrix(np.zeros((len(RhoCort)*3, 3)))
    Y2 = np.matrix(np.zeros((len(RhoCort)*3, 1)))
    for f in range(len(RhoCort)):
        
        Start, Stop = 9*f, 9*(f+1)
        m1, m2, m3 = FabricCort[f]


        # Build system
        X1[Start:Stop] = np.array([[1, 0, np.log(RhoCort[f]), np.log(m1 ** 2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m1 * m2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m1 * m3)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m2 * m1)],
                                   [1, 0, np.log(RhoCort[f]), np.log(m2 ** 2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m2 * m3)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m3 * m1)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m3 * m2)],
                                   [1, 0, np.log(RhoCort[f]), np.log(m3 ** 2)]])
        
        Y1[Start:Stop] = np.log([[Transverse[f][0,0]],
                                 [Transverse[f][0,1]],
                                 [Transverse[f][0,2]],
                                 [Transverse[f][1,0]],
                                 [Transverse[f][1,1]],
                                 [Transverse[f][1,2]],
                                 [Transverse[f][2,0]],
                                 [Transverse[f][2,1]],
                                 [Transverse[f][2,2]]])
        
        Start, Stop = 3*f, 3*(f+1)
        
        X2[Start:Stop] = np.array([[1, np.log(RhoCort[f]), np.log(m2 * m3)],
                                   [1, np.log(RhoCort[f]), np.log(m3 * m1)],
                                   [1, np.log(RhoCort[f]), np.log(m1 * m2)]])
        
        Y2[Start:Stop] = np.log([[Transverse[f][3,3]],
                                 [Transverse[f][4,4]],
                                 [Transverse[f][5,5]/2]])
    
    FName = Path(__file__).parents[1] / 'Results/Panyasantisuk_Transverse.png'
    TransParameters, TransR2adj, TransNE = SplitHomogenization_OLS(X1, Y1, X2, Y2, FName=str(FName))
    print(TransParameters)

    # Splitted fitting for shear and bulk components
    X1 = np.matrix(np.zeros((len(RhoCort)*9, 4)))
    Y1 = np.matrix(np.zeros((len(RhoCort)*9, 1)))
    X2 = np.matrix(np.zeros((len(RhoCort)*3, 3)))
    Y2 = np.matrix(np.zeros((len(RhoCort)*3, 1)))
    for f in range(len(RhoCort)):
        
        Start, Stop = 9*f, 9*(f+1)
        m1, m2, m3 = FabricCort[f]


        # Build system
        X1[Start:Stop] = np.array([[1, 0, np.log(RhoCort[f]), np.log(m1 ** 2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m1 * m2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m1 * m3)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m2 * m1)],
                                   [1, 0, np.log(RhoCort[f]), np.log(m2 ** 2)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m2 * m3)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m3 * m1)],
                                   [0, 1, np.log(RhoCort[f]), np.log(m3 * m2)],
                                   [1, 0, np.log(RhoCort[f]), np.log(m3 ** 2)]])
        
        Y1[Start:Stop] = np.log([[Isotropic[f][0,0]],
                                 [Isotropic[f][0,1]],
                                 [Isotropic[f][0,2]],
                                 [Isotropic[f][1,0]],
                                 [Isotropic[f][1,1]],
                                 [Isotropic[f][1,2]],
                                 [Isotropic[f][2,0]],
                                 [Isotropic[f][2,1]],
                                 [Isotropic[f][2,2]]])
        
        Start, Stop = 3*f, 3*(f+1)
        
        X2[Start:Stop] = np.array([[1, np.log(RhoCort[f]), np.log(m2 * m3)],
                                   [1, np.log(RhoCort[f]), np.log(m3 * m1)],
                                   [1, np.log(RhoCort[f]), np.log(m1 * m2)]])
        
        Y2[Start:Stop] = np.log([[Isotropic[f][3,3]],
                                 [Isotropic[f][4,4]],
                                 [Isotropic[f][5,5]/2]])
        
    FName = Path(__file__).parents[1] / 'Results/Panyasantisuk_Isotropic.png'
    IsoParameters, IsoR2adj, IsoNE = SplitHomogenization_OLS(X1, Y1, X2, Y2, FName=str(FName))
    print(IsoParameters)
    
    # Splitted fitting for shear and bulk components
    X1 = np.matrix(np.zeros((len(RhoTrab)*9, 4)))
    Y1 = np.matrix(np.zeros((len(RhoTrab)*9, 1)))
    X2 = np.matrix(np.zeros((len(RhoTrab)*3, 3)))
    Y2 = np.matrix(np.zeros((len(RhoTrab)*3, 1)))
    for f in range(len(RhoTrab)):
        
        Start, Stop = 9*f, 9*(f+1)
        m1, m2, m3 = FabricTrab[f]


        # Build system
        X1[Start:Stop] = np.array([[1, 0, np.log(RhoTrab[f]), np.log(m1 ** 2)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m1 * m2)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m1 * m3)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m2 * m1)],
                                   [1, 0, np.log(RhoTrab[f]), np.log(m2 ** 2)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m2 * m3)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m3 * m1)],
                                   [0, 1, np.log(RhoTrab[f]), np.log(m3 * m2)],
                                   [1, 0, np.log(RhoTrab[f]), np.log(m3 ** 2)]])
        
        Y1[Start:Stop] = np.log([[Trabecular[f][0,0]],
                                 [Trabecular[f][0,1]],
                                 [Trabecular[f][0,2]],
                                 [Trabecular[f][1,0]],
                                 [Trabecular[f][1,1]],
                                 [Trabecular[f][1,2]],
                                 [Trabecular[f][2,0]],
                                 [Trabecular[f][2,1]],
                                 [Trabecular[f][2,2]]])
        
        Start, Stop = 3*f, 3*(f+1)
        
        X2[Start:Stop] = np.array([[1, np.log(RhoTrab[f]), np.log(m2 * m3)],
                                   [1, np.log(RhoTrab[f]), np.log(m3 * m1)],
                                   [1, np.log(RhoTrab[f]), np.log(m1 * m2)]])
        
        Y2[Start:Stop] = np.log([[Trabecular[f][3,3]],
                                 [Trabecular[f][4,4]],
                                 [Trabecular[f][5,5]/2]])
    
    FName = Path(__file__).parents[1] / 'Results/Panyasantisuk_Trabecular.png'
    TrabParameters, TrabR2adj, TrabNE = SplitHomogenization_OLS(X1, Y1, X2, Y2, FName=str(FName))
    print(TrabParameters)

    # Write Latex table
    Parameters = pd.DataFrame()
    Parameters.loc['Transverse','$\lambda_0$'] = str(TransParameters.loc['Value', 'Lambda0'])
    Parameters.loc['Transverse','$\lambda_0\'$'] = str(TransParameters.loc['Value', 'Lambda0p'])
    Parameters.loc['Transverse','$\mu_0$'] = str(TransParameters.loc['Value', 'Mu0'])
    Parameters.loc['Transverse','$k$ Bulk'] = str(TransParameters.loc['Value', 'kb'])
    Parameters.loc['Transverse','$l$ Bulk'] = str(TransParameters.loc['Value', 'lb'])
    Parameters.loc['Transverse','$k$ Shear'] = str(TransParameters.loc['Value', 'ks'])
    Parameters.loc['Transverse','$l$ Shear'] = str(TransParameters.loc['Value', 'ls'])
    Parameters.loc['Transverse','$R^2_{adj}$'] = str(round(TransR2adj, 2))
    Parameters.loc['Transverse','NE'] = str(round(np.mean(TransNE),2))

    Parameters.loc['Isotropic','$\lambda_0$'] = str(IsoParameters.loc['Value', 'Lambda0'])
    Parameters.loc['Isotropic','$\lambda_0\'$'] = str(IsoParameters.loc['Value', 'Lambda0p'])
    Parameters.loc['Isotropic','$\mu_0$'] = str(IsoParameters.loc['Value', 'Mu0'])
    Parameters.loc['Isotropic','$k$ Bulk'] = str(IsoParameters.loc['Value', 'kb'])
    Parameters.loc['Isotropic','$l$ Bulk'] = str(IsoParameters.loc['Value', 'lb'])
    Parameters.loc['Isotropic','$k$ Shear'] = str(IsoParameters.loc['Value', 'ks'])
    Parameters.loc['Isotropic','$l$ Shear'] = str(IsoParameters.loc['Value', 'ls'])
    Parameters.loc['Isotropic','$R^2_{adj}$'] = str(round(IsoR2adj,2))
    Parameters.loc['Isotropic','NE'] = str(round(np.mean(IsoNE),2))

    Parameters.loc['Trabecular','$\lambda_0$'] = str(TrabParameters.loc['Value', 'Lambda0'])
    Parameters.loc['Trabecular','$\lambda_0\'$'] = str(TrabParameters.loc['Value', 'Lambda0p'])
    Parameters.loc['Trabecular','$\mu_0$'] = str(TrabParameters.loc['Value', 'Mu0'])
    Parameters.loc['Trabecular','$k$ Bulk'] = str(TrabParameters.loc['Value', 'kb'])
    Parameters.loc['Trabecular','$l$ Bulk'] = str(TrabParameters.loc['Value', 'lb'])
    Parameters.loc['Trabecular','$k$ Shear'] = str(TrabParameters.loc['Value', 'ks'])
    Parameters.loc['Trabecular','$l$ Shear'] = str(TrabParameters.loc['Value', 'ls'])
    Parameters.loc['Trabecular','$R^2_{adj}$'] = str(round(TrabR2adj,2))
    Parameters.loc['Trabecular','NE'] = str(round(np.mean(TrabNE),2))

    Parameters.index.name = 'Model'
    Parameters = Parameters.reset_index()

    Table = Path(__file__).parents[1] / 'Results' / 'Panyasantisuk.tex'
    Caption = 'Parameters and fit quality coefficient obtained with the modified Zysset-Curnier model.'
    Label = 'TabPanyasantisuk'
    cFormat  = 'l|c|c|c|c|c|c|c'
    Parameters.to_latex(Table, caption=Caption,
                        column_format=cFormat,
                        label=Label, escape=False,
                        index=False, position='!h')

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
