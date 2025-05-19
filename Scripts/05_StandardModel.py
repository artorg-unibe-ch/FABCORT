#%% !/usr/bin/env python3

Description = """
Read results from homogenisation and fit
them to standard Zysset-Curnier model
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
    Parameters = pd.DataFrame(columns=['Lambda0','Lambda0p','Mu0','k','l'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0])/2, B[3,0], B[4,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - np.exp(B_CI_Top[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2])/2, B_CI_Low[0,3], B_CI_Low[0,4]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - np.exp(B_CI_Low[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2])/2, B_CI_Top[0,3], B_CI_Top[0,4]]

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

    # Read CV data
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
    Isotropic = np.zeros((len(Folders), len(ROIs), 6, 6))
    Transverse = np.zeros((len(Folders), len(ROIs), 6, 6))
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

    # Reshape arrays
    Isotropic = np.reshape(Isotropic, (-1, 6, 6))
    Transverse = np.reshape(Transverse, (-1, 6, 6))
    Fabric = np.reshape(Fabric, (-1, 3))
    Rho = np.reshape(Rho, (-1))
    CV = np.reshape(CV, (-1))

    # Remove ROIs with too low volume fraction
    Filter = CV < 0.263
    Rho = Rho[Filter]
    Fabric = Fabric[Filter]
    Transverse = Transverse[Filter]
    Isotropic = Isotropic[Filter]

    # Project onto transverse isotropy
    pIsotropic, NEIso = Tensor.TransProj(Isotropic)
    pTransverse, NETra = Tensor.TransProj(Transverse)

    # Fit homogenization with orthotropic theorical model
    X = np.matrix(np.zeros((len(Rho)*12, 5)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    for f in range(len(Rho)):
        
        Start, Stop = 12*f, 12*(f+1)
        m1, m2, m3 = Fabric[f]

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]), np.log(m1 ** 2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m1)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m2 ** 2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m2)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m3 ** 2)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m1 * m2)]])
        
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
    
    FName = Path(__file__).parents[1] / 'Results/StandardOrthoModel.png'
    OrthoParameters, OrthoR2adj, OrthoNE = Homogenization_OLS(X, Y, FName=str(FName))
    print(OrthoParameters)
    
    # Fit homogenization with theorical model using average transverse M
    X = np.matrix(np.zeros((len(Rho)*12, 5)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    m1, m2, m3 = np.mean(Fabric, axis=0)
    m1 = (m1 + m2) / 2
    m2 = m1
    for f in range(len(Rho)):
        
        Start, Stop = 12*f, 12*(f+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]), np.log(m1 ** 2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m1)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m2 ** 2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m2)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m3 ** 2)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m1 * m2)]])
        
        Y[Start:Stop] = np.log([[pTransverse[f][0,0]],
                                [pTransverse[f][0,1]],
                                [pTransverse[f][0,2]],
                                [pTransverse[f][1,0]],
                                [pTransverse[f][1,1]],
                                [pTransverse[f][1,2]],
                                [pTransverse[f][2,0]],
                                [pTransverse[f][2,1]],
                                [pTransverse[f][2,2]],
                                [pTransverse[f][3,3]],
                                [pTransverse[f][4,4]],
                                [pTransverse[f][5,5]]])
    
    FName = Path(__file__).parents[1] / 'Results/StandardTransModel.png'
    TransParameters, TransR2adj, TransNE = Homogenization_OLS(X, Y, FName=str(FName))
    print(TransParameters)

    # Write Latex table
    Parameters = pd.DataFrame()
    Parameters.loc['Orthotropic','$\lambda_0$'] = str(OrthoParameters.loc['Value', 'Lambda0'])
    Parameters.loc['Orthotropic','$\lambda_0\'$'] = str(OrthoParameters.loc['Value', 'Lambda0p'])
    Parameters.loc['Orthotropic','$\mu_0$'] = str(OrthoParameters.loc['Value', 'Mu0'])
    Parameters.loc['Orthotropic','$k$'] = str(OrthoParameters.loc['Value', 'k'])
    Parameters.loc['Orthotropic','$l$'] = str(OrthoParameters.loc['Value', 'l'])
    Parameters.loc['Orthotropic','$R^2_{adj}$'] = str(round(OrthoR2adj, 2))
    Parameters.loc['Orthotropic','NE'] = str(round(np.mean(OrthoNE),2))

    Parameters.loc['Transverse','$\lambda_0$'] = str(TransParameters.loc['Value', 'Lambda0'])
    Parameters.loc['Transverse','$\lambda_0\'$'] = str(TransParameters.loc['Value', 'Lambda0p'])
    Parameters.loc['Transverse','$\mu_0$'] = str(TransParameters.loc['Value', 'Mu0'])
    Parameters.loc['Transverse','$k$'] = str(TransParameters.loc['Value', 'k'])
    Parameters.loc['Transverse','$l$'] = str(TransParameters.loc['Value', 'l'])
    Parameters.loc['Transverse','$R^2_{adj}$'] = str(round(TransR2adj,2))
    Parameters.loc['Transverse','NE'] = str(round(np.mean(TransNE),2))

    Parameters.index.name = 'Model'
    Parameters = Parameters.reset_index()

    Table = Path(__file__).parents[1] / 'Article' / 'TabZysset.tex'
    Caption = 'Parameters and fit quality coefficient obtained with the standard Zysset-Curnier model in orthotropic and transverse isotropic space.'
    Label = 'TabZysset'
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
