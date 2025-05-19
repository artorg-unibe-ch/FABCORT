#%% !/usr/bin/env python3

Description = """
Reproduce healthy group results from:
Simon et al. 2022
Fabric-elasticity relationships of tibial trabecular bone are similar
in osteogenesis imperfecta and healthy individuals
Bone 155
https://doi.org/10.1016/j.bone.2021.116282
"""

__author__ = ['Mathieu Simon']
__date_created__ = '29-11-2024'
__date__ = '29-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats.distributions import norm, t

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

#%% Functions

def BuildSystem(Data):

    Y = np.zeros((len(Data)*12, 1))
    Y[::12,0] = np.log(Data['S11'].values)
    Y[1::12,0] = np.log(Data['S12'].values)
    Y[2::12,0] = np.log(Data['S13'].values)
    Y[3::12,0] = np.log(Data['S21'].values)
    Y[4::12,0] = np.log(Data['S22'].values)
    Y[5::12,0] = np.log(Data['S23'].values)
    Y[6::12,0] = np.log(Data['S31'].values)
    Y[7::12,0] = np.log(Data['S32'].values)
    Y[8::12,0] = np.log(Data['S33'].values)
    Y[9::12,0] = np.log(Data['S44'].values)
    Y[10::12,0] = np.log(Data['S55'].values)
    Y[11::12,0] = np.log(Data['S66'].values)

    X = np.zeros((len(Data)*12, 5))
    X[::12,0] = np.ones(len(Data))
    X[1::12,1] = np.ones(len(Data))
    X[2::12,1] = np.ones(len(Data))
    X[3::12,1] = np.ones(len(Data))
    X[4::12,0] = np.ones(len(Data))
    X[5::12,1] = np.ones(len(Data))
    X[6::12,1] = np.ones(len(Data))
    X[7::12,1] = np.ones(len(Data))
    X[8::12,0] = np.ones(len(Data))
    X[9::12,2] = np.ones(len(Data))
    X[10::12,2] = np.ones(len(Data))
    X[11::12,2] = np.ones(len(Data))
    X[:,3] = np.repeat(np.log(Data['BV/TV']),12)
    X[::12,4] = np.log(Data['m1'] ** 2)
    X[1::12,4] = np.log(Data['m1'] * Data['m2'])
    X[2::12,4] = np.log(Data['m1'] * Data['m3'])
    X[3::12,4] = np.log(Data['m2'] * Data['m1'])
    X[4::12,4] = np.log(Data['m2'] ** 2)
    X[5::12,4] = np.log(Data['m2'] * Data['m3'])
    X[6::12,4] = np.log(Data['m3'] * Data['m1'])
    X[7::12,4] = np.log(Data['m3'] * Data['m2'])
    X[8::12,4] = np.log(Data['m3'] ** 2)
    X[9::12,4] = np.log(Data['m2'] * Data['m3'])
    X[10::12,4] = np.log(Data['m3'] * Data['m1'])
    X[11::12,4] = np.log(Data['m1'] * Data['m2'])

    return np.matrix(X), np.matrix(Y)

def OLS(X, Y, Alpha=0.95):

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
    # B_0 = np.sort(np.sqrt(np.diag(X * Cov * X.T)))
    # CI_Line_u = np.exp(Line + t_Alpha[0] * B_0)
    # CI_Line_o = np.exp(Line + t_Alpha[1] * B_0)

    # Plots
    DPI = 500
    SMax = max([Y_Obs.max(), Y_Fit.max()]) * 5
    SMin = min([Y_Obs.min(), Y_Fit.min()]) / 5
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    # Set boundaries of fabtib
    SMax = 1e4
    SMin = 1e-3

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
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
    Axes.set_xlim([SMin, SMax])
    Axes.set_ylim([SMin, SMax])
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.show()

    return Parameters, R2adj, NE


#%% Main

def Main():

    # Read data
    Data = pd.read_csv('Data.csv')

    # Build matrices
    X, Y = BuildSystem(Data)

    # Perform linear regression
    Parameters, R2adj, NE = OLS(X,Y)

    # Filter data set
    Threshold = 0.263
    Filter = Data['Variation Coefficient'] < Threshold
    FData = Data[Filter]

    # Build matrices
    F_X, F_Y = BuildSystem(FData)

    # Perform linear regression
    F_Parameters, F_R2adj, F_NE = OLS(F_X,F_Y)
    print(F_Parameters.loc['Value'].round(2))

    # Compute mean and standard deviations
    Mean = round(FData['BV/TV'].mean(),3)
    Std = round(FData['BV/TV'].std(),3)
    print(f'BV/TV Mean $\pm$ SD: {Mean} $\pm$ {Std}')

    Mean = round(FData['m1'].mean(),3)
    Std = round(FData['m1'].std(),3)
    print(f'm1 Mean $\pm$ SD: {Mean} $\pm$ {Std}')

    Mean = round(FData['m2'].mean(),3)
    Std = round(FData['m2'].std(),3)
    print(f'm1 Mean $\pm$ SD: {Mean} $\pm$ {Std}')

    Mean = round(FData['m3'].mean(),3)
    Std = round(FData['m3'].std(),3)
    print(f'm1 Mean $\pm$ SD: {Mean} $\pm$ {Std}')

    Mean = round((FData['m3'] / FData['m1']).mean(),3)
    Std = round((FData['m3'] / FData['m1']).std(),3)
    print(f'DA Mean $\pm$ SD: {Mean} $\pm$ {Std}')



if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)
