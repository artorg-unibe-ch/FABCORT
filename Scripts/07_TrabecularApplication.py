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
from scipy.stats import norm
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def PlotHistogram(Variable:np.array, FName=''):

    # Get data attributes
    SortedValues = np.sort(Variable).astype(float)
    N = len(Variable)
    X_Bar = np.mean(Variable)
    S_X = np.std(Variable, ddof=1)
    Q025 = np.quantile(Variable, 0.25)
    Q075 = np.quantile(Variable, 0.75)

    Histogram, Edges = np.histogram(Variable, bins=20)
    Width = (Edges[1] - Edges[0])
    Center = (Edges[:-1] + Edges[1:]) / 2

    # Kernel density estimation (Gaussian kernel)
    KernelEstimator = np.zeros(N)
    NormalIQR = np.sum(np.abs(norm.ppf(np.array([0.25,0.75]), 0, 1)))
    DataIQR = np.abs(Q075) - np.abs(Q025)
    KernelHalfWidth = 0.9*N**(-1/5) * min(S_X,DataIQR/NormalIQR)
    for Value in SortedValues:
        Norm = norm.pdf(SortedValues-Value,loc=0,scale=KernelHalfWidth*2)
        KernelEstimator += Norm / max(Norm)
    # KernelEstimator = KernelEstimator/N*max(Histogram)

    # Scale values to sum up to 1
    Histogram = Histogram / np.sum(Histogram * Width) * Width
    Deltas = SortedValues[1:] - SortedValues[:-1]
    KernelMax = np.max([KernelEstimator[1:], KernelEstimator[:-1]], axis=0)
    KernelEstimator = KernelEstimator / np.sum(KernelMax * Deltas) * Width

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=200)
    Axes.plot(SortedValues, KernelEstimator,color=(1,0,0))
    Axes.bar(Center, Histogram, align='center', width=Width,edgecolor=(0,0,0),color=(1,1,1,0))
    plt.xlabel(r'$\rho$ (-)')
    plt.ylabel('Frequency (-)')
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return

def PlotOldFit(Y, X, B, FName=''):

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

def PlotFit(S, Lambda0, Lambda0p, Mu0, kl, km, l, Rho, eValues, FName=''):

    # Compute residuals, variance, and covariance matrix
    Y_Obs = S
    Y_Fit = np.zeros(S.shape)
    for i, (rho, eVals) in enumerate(zip(Rho, eValues)):
        S4 = ProposedModel(Lambda0, Lambda0p, Mu0, kl, km, l, rho, eVals)
        Y_Fit[i] = Tensor.IsoMorphism3333_66(S4)
    
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

    # Read ROIs data
    Data = pd.read_csv( Path(__file__).parents[1] / '_Data' / 'Trabecular' / 'ROIsData.csv')

    # Keep ROIs with CV < 0.263
    F = Data['Variation Coefficient'] < 0.263
    Files = []
    for I, Row in Data[F].iterrows():
        N = Row['ROI Number']
        S = Row['Scan']
        Files.append(f'{N}_{S}')

    # Define paths
    ResultsPath =Path(__file__).parents[1] / '_Results' / 'Trabecular' / 'Homogenisation'
    FabricPath = Path(__file__).parents[1] / '_Results' / 'Trabecular' / 'Fabric'

    # Read fabric data
    eValues = np.zeros((len(Files),3))
    Rho = np.zeros((len(Files)))
    eVectors = []
    for i, F in enumerate(Files):
        Fabric = Read.Fabric(FabricPath / (F+'.fab'))
        eValues[i] = Fabric[0]
        eVectors.append(Fabric[1])
        Rho[i] = Fabric[2]

    # Read simulations results
    Stiffness = np.zeros((len(Files),6,6))
    Compliance = np.zeros((len(Files),6,6))
    for i, F in enumerate(Files):

        # Read compliance results from medtool
        Compliance[i] = Read.ComplianceDat(ResultsPath / (F+'.dat'))
        Compliance[i] = (Compliance[i]+Compliance[i].T)/2
        
        # Compute corresponding stiffness
        Stiffness[i] = np.linalg.inv(Compliance[i])
        Stiffness[i] = (Stiffness[i]+Stiffness[i].T)/2

        # Write tensor into mandel notation
        Mandel = Tensor.Engineering2MandelNotation(Stiffness[i])

        # Transform tensor into fabric coordinate system
        I = np.eye(3)
        Q = eVectors[i]
        Transformed = Tensor.TransformTensor(Mandel, I, Q)

        # Project onto orthotropy
        for j in range(6):
            for k in range(6):
                if j>2 or k>2:
                    if j != k:
                        Transformed[j,k] = 0

        # Get tensor back to engineering notation
        Stiffness[i] = Tensor.Mandel2EngineeringNotation(Transformed)

        # Compute back corresponding compliance
        Compliance[i] = np.linalg.inv(Stiffness[i])


    # Reproduce old results
    Step = 12
    X = np.matrix(np.zeros((len(Rho)*Step, 5)))
    Y = np.matrix(np.zeros((len(Rho)*Step, 1)))
    # m1, m2, m3 = np.mean(eValues, axis=0)
    for f in range(len(Rho)):

        m1, m2, m3 = eValues[f]
        
        Start, Stop = Step*f, Step*(f+1)
        
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]), np.log(m1 * m1)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m1)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m2 * m2)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 1, 0, np.log(Rho[f]), np.log(m3 * m2)],
                                  [1, 0, 0, np.log(Rho[f]), np.log(m3 * m3)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m2 * m3)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m3 * m1)],
                                  [0, 0, 1, np.log(Rho[f]), np.log(m1 * m2)]])
        
        Y[Start:Stop] = np.log([[Stiffness[f][0,0]],
                                [Stiffness[f][0,1]],
                                [Stiffness[f][0,2]],
                                [Stiffness[f][1,0]],
                                [Stiffness[f][1,1]],
                                [Stiffness[f][1,2]],
                                [Stiffness[f][2,0]],
                                [Stiffness[f][2,1]],
                                [Stiffness[f][2,2]],
                                [Stiffness[f][3,3]],
                                [Stiffness[f][4,4]],
                                [Stiffness[f][5,5]]])
    
    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    Lambda0 = np.exp(B[0,0]) - np.exp(B[2,0])
    Lambda0p = np.exp(B[1,0])
    Mu0 = np.exp(B[2,0]) / 2
    k = B[3,0]
    l = B[4,0]

    FName = Path(__file__).parents[1] / '_Results/StandardModel_Trabecular.png'
    R2adj, NE = PlotOldFit(Y, X, B, FName=str(FName))
    
    # Test new model
    Step = 12
    X = np.matrix(np.zeros((len(Rho)*Step, 6)))
    Y = np.matrix(np.zeros((len(Rho)*Step, 1)))
    for f in range(len(Rho)):

        m1, m2, m3 = eValues[f]
        
        Start, Stop = Step*f, Step*(f+1)
        
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[f]) * (m1+m1)/(2*m1*m1), 0, np.log(m1 * m1)],
                                  [0, 1, 0, np.log(Rho[f]) * (m1+m2)/(2*m1*m2), 0, np.log(m1 * m2)],
                                  [0, 1, 0, np.log(Rho[f]) * (m1+m3)/(2*m1*m3), 0, np.log(m1 * m3)],
                                  [0, 1, 0, np.log(Rho[f]) * (m2+m1)/(2*m2*m1), 0, np.log(m2 * m1)],
                                  [1, 0, 0, np.log(Rho[f]) * (m2+m2)/(2*m2*m2), 0, np.log(m2 * m2)],
                                  [0, 1, 0, np.log(Rho[f]) * (m2+m3)/(2*m2*m3), 0, np.log(m2 * m3)],
                                  [0, 1, 0, np.log(Rho[f]) * (m3+m1)/(2*m3*m1), 0, np.log(m3 * m1)],
                                  [0, 1, 0, np.log(Rho[f]) * (m3+m2)/(2*m3*m2), 0, np.log(m3 * m2)],
                                  [1, 0, 0, np.log(Rho[f]) * (m3+m3)/(2*m3*m3), 0, np.log(m3 * m3)],
                                  [0, 0, 1, 0, np.log(Rho[f]) * (m2+m3)/(2*m2*m3), np.log(m2 * m3)],
                                  [0, 0, 1, 0, np.log(Rho[f]) * (m3+m1)/(2*m3*m1), np.log(m3 * m1)],
                                  [0, 0, 1, 0, np.log(Rho[f]) * (m1+m2)/(2*m1*m2), np.log(m1 * m2)]])
        
        Y[Start:Stop] = np.log([[Stiffness[f][0,0]],
                                [Stiffness[f][0,1]],
                                [Stiffness[f][0,2]],
                                [Stiffness[f][1,0]],
                                [Stiffness[f][1,1]],
                                [Stiffness[f][1,2]],
                                [Stiffness[f][2,0]],
                                [Stiffness[f][2,1]],
                                [Stiffness[f][2,2]],
                                [Stiffness[f][3,3]],
                                [Stiffness[f][4,4]],
                                [Stiffness[f][5,5]]])
    
    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    Lambda0 = np.exp(B[0,0]) - np.exp(B[2,0])
    Lambda0p = np.exp(B[1,0])
    Mu0 = np.exp(B[2,0]) / 2
    kl = B[3,0]
    km = B[4,0]
    l  = B[5,0]

    FName = Path(__file__).parents[1] / '_Results/ProposedModel_Trabecular.png'
    R2adj, NE = PlotFit(Stiffness, Lambda0, Lambda0p, Mu0, kl, km, l, Rho, eValues, FName=str(FName))


    # Print resulting values
    print(np.round([Lambda0, Lambda0p, Mu0]))
    print(np.round([kl, km],2))
    print(np.round(l,2))
    m1, m2, m3 = np.mean(eValues, axis=0)
    kl1 = kl*(m1+m1)/(2*m1*m1)
    kl2 = kl*(m2+m2)/(2*m2*m2)
    kl3 = kl*(m3+m3)/(2*m3*m3)
    kl12 = kl*(m1+m2)/(2*m1*m2)
    kl23 = kl*(m2+m3)/(2*m2*m3)
    kl31 = kl*(m3+m1)/(2*m3*m1)
    km23 = km*(m2+m3)/(2*m2*m3)
    km31 = km*(m3+m1)/(2*m3*m1)
    km12 = km*(m1+m2)/(2*m1*m2)
    print('   kl11   kl22   kl33   kl12   kl31   kl23   km12   km31   km23')
    print(np.round([kl1,kl2,kl3,kl12,kl31,kl23,km12,km31,km23],2))

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
