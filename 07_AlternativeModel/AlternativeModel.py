#%% !/usr/bin/env python3

Description = """
Description
"""

__author__ = ['Mathieu Simon']
__date_created__ = '14-03-2025'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def IsoTensor(Lambda0, Mu0):

    I = np.eye(3)
    I4 = Tensor.Dyadic(I,I)
    I4s = Tensor.Symmetric(I,I)
    S4 = Lambda0 * I4 + 2*Mu0 * I4s

    return S4

def IsoResiduals(Parameters, S):

    S4 = IsoTensor(Parameters[0], Parameters[1])
    S66 = Tensor.IsoMorphism3333_66(S4)
    SIso = np.zeros(S.shape)
    SIso[:,None] = S66
    return np.ravel(SIso - S)

def FitIso(S, FName=''):

    # Solve non linear system
    Results = least_squares(IsoResiduals,np.zeros(2), args=[S])
    Lambda0, Mu0 = Results.x

    # Build predicted tensor
    S4 = IsoTensor(Lambda0, Mu0)
    S66 = Tensor.IsoMorphism3333_66(S4)
    SPred = np.zeros(S.shape)
    SPred[:,None] = S66

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.concatenate([S[:,0,0],S[:,0,1],S[:,0,2],
                            S[:,1,0],S[:,1,1],S[:,1,2],
                            S[:,2,0],S[:,2,1],S[:,2,2],
                            S[:,3,3],S[:,4,4],S[:,5,5]])
    Y_Fit = np.concatenate([SPred[:,0,0],SPred[:,0,1],SPred[:,0,2],
                            SPred[:,1,0],SPred[:,1,1],SPred[:,1,2],
                            SPred[:,2,0],SPred[:,2,1],SPred[:,2,2],
                            SPred[:,3,3],SPred[:,4,4],SPred[:,5,5]])
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y_Obs) - len([Lambda0, Mu0])

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-len(Results.x)-1)

    NE = []
    for i, Si in enumerate(S):
        T_Obs = Si
        T_Fit = SPred[i]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), (Y_Fit).min()),
                       max(Y_Obs.max(), (Y_Fit).max()), len(Y_Obs))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(S[:,0,0], SPred[:,0,0], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,1,1], SPred[:,1,1], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,2,2], SPred[:,2,2], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,0,1], SPred[:,0,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,0,2], SPred[:,0,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,2], SPred[:,1,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,0], SPred[:,1,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,0], SPred[:,2,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,1], SPred[:,2,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,3,3], SPred[:,3,3], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,4,4], SPred[:,4,4], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,5,5], SPred[:,5,5], color=Colors[2], linestyle='none', marker='^')
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

    return [Lambda0, Mu0], R2adj, NE

def IsoRhoTensor(Lambda0, Mu0, k, Rho):

    I = np.eye(3)
    I4 = Tensor.Dyadic(I,I)
    I4s = Tensor.Symmetric(I,I)
    S4 = Lambda0 * Rho**k * I4 + 2*Mu0 * Rho**k * I4s

    return S4

def IsoRhoResiduals(Parameters, S, Rho):

    SIso = np.zeros(S.shape)
    for i, rho in enumerate(Rho):
        S4 = IsoRhoTensor(Parameters[0], Parameters[1], Parameters[2], rho)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SIso[i] = S66

    return np.ravel(SIso - S)

def FitIsoRho(S, Rho, FName=''):

    Results = least_squares(IsoRhoResiduals, np.zeros(3), args=[S, Rho])
    Lambda0, Mu0, k = Results.x

    # Build predicted tensors
    SPred = np.zeros(S.shape)
    for i, rho in enumerate(Rho):
        S4 = IsoRhoTensor(Lambda0, Mu0, k, rho)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SPred[i] = S66

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.concatenate([S[:,0,0],S[:,0,1],S[:,0,2],
                            S[:,1,0],S[:,1,1],S[:,1,2],
                            S[:,2,0],S[:,2,1],S[:,2,2],
                            S[:,3,3],S[:,4,4],S[:,5,5]])
    Y_Fit = np.concatenate([SPred[:,0,0],SPred[:,0,1],SPred[:,0,2],
                            SPred[:,1,0],SPred[:,1,1],SPred[:,1,2],
                            SPred[:,2,0],SPred[:,2,1],SPred[:,2,2],
                            SPred[:,3,3],SPred[:,4,4],SPred[:,5,5]])
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y_Obs) - len([Lambda0, Mu0])

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-len(Results.x)-1)

    NE = []
    for i, Si in enumerate(S):
        T_Obs = Si
        T_Fit = SPred[i]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), (Y_Fit).min()),
                       max(Y_Obs.max(), (Y_Fit).max()), len(Y_Obs))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(S[:,0,0], SPred[:,0,0], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,1,1], SPred[:,1,1], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,2,2], SPred[:,2,2], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,0,1], SPred[:,0,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,0,2], SPred[:,0,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,2], SPred[:,1,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,0], SPred[:,1,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,0], SPred[:,2,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,1], SPred[:,2,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,3,3], SPred[:,3,3], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,4,4], SPred[:,4,4], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,5,5], SPred[:,5,5], color=Colors[2], linestyle='none', marker='^')
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

    return [Lambda0, Mu0, k], R2adj, NE

def FabricTensor(Lambda0, Mu0, k, l, Rho, eValues, eVectors):

    m1, m2, m3 = eVectors
    M1, M2, M3 = np.outer(m1,m1), np.outer(m2,m2), np.outer(m3,m3)
    M = [M1, M2, M3]
    m1, m2, m3 = eValues
    m = [m1,m2,m3]

    S4 = np.zeros((3,3,3,3))
    for i in range(3):
        S4 += Lambda0 * Rho**k   * m[i]**l * m[i]**l   * Tensor.Dyadic(M[i],M[i])
        S4 += Lambda0 * Rho**k * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i],M[i-1])
        S4 += Lambda0 * Rho**k * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i-1],M[i])
        S4 +=   2*Mu0 * Rho**k   * m[i]**l * m[i]**l   * Tensor.Symmetric(M[i],M[i])
        S4 +=   2*Mu0 * Rho**k * m[i]**l * m[i-1]**l * Tensor.Symmetric(M[i-1],M[i])
    
    return S4

def FabricResiduals(Parameters, S, Rho, eValues, eVectors):

    Lambda0, Mu0, k, l = Parameters

    SFabricK = np.zeros(S.shape)
    for i, (rho, m, M) in enumerate(zip(Rho, eValues, eVectors)):
        S4 = FabricTensor(Lambda0, Mu0, k, l, rho, m, np.eye(3))
        # St = Tensor.TransformTensor(S4, np.eye(3), M)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SFabricK[i] = S66

    return np.ravel(SFabricK - S)

def FitFabric(S, Rho, eValues, eVectors, FName=''):

    # Estimate parameters
    X = np.matrix(np.zeros((len(Rho)*12, 5)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    m1 = eValues[:,0]
    m2 = eValues[:,1]
    m3 = eValues[:,2]
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[i]), np.log(m1[i] * m1[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m2[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m2[i])],
                                  [1, 0, 0, np.log(Rho[i]), np.log(m2[i] * m2[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m2[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m2[i] * m3[i])],
                                  [1, 0, 0, np.log(Rho[i]), np.log(m3[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m2[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m1[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m1[i] * m2[i])]])

        Y[Start:Stop] = np.log([[Si[0,0]],
                                [Si[0,1]],
                                [Si[0,2]],
                                [Si[1,0]],
                                [Si[1,1]],
                                [Si[1,2]],
                                [Si[2,0]],
                                [Si[2,1]],
                                [Si[2,2]],
                                [Si[3,3]],
                                [Si[4,4]],
                                [Si[5,5]]])
    
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    Estimate = [np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0]]
    Bounds = [0.0, np.inf]

    Results = least_squares(FabricResiduals, Estimate, bounds=Bounds, args=[S, Rho, eValues, eVectors])
    Lambda0, Mu0, k, l = Results.x

    # Build predicted tensors
    SPred = np.zeros(S.shape)
    for i, (rho, m, M) in enumerate(zip(Rho, eValues, eVectors)):
        S4 = FabricTensor(Lambda0, Mu0, k, l, rho, m, np.eye(3))
        # St = Tensor.TransformTensor(S4, np.eye(3), M)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SPred[i] = S66

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.concatenate([S[:,0,0],S[:,0,1],S[:,0,2],
                            S[:,1,0],S[:,1,1],S[:,1,2],
                            S[:,2,0],S[:,2,1],S[:,2,2],
                            S[:,3,3],S[:,4,4],S[:,5,5]])
    Y_Fit = np.concatenate([SPred[:,0,0],SPred[:,0,1],SPred[:,0,2],
                            SPred[:,1,0],SPred[:,1,1],SPred[:,1,2],
                            SPred[:,2,0],SPred[:,2,1],SPred[:,2,2],
                            SPred[:,3,3],SPred[:,4,4],SPred[:,5,5]])
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y_Obs) - len([Lambda0, Mu0])

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-len(Results.x)-1)

    NE = []
    for i, Si in enumerate(S):
        T_Obs = Si
        T_Fit = SPred[i]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), (Y_Fit).min()),
                       max(Y_Obs.max(), (Y_Fit).max()), len(Y_Obs))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(S[:,0,0], SPred[:,0,0], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,1,1], SPred[:,1,1], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,2,2], SPred[:,2,2], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,0,1], SPred[:,0,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,0,2], SPred[:,0,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,2], SPred[:,1,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,0], SPred[:,1,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,0], SPred[:,2,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,1], SPred[:,2,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,3,3], SPred[:,3,3], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,4,4], SPred[:,4,4], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,5,5], SPred[:,5,5], color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(Y)//12), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y)), xy=(0.3, 0.025), xycoords='axes fraction')
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

    return [Lambda0, Mu0, k, l], R2adj, NE

def FabricKTensor(Lambda0, Mu0, k1, k2, k3, l, Rho, eValues, eVectors):

    k = [k1, k2, k3]
    m1, m2, m3 = eVectors
    M1, M2, M3 = np.outer(m1,m1), np.outer(m2,m2), np.outer(m3,m3)
    M = [M1, M2, M3]
    m1, m2, m3 = eValues
    m = [m1,m2,m3]

    S4 = np.zeros((3,3,3,3))
    for i in range(3):
        S4 += Lambda0 * Rho**((k[i]+k[i])/2)   * m[i]**l * m[i]**l   * Tensor.Dyadic(M[i],M[i])
        S4 += Lambda0 * Rho**((k[i]+k[i-1])/2) * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i],M[i-1])
        S4 += Lambda0 * Rho**((k[i]+k[i-1])/2) * m[i]**l * m[i-1]**l * Tensor.Dyadic(M[i-1],M[i])
        S4 +=   2*Mu0 * Rho**((k[i]+k[i])/2)   * m[i]**l * m[i]**l   * Tensor.Symmetric(M[i],M[i])
        S4 +=   2*Mu0 * Rho**((k[i]+k[i-1])/2) * m[i]**l * m[i-1]**l * Tensor.Symmetric(M[i-1],M[i])
    
    return S4

def FabricKResiduals(Parameters, S, Rho, eValues, eVectors):

    Lambda0, Mu0, k1, k2, k3, l = Parameters

    SFabricK = np.zeros(S.shape)
    for i, (rho, m, M) in enumerate(zip(Rho, eValues, eVectors)):
        S4 = FabricKTensor(Lambda0, Mu0, k1, k2, k3, l, rho, m, np.eye(3))
        # St = Tensor.TransformTensor(S4, np.eye(3), M)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SFabricK[i] = S66

    return np.ravel(SFabricK - S)

def FitFabricK(S, Rho, eValues, eVectors, FName=''):

    # Estimate parameters
    X = np.matrix(np.zeros((len(Rho)*12, 7)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    m1 = eValues[:,0]
    m2 = eValues[:,1]
    m3 = eValues[:,2]
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0,   np.log(Rho[i]),             0,             0, np.log(m1[i] * m1[i])],
                                  [0, 1, 0, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])],
                                  [0, 1, 0, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])],
                                  [1, 0, 0,             0,   np.log(Rho[i]),             0, np.log(m2[i] * m2[i])],
                                  [0, 1, 0,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 1, 0,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [1, 0, 0,             0,             0,   np.log(Rho[i]), np.log(m3[i] * m3[i])],
                                  [0, 0, 1,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])]])

        Y[Start:Stop] = np.log([[Si[0,0]],
                                [Si[0,1]],
                                [Si[0,2]],
                                [Si[1,0]],
                                [Si[1,1]],
                                [Si[1,2]],
                                [Si[2,0]],
                                [Si[2,1]],
                                [Si[2,2]],
                                [Si[3,3]],
                                [Si[4,4]],
                                [Si[5,5]]])
    
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    Estimate = [np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0], B[5,0], B[6,0]]
    Bounds = [0.0, np.inf]

    Results = least_squares(FabricKResiduals, Estimate, bounds=Bounds, args=[S, Rho, eValues, eVectors])
    Lambda0, Mu0, k1, k2, k3, l = Results.x

    # Build predicted tensors
    SPred = np.zeros(S.shape)
    for i, (rho, m, M) in enumerate(zip(Rho, eValues, eVectors)):
        S4 = FabricKTensor(Lambda0, Mu0, k1, k2, k3, l, rho, m, np.eye(3))
        # St = Tensor.TransformTensor(S4, np.eye(3), M)
        S66 = Tensor.IsoMorphism3333_66(S4)
        SPred[i] = S66

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.concatenate([S[:,0,0],S[:,0,1],S[:,0,2],
                            S[:,1,0],S[:,1,1],S[:,1,2],
                            S[:,2,0],S[:,2,1],S[:,2,2],
                            S[:,3,3],S[:,4,4],S[:,5,5]])
    Y_Fit = np.concatenate([SPred[:,0,0],SPred[:,0,1],SPred[:,0,2],
                            SPred[:,1,0],SPred[:,1,1],SPred[:,1,2],
                            SPred[:,2,0],SPred[:,2,1],SPred[:,2,2],
                            SPred[:,3,3],SPred[:,4,4],SPred[:,5,5]])
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y_Obs) - len([Lambda0, Mu0])

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y_Obs)-1)/(len(Y_Obs)-len(Results.x)-1)

    NE = []
    for i, Si in enumerate(S):
        T_Obs = Si
        T_Fit = SPred[i]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), (Y_Fit).min()),
                       max(Y_Obs.max(), (Y_Fit).max()), len(Y_Obs))

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(S[:,0,0], SPred[:,0,0], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,1,1], SPred[:,1,1], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,2,2], SPred[:,2,2], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(S[:,0,1], SPred[:,0,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,0,2], SPred[:,0,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,2], SPred[:,1,2], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,1,0], SPred[:,1,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,0], SPred[:,2,0], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,2,1], SPred[:,2,1], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(S[:,3,3], SPred[:,3,3], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,4,4], SPred[:,4,4], color=Colors[2], linestyle='none', marker='^')
    Axes.plot(S[:,5,5], SPred[:,5,5], color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
    Axes.annotate(r'N ROIs   : ' + str(len(Y)//12), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y)), xy=(0.3, 0.025), xycoords='axes fraction')
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

    return [Lambda0, Mu0, k1, k2, k3, l], R2adj, NE



def FitRhoModel(X, Y, Parameters, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute observed and fitted values
    Lambda0 = Parameters.loc['Value','Lambda0']
    Mu0 = Parameters.loc['Value','Mu0']
    Y_Iso = Y.copy()
    Y_Iso[0::12] = np.log(Lambda0 + 2*Mu0)
    Y_Iso[1::12] = np.log(Lambda0)
    Y_Iso[2::12] = np.log(Lambda0)
    Y_Iso[3::12] = np.log(Lambda0)
    Y_Iso[4::12] = np.log(Lambda0 + 2*Mu0)
    Y_Iso[5::12] = np.log(Lambda0)
    Y_Iso[6::12] = np.log(Lambda0)
    Y_Iso[7::12] = np.log(Lambda0)
    Y_Iso[8::12] = np.log(Lambda0 + 2*Mu0)
    Y_Iso[9::12] = np.log(2*Mu0)
    Y_Iso[10::12] = np.log(2*Mu0)
    Y_Iso[11::12] = np.log(2*Mu0)
    Y = Y + Y_Iso

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.exp(Y)
    Y_Fit = np.exp(X * B + Y_Iso)
    Residuals = Y_Obs - Y_Fit
    DOFs = len(Y) - X.shape[1]
    Sigma = Residuals.T * Residuals / DOFs
    Cov = Sigma[0,0] * XTXi

    # Compute B confidence interval
    t_Alpha = t.interval(Alpha, DOFs)
    B_CI_Low = B.T + t_Alpha[0] * np.sqrt(np.diag(Cov))
    B_CI_Top = B.T + t_Alpha[1] * np.sqrt(np.diag(Cov))

    # Store parameters in data frame
    Parameters.loc['Value', 'k'] = B[0,0]
    Parameters.loc['95% CI Low', 'k'] = B_CI_Low[0,0]
    Parameters.loc['95% CI Top', 'k'] = B_CI_Top[0,0]

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y_Obs - Y_Obs.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y)-1)/(len(Y)-3-1)

    NE = []
    for i in range(0,len(Y),12):
        T_Obs = Y_Obs[i:i+12]
        T_Fit = Y_Fit[i:i+12]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(min(Y_Obs.min(), (Y_Fit).min()),
                       max(Y_Obs.max(), (Y_Fit).max()), len(Y))
    # B_0 = np.sort(np.sqrt(np.diag(X * Cov * X.T)))
    # CI_Line_u = np.exp(Line + t_Alpha[0] * B_0)
    # CI_Line_o = np.exp(Line + t_Alpha[1] * B_0)

    # Plots
    DPI = 500
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    # SMax = 3.5*1E4
    # SMin = 1.75*1E3
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = np.tile([True,False,False,False,True,False,False,False,True,False,False,False], len(Y)//12)
    Lambda_ij = np.tile([False,True,True,True,False,True,True,True,False,False,False,False], len(Y)//12)
    Mu_ij = np.tile([False,False,False,False,False,False,False,False,False,True,True,True], len(Y)//12)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
              color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
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
    # if len(FName) > 0:
    #     plt.savefig(FName)
    plt.show()

    return Parameters, R2adj, NE

def FitNewModel(X, Y, Alpha=0.95, FName=''):

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
    Parameters = pd.DataFrame(columns=['Lambda0', 'Lambda0p', 'Mu0', 'k1', 'k2', 'k3', 'l'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0], B[5,0], B[6,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Low[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4], B_CI_Low[0,5], B_CI_Low[0,6]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Top[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4], B_CI_Top[0,5], B_CI_Top[0,6]]

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
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = (X[:, 0] == 1)
    Lambda_ij = (X[:, 1] == 1)
    Mu_ij = (X[:, 2] == 1)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
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
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return Parameters, R2adj, NE

def FitZysset(X, Y, Alpha=0.95, FName=''):

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
    Parameters = pd.DataFrame(columns=['Lambda0', 'Lambda0', 'Mu0', 'k','l'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[1,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Top[0,1]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Low[0,1]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4]]

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
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    # SMax = 3.5*1E4
    # SMin = 1.75*1E3
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = (X[:, 0] == 1)
    Lambda_ij = (X[:, 1] == 1)
    Mu_ij = (X[:, 2] == 1)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
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

def FitkModel(X, Y, Alpha=0.95, FName=''):

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
    Parameters = pd.DataFrame(columns=['Lambda0', 'Lambda0p','Mu0','k12', 'k3'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Top[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Low[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4]]

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
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    # SMax = 3.5*1E4
    # SMin = 1.75*1E3
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = np.tile([True,False,False,False,True,False,False,False,True,False,False,False], len(Y)//12)
    Lambda_ij = np.tile([False,True,True,True,False,True,True,True,False,False,False,False], len(Y)//12)
    Mu_ij = np.tile([False,False,False,False,False,False,False,False,False,True,True,True], len(Y)//12)


    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
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

def FitRhoModel(X, Y, Alpha=0.95, FName=''):

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
    Parameters = pd.DataFrame(columns=['Lambda0', 'Lambda0p', 'Mu0', 'k11', 'k22', 'k33', 'k12', 'k23', 'k31', 'k44', 'k55', 'k66'])
    Parameters.loc['Value'] = [np.exp(B[0,0]) - 2*np.exp(B[2,0]), np.exp(B[1,0]), np.exp(B[2,0]), B[3,0], B[4,0], B[5,0], B[6,0], B[7,0], B[8,0], B[9,0], B[10,0], B[11,0]]
    Parameters.loc['95% CI Low'] = [np.exp(B_CI_Low[0,0]) - 2*np.exp(B_CI_Low[0,2]), np.exp(B_CI_Low[0,1]), np.exp(B_CI_Low[0,2]), B_CI_Low[0,3], B_CI_Low[0,4], B_CI_Low[0,5], B_CI_Low[0,6], B_CI_Low[0,7], B_CI_Low[0,8], B_CI_Low[0,9], B_CI_Low[0,10], B_CI_Low[0,11]]
    Parameters.loc['95% CI Top'] = [np.exp(B_CI_Top[0,0]) - 2*np.exp(B_CI_Top[0,2]), np.exp(B_CI_Top[0,1]), np.exp(B_CI_Top[0,2]), B_CI_Top[0,3], B_CI_Top[0,4], B_CI_Top[0,5], B_CI_Top[0,6], B_CI_Top[0,7], B_CI_Top[0,8], B_CI_Top[0,9], B_CI_Top[0,10], B_CI_Top[0,11]]

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
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    # SMax = 3.5*1E4
    # SMin = 1.75*1E3
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = (X[:, 0] == 1)
    Lambda_ij = (X[:, 1] == 1)
    Mu_ij = (X[:, 2] == 1)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
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

def OLS(X,Y):
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    return B


def FitIsotropicModel(X, Y, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute residuals, variance, and covariance matrix
    Y_Obs = np.array(Y)
    Y_Fit = np.array(X * B)
    Residuals = Y - X*B
    DOFs = len(Y) - X.shape[1]
    Sigma = Residuals.T * Residuals / DOFs
    Cov = Sigma[0,0] * XTXi

    # Compute B confidence interval
    t_Alpha = t.interval(Alpha, DOFs)
    B_CI_Low = B.T + t_Alpha[0] * np.sqrt(np.diag(Cov))
    B_CI_Top = B.T + t_Alpha[1] * np.sqrt(np.diag(Cov))

    # Store parameters in data frame
    Parameters = pd.DataFrame(columns=['Lambda0','Mu0'])
    Parameters.loc['Value'] = [B[1,0], B[2,0]]
    Parameters.loc['95% CI Low'] = [B_CI_Low[0,1], B_CI_Low[0,2]]
    Parameters.loc['95% CI Top'] = [B_CI_Top[0,1], B_CI_Top[0,2]]

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
    # SMax = max([Y_Obs.max(), Y_Fit.max()]) * 1.2
    # SMin = min([Y_Obs.min(), Y_Fit.min()]) / 1.2
    # SMax = 3.5*1E4
    # SMin = 1.75*1E3
    Colors=[(0,0,1),(0,1,0),(1,0,0)]
    Lambda_ii = np.tile([True,False,False,False,True,False,False,False,True,False,False,False], len(Y)//12)
    Lambda_ij = np.tile([False,True,True,True,False,True,True,True,False,False,False,False], len(Y)//12)
    Mu_ij = np.tile([False,False,False,False,False,False,False,False,False,True,True,True], len(Y)//12)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(Y_Obs[Lambda_ii], Y_Fit[Lambda_ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(Y_Obs[Lambda_ij], Y_Fit[Lambda_ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(Y_Obs[Mu_ij], Y_Fit[Mu_ij],
              color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0, 0, 0), linestyle='--')
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

    # Define paths to fabric and homogenisation data of cortical ROIs
    FabPath =Path(__file__).parents[1] / '05_Homogenization/Fabric'
    ElaPath =Path(__file__).parents[1] / '05_Homogenization/Elasticity'

    # List folders
    Folders = [F.name for F in FabPath.iterdir() if F.is_dir()]

    # Define ROIs numbering
    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    # Read fabric data
    CortVal = np.zeros((len(Folders), 16, 3))
    CortVec = np.zeros((len(Folders), 16, 3, 3))
    CortRho = np.zeros((len(Folders), 16))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get homogenization stress results
            Fabric = Read.Fabric(FabPath / Folder / (ROI+'.fab'))
            CortVal[f,r] = Fabric[0]
            CortVec[f,r] = Fabric[1]
            CortRho[f,r] = Fabric[2]

    # Read homogenisation results
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    CortStiff = np.zeros((len(Folders), 16, 6, 6))
    CortTraStiff = np.zeros((len(Folders), 16, 6, 6))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get homogenization stress results for isotropic matrix
            File = open(ElaPath / Folder / (ROI + f'_Isotropic.out'), 'r').readlines()
            Stress = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stress[i,j] = float(File[i+4].split()[j+1])

            # Compute stiffness
            Stiffness = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stiffness[i,j] = Stress[i,j] / Strain[i]

            # Symmetrize matrix
            Stiffness = 1/2 * (Stiffness + Stiffness.T)

            # Project onto orthotropy and store
            CortStiff[f,r] = Tensor.OrthoProj(Stiffness)

            # Get homogenization stress results for transverse isotropic matrix
            File = open(ElaPath / Folder / (ROI + f'_Transverse.out'), 'r').readlines()
            Stress = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stress[i,j] = float(File[i+4].split()[j+1])

            # Compute stiffness
            Stiffness = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stiffness[i,j] = Stress[i,j] / Strain[i]

            # Symmetrize matrix
            Stiffness = 1/2 * (Stiffness + Stiffness.T)

            # Project onto orthotropy and store
            CortTraStiff[f,r] = Tensor.OrthoProj(Stiffness)

    # Reduce dimensionality of cortical data
    CortVal = CortVal.reshape(-1,3)
    CortVec = CortVec.reshape(-1,3,3)
    CortRho = CortRho.reshape(-1)
    CortStiff = CortStiff.reshape(-1,6,6)
    CortTraStiff = CortTraStiff.reshape(-1,6,6)

    # Compute corresponding compliance 
    CortComp = np.zeros(CortStiff.shape)
    CortTraComp = np.zeros(CortTraStiff.shape)
    for i, (S, St) in enumerate(zip(CortStiff, CortTraStiff)):
        CortComp[i] = np.linalg.inv(S)
        CortTraComp[i] = np.linalg.inv(St)


    # Plot fabric as function of porosity
    MA, RA = [], []
    Delta = 0.05
    for i in range(int(np.ceil((1-CortRho.min())/Delta))):
        Min = (1-CortRho) > (i)*Delta - 0.05
        Max = (1-CortRho) < (i+1)*Delta + 0.05
        RA.append((i+0.5)*Delta)
        MA.append(np.mean(CortVal[Min & Max], axis=0))

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-CortRho, CortVal[:,2] / CortVal[:,0], marker='o',
            linestyle='none', fillstyle='none', color=(1,0,0))
    Axis.plot(RA, np.array(MA)[:,2] / np.array(MA)[:,0], color=(0,0,1))
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('DA')
    Axis.set_ylim([1.0,6])
    plt.show(Figure)

    # Plot stiffness as function of porosity
    Labels = ['Stiffness', 'Modulus']
    for i, (Iso, Tra) in enumerate(zip([CortStiff, CortComp], [CortTraStiff, CortTraComp])):

        if i == 0:

            Figure, Axis = plt.subplots(1,1)
            Axis.plot(1-CortRho, Iso[:,0,0], label=r'$e_{11}$', marker='o',
                    linestyle='none', fillstyle='none', color=(1,0,0))
            Axis.plot(1-CortRho, Iso[:,1,1], label=r'$e_{22}$', marker='o',
                        linestyle='none', fillstyle='none', color=(1,0,1))
            Axis.plot(1-CortRho, Iso[:,2,2], label=r'$e_{33}$', marker='o',
                        linestyle='none', fillstyle='none', color=(0,0,1))
            Axis.set_xlabel(r'1 - $\rho$ (-)')
            Axis.set_ylabel(Labels[i] + ' (MPa)')
            Axis.legend(loc='upper right')
            plt.show(Figure)

            Figure, Axis = plt.subplots(1,1)
            Axis.plot(1-CortRho, Tra[:,0,0], label=r'$e_{11}$', marker='o',
                    linestyle='none', fillstyle='none', color=(1,0,0))
            Axis.plot(1-CortRho, Tra[:,1,1], label=r'$e_{22}$', marker='o',
                        linestyle='none', fillstyle='none', color=(1,0,1))
            Axis.plot(1-CortRho, Tra[:,2,2], label=r'$e_{33}$', marker='o',
                        linestyle='none', fillstyle='none', color=(0,0,1))
            Axis.set_xlabel(r'1 - $\rho$ (-)')
            Axis.set_ylabel(Labels[i] + ' (MPa)')
            Axis.legend(loc='upper right')
            plt.show(Figure)
        
        else:

            Figure, Axis = plt.subplots(1,1)
            Axis.plot(1-CortRho, 1/Iso[:,0,0], label=r'$e_{11}$', marker='o',
                    linestyle='none', fillstyle='none', color=(1,0,0))
            Axis.plot(1-CortRho, 1/Iso[:,1,1], label=r'$e_{22}$', marker='o',
                        linestyle='none', fillstyle='none', color=(1,0,1))
            Axis.plot(1-CortRho, 1/Iso[:,2,2], label=r'$e_{33}$', marker='o',
                        linestyle='none', fillstyle='none', color=(0,0,1))
            Axis.set_xlabel(r'1 - $\rho$ (-)')
            Axis.set_ylabel(Labels[i] + ' (MPa)')
            Axis.legend(loc='upper right')
            plt.show(Figure)

            Figure, Axis = plt.subplots(1,1)
            Axis.plot(1-CortRho, 1/Tra[:,0,0], label=r'$e_{11}$', marker='o',
                    linestyle='none', fillstyle='none', color=(1,0,0))
            Axis.plot(1-CortRho, 1/Tra[:,1,1], label=r'$e_{22}$', marker='o',
                        linestyle='none', fillstyle='none', color=(1,0,1))
            Axis.plot(1-CortRho, 1/Tra[:,2,2], label=r'$e_{33}$', marker='o',
                        linestyle='none', fillstyle='none', color=(0,0,1))
            Axis.set_xlabel(r'1 - $\rho$ (-)')
            Axis.set_ylabel(Labels[i] + ' (MPa)')
            Axis.legend(loc='upper right')
            plt.show(Figure)

        # Compare stiffness ratio
        Figure, Axis = plt.subplots(1,1)
        Axis.plot([0,1],[0,1], linestyle='--', color=(0,0,0), label='1:1', linewidth=1)
        Axis.plot(Iso[:,0,0] / Iso[:,0,0].max(), Tra[:,0,0] / Tra[:,0,0].max(), label=r'$e_{11}$', marker='o',
                linestyle='none', fillstyle='none', color=(1,0,0))
        Axis.plot(Iso[:,1,1] / Iso[:,1,1].max(), Tra[:,1,1] / Tra[:,1,1].max(), label=r'$e_{22}$', marker='o',
                    linestyle='none', fillstyle='none', color=(1,0,1))
        Axis.plot(Iso[:,2,2] / Iso[:,2,2].max(), Tra[:,2,2] / Tra[:,2,2].max(), label=r'$e_{33}$', marker='o',
                    linestyle='none', fillstyle='none', color=(0,0,1))
        Axis.set_xlabel('Isotropic Matrix ' + Labels[i] + ' Ratio (-)')
        Axis.set_ylabel('Transverse Isotropic Matrix ' + Labels[i] + ' Ratio (-)')
        Axis.legend()
        plt.show(Figure)


    # Fit homogenization with simple isotropic model (2 constants)
    S = CortStiff
    FName = str(Path(__file__).parent / 'RegressionIso_Isotropic.png')
    ParametersIso, R2adjIso, NEIso = FitIso(S, FName)

    # Define isotropic model scaled with rho
    S = CortStiff
    Rho = CortRho
    FName = str(Path(__file__).parent / 'RegressionRho_Isotropic.png')
    ParametersRho, R2adjRho, NERho = FitIsoRho(S, Rho, FName)

    # Define fabric-based model scaled with rho^k and m^l 
    S = CortStiff
    Rho = CortRho
    eValues = CortVal
    eVectors = CortVec
    FName = str(Path(__file__).parent / 'RegressionFabric_Isotropic.png')
    ParametersFabric, R2adjFabric, NEFabric = FitFabric(S, Rho, eValues, eVectors)

    # Define fabric-based model scaled with rho and k according to direction
    S = CortStiff
    Rho = CortRho
    eValues = CortVal
    eVectors = CortVec
    FName = str(Path(__file__).parent / 'RegressionFabricK_Isotropic.png')
    ParametersFabricK, R2adjFabricK, NEFabricK = FitFabricK(S, Rho, eValues, eVectors)

    # Build theoric stiffness and compliance tensors
    Lambda0 = ParametersFabric[0]
    Mu0 = ParametersFabric[1]
    k = ParametersFabric[2]
    l = ParametersFabric[3]
    N = 100
    Rho = np.linspace(CortRho.min(), CortRho.max(), N)
    Values = np.mean(CortVal, axis=0)
    Vectors = np.eye(3)
    FabStiffness = np.zeros((N, 6, 6))
    FabCompliance = np.zeros((N, 6, 6))
    for i, rho in enumerate(Rho):
        T = FabricTensor(Lambda0, Mu0, k, l, rho, Values, Vectors)
        T = Tensor.IsoMorphism3333_66(T)
        T = Tensor.OrthoProj(T)
        FabStiffness[i] = T
        FabCompliance[i] = np.linalg.inv(T)

    S = CortStiff
    C = CortComp

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-CortRho, S[:,0,0], label=r'$e_{11}$', color=(1,0,0), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabStiffness[:,0,0], label=r'Fit $e_{11}$', color=(1,0,0))
    Axis.plot(1-CortRho, S[:,1,1], label=r'$e_{22}$', color=(1,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabStiffness[:,1,1], label=r'Fit $e_{22}$', color=(1,0,1))
    Axis.plot(1-CortRho, S[:,2,2], label=r'$e_{33}}$', color=(0,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabStiffness[:,2,2], label=r'Fit $e_{33}$', color=(0,0,1))
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Stiffness (MPa)')
    Axis.legend(loc='upper right')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-CortRho, 1/C[:,0,0], label=r'$e_{11}$', color=(1,0,0), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabCompliance[:,0,0], label=r'Fit $e_{11}$', color=(1,0,0))
    Axis.plot(1-CortRho, 1/C[:,1,1], label=r'$e_{22}$', color=(1,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabCompliance[:,1,1], label=r'Fit $e_{22}$', color=(1,0,1))
    Axis.plot(1-CortRho, 1/C[:,2,2], label=r'$e_{33}$', color=(0,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabCompliance[:,2,2], label=r'Fit $e_{33}$', color=(0,0,1))
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Modulus (MPa)')
    Axis.legend(loc='upper right')
    plt.show(Figure)
    

    Lambda0 = ParametersFabricK[0]
    Mu0 = ParametersFabricK[1]
    k1, k2, k3 = ParametersFabricK[2:5]
    l = ParametersFabricK[5]
    N = 100
    Rho = np.linspace(CortRho.min(), CortRho.max(), N)
    Values = np.mean(CortVal, axis=0)
    Vectors = np.eye(3)
    FabKStiffness = np.zeros((N, 6, 6))
    FabKCompliance = np.zeros((N, 6, 6))
    for i, rho in enumerate(Rho):
        T = FabricKTensor(Lambda0, Mu0, k1, k2, k3, l, rho, Values, Vectors)
        T = Tensor.IsoMorphism3333_66(T)
        T = Tensor.OrthoProj(T)
        FabKStiffness[i] = T
        FabKCompliance[i] = np.linalg.inv(T)

    S = CortStiff
    C = CortComp

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-CortRho, S[:,0,0], label=r'$e_{11}$', color=(1,0,0), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabKStiffness[:,0,0], label=r'Fit $e_{11}$', color=(1,0,0))
    Axis.plot(1-CortRho, S[:,1,1], label=r'$e_{22}$', color=(1,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabKStiffness[:,1,1], label=r'Fit $e_{22}$', color=(1,0,1))
    Axis.plot(1-CortRho, S[:,2,2], label=r'$e_{33}$', color=(0,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, FabKStiffness[:,2,2], label=r'Fit $e_{33}$', color=(0,0,1))
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Stiffness (MPa)')
    Axis.legend(loc='upper right')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(1-CortRho, 1/C[:,0,0], label=r'$e_{11}$', color=(1,0,0), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabKCompliance[:,0,0], label=r'Fit $e_{11}$', color=(1,0,0))
    Axis.plot(1-CortRho, 1/C[:,1,1], label=r'$e_{22}$', color=(1,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabKCompliance[:,1,1], label=r'Fit $e_{22}$', color=(1,0,1))
    Axis.plot(1-CortRho, 1/C[:,2,2], label=r'$e_{33}$', color=(0,0,1), 
              linestyle='none', marker='o', fillstyle='none')
    Axis.plot(1-Rho, 1/FabKCompliance[:,2,2], label=r'Fit $e_{33}$', color=(0,0,1))
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Modulus (MPa)')
    Axis.legend(loc='upper right')
    plt.show(Figure)



    # Build stiffness and compliance tensors
    Lambda0 = ParametersFabric[0]
    Mu0 = ParametersFabric[1]
    k = ParametersFabric[2]
    l = ParametersFabric[3]
    FabStiffness = np.zeros((len(CortRho), 6, 6))
    FabCompliance = np.zeros((len(CortRho), 6, 6))
    for i, (Rho, Values, Vectors) in enumerate(zip(CortRho, CortVal, CortVec)):
        T = FabricTensor(Lambda0, Mu0, k, l, Rho, Values, Vectors)
        T = Tensor.IsoMorphism3333_66(T)
        T = Tensor.OrthoProj(T)
        FabStiffness[i] = T
        FabCompliance[i] = np.linalg.inv(T)
    

    Lambda0 = ParametersFabricK[0]
    Mu0 = ParametersFabricK[1]
    k1, k2, k3 = ParametersFabricK[2:5]
    l = ParametersFabricK[5]
    FabKStiffness = np.zeros((len(CortRho), 6, 6))
    FabKCompliance = np.zeros((len(CortRho), 6, 6))
    for i, (Rho, Values, Vectors) in enumerate(zip(CortRho, CortVal, CortVec)):
        T = FabricKTensor(Lambda0, Mu0, k1, k2, k3, l, Rho, Values, Vectors)
        T = Tensor.IsoMorphism3333_66(T)
        T = Tensor.OrthoProj(T)
        FabKStiffness[i] = T
        FabKCompliance[i] = np.linalg.inv(T)

    # Analyze anisotropy
    S = CortStiff
    X = np.concatenate([[np.ones(len(S))],[1-CortRho]], axis=0).T
    FName = str(Path(__file__).parent / 'Anisotropy_SIsotropic.png')
    Colors = [(0,0,1),(1,0,0),(0,0,0)]
    Labels = [r'Standard Model $S_{33}/S_{11}$', r'Proposed Model $S_{33}/S_{11}$',r'Simulation $S_{33}/S_{11}$']

    Figure, Axis = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=200)
    # Y = np.array([m[2] / m[0] for m in CortVal])
    # B = OLS(np.matrix(X), np.matrix(Y).T)
    # xLine = 1-np.linspace(CortRho.max(), CortRho.min(), 10)
    # yLine = B[0,0] + B[1,0] * xLine
    # Axis.plot(X[:,1], Y, label=r'Fabric $m_{3}/m_{1}$', linestyle='none',
    #             marker='o', color=(0.5,0.5,0.5), fillstyle='none')
    # Axis.plot(xLine, yLine, linestyle='--', color=(0.5,0.5,0.5), linewidth=1)
    # Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x',
    #                 (0.05,0.65), color=(0.5,0.5,0.5), xycoords='axes fraction')
    # Axis.plot(1-CortRho, [v[2] / v[0] for v in CortVal], linestyle='none',
    #           marker='o', color=(0.5,0.5,0.5), fillstyle='none')
    for i, Stiffness in enumerate([FabStiffness, FabKStiffness, S]):
            Y = np.array([s[2,2] / s[0,0] for s in Stiffness])
            B = OLS(np.matrix(X), np.matrix(Y).T)
            xLine = 1-np.linspace(CortRho.max(), CortRho.min(), 10)
            yLine = B[0,0] + B[1,0] * xLine
            Axis.plot(X[:,1], Y, label=Labels[i], linestyle='none',
                      marker='o', color=Colors[i], fillstyle='none')
            Axis.plot(xLine, yLine, linestyle='--', color=Colors[i], linewidth=1)
            Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x',
                          (0.05,0.65-0.075*(i+1)), color=Colors[i], xycoords='axes fraction')
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Anisotropy (-)')
    Axis.legend(loc='upper left')
    plt.savefig(fname=FName)
    plt.show(Figure)

    # Plot compliance anisotropy
    E = CortComp
    X = np.concatenate([[np.ones(len(E))],[1-CortRho]], axis=0).T
    FName = str(Path(__file__).parent / 'Anisotropy_EIsotropic.png')
    Colors = [(0,0,1),(1,0,0),(0,0,0)]
    Labels = [r'Standard Model $E_{33}/E_{11}$', r'Proposed Model $E_{33}/E_{11}$',r'Simulation $E_{33}/E_{11}$']

    Figure, Axis = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=200)
    # Y = np.array([m[2] / m[0] for m in CortVal])
    # B = OLS(np.matrix(X), np.matrix(Y).T)
    # xLine = 1-np.linspace(CortRho.max(), CortRho.min(), 10)
    # yLine = B[0,0] + B[1,0] * xLine
    # Axis.plot(X[:,1], Y, label=r'Fabric $m_{3}/m_{1}$', linestyle='none',
    #             marker='o', color=(0.5,0.5,0.5), fillstyle='none')
    # Axis.plot(xLine, yLine, linestyle='--', color=(0.5,0.5,0.5), linewidth=1)
    # Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x',
    #                 (0.05,0.65), color=(0.5,0.5,0.5), xycoords='axes fraction')
    # Axis.plot(1-CortRho, [v[2] / v[0] for v in CortVal], linestyle='none',
    #           marker='o', color=(0.5,0.5,0.5), fillstyle='none')
    for i, Compliance in enumerate([FabCompliance, FabKCompliance, E]):
            Y = np.array([c[0,0] / c[2,2] for c in Compliance])
            B = OLS(np.matrix(X), np.matrix(Y).T)
            xLine = 1-np.linspace(CortRho.max(), CortRho.min(), 10)
            yLine = B[0,0] + B[1,0] * xLine
            Axis.plot(X[:,1], Y, label=Labels[i], linestyle='none',
                      marker='o', color=Colors[i], fillstyle='none')
            Axis.plot(xLine, yLine, linestyle='--', color=Colors[i], linewidth=1)
            Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x',
                          (0.05,0.80-0.07*(i+1)), color=Colors[i], xycoords='axes fraction')
    Axis.set_xlabel(r'1 - $\rho$ (-)')
    Axis.set_ylabel('Anisotropy (-)')
    Axis.legend(loc='upper left')
    plt.savefig(fname=FName)
    plt.show(Figure)




    





    # Fit homogenization with Zysset-Curnier theorical model
    Rho = CortRho
    m1 = CortVal[:,0]
    m2 = CortVal[:,1]
    m3 = CortVal[:,2]
    S = CortStiff

    X = np.matrix(np.zeros((len(Rho)*12, 5)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[i]), np.log(m1[i] ** 2)],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m2[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m2[i] * m1[i])],
                                  [1, 0, 0, np.log(Rho[i]), np.log(m2[i] ** 2)],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m2[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m3[i] * m1[i])],
                                  [0, 1, 0, np.log(Rho[i]), np.log(m3[i] * m2[i])],
                                  [1, 0, 0, np.log(Rho[i]), np.log(m3[i] ** 2)],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m2[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m3[i] * m1[i])],
                                  [0, 0, 1, np.log(Rho[i]), np.log(m1[i] * m2[i])]])
        
        Y[Start:Stop] = np.log([[Si[0,0]],
                         [Si[0,1]],
                         [Si[0,2]],
                         [Si[1,0]],
                         [Si[1,1]],
                         [Si[1,2]],
                         [Si[2,0]],
                         [Si[2,1]],
                         [Si[2,2]],
                         [Si[3,3]],
                         [Si[4,4]],
                         [Si[5,5]]])
    
    FName = Path(__file__).parent / 'RegressionZysset.png'
    Parameters, R2adj, NE = FitZysset(X, Y, FName=str(FName))

    # Fit homogenization with alternative theorical model
    Rho = CortRho
    m1 = CortVal[:,0]
    m2 = CortVal[:,1]
    m3 = CortVal[:,2]
    S = CortTraStiff

    X = np.matrix(np.zeros((len(Rho)*12, 5)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [1, 0, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [0, 1, 0, np.log(Rho[i]), 0],
                                  [1, 0, 0, 0, np.log(Rho[i])],
                                  [0, 0, 1, 0, np.log(Rho[i])],
                                  [0, 0, 1, np.log(Rho[i]), 0],
                                  [0, 0, 1, np.log(Rho[i]), 0]])
        
        Y[Start:Stop] = np.log([[Si[0,0]],
                         [Si[0,1]],
                         [Si[0,2]],
                         [Si[1,0]],
                         [Si[1,1]],
                         [Si[1,2]],
                         [Si[2,0]],
                         [Si[2,1]],
                         [Si[2,2]],
                         [Si[3,3]],
                         [Si[4,4]],
                         [Si[5,5]]])
    
    FName = Path(__file__).parent / 'RegressionNew.png'
    Parameters, R2adj, NE = FitkModel(X, Y, FName=str(FName))
    L0 = Parameters.loc['Value', 'Lambda0']
    L0p = Parameters.loc['Value', 'Lambda0p']
    Mu0 = Parameters.loc['Value', 'Mu0']
    k12 = Parameters.loc['Value', 'k12']
    k3 = Parameters.loc['Value', 'k3']






    # Isotropic tensor
    S = CortStiff
    X = np.matrix(np.zeros((len(S)*12, 3)))
    Y = np.matrix(np.zeros((len(S)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [1, 0, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [1, 0, 0],
                                  [0, 0, 1],
                                  [0, 0, 1],
                                  [0, 0, 1]])
        
        Y[Start:Stop] = np.array([[Si[0,0]],
                                  [Si[0,1]],
                                  [Si[0,2]],
                                  [Si[1,0]],
                                  [Si[1,1]],
                                  [Si[1,2]],
                                  [Si[2,0]],
                                  [Si[2,1]],
                                  [Si[2,2]],
                                  [Si[3,3]],
                                  [Si[4,4]],
                                  [Si[5,5]]])
    
    FName = Path(__file__).parent / 'RegressionIso.png'
    Parameters, R2adj, NE = FitIsotropicModel(X, Y, FName=str(FName))

    # Fit homogenization with isotropic model scaled with rho
    X = np.matrix(np.zeros((len(CortRho)*12, 1)))
    Y = np.matrix(np.zeros((len(CortRho)*12, 1)))
    Rho = CortRho
    m1 = CortVal[:,0]
    m2 = CortVal[:,1]
    m3 = CortVal[:,2]
    k = [2.9, 3.6, 1.8]
    S = CortStiff
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.log([[Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]],
                                [Rho[i]]])
        
        Y[Start:Stop] = np.array([[np.log(Si[0,0]) - np.log(Lambda0 + 2*Mu0)],
                                  [np.log(Si[0,1]) - np.log(Lambda0)],
                                  [np.log(Si[0,2]) - np.log(Lambda0)],
                                  [np.log(Si[1,0]) - np.log(Lambda0)],
                                  [np.log(Si[1,1]) - np.log(Lambda0 + 2*Mu0)],
                                  [np.log(Si[1,2]) - np.log(Lambda0)],
                                  [np.log(Si[2,0]) - np.log(Lambda0)],
                                  [np.log(Si[2,1]) - np.log(Lambda0)],
                                  [np.log(Si[2,2]) - np.log(Lambda0 + 2*Mu0)],
                                  [np.log(Si[3,3]) - np.log(2*Mu0)],
                                  [np.log(Si[4,4]) - np.log(2*Mu0)],
                                  [np.log(Si[5,5]) - np.log(2*Mu0)]])
    
    FName = Path(__file__).parent / 'RegressionRho.png'
    ParametersRho, R2adj, NE = FitRhoModel(X, Y, Parameters, FName=str(FName))

    L = np.array([[m[0]**(2*l) * Rho**k[0] * (Lambda0 + 2*Mu0), m[0]**l * m[1]**l * Rho**((k[0]+k[1])/2) * Lambda0, m[0]**l * m[2]**l * Rho**((k[0]+k[2])/2) * Lambda0, 0, 0, 0],
                  [m[0]**l * m[1]**l * Rho**((k[0]+k[1])/2) * Lambda0, m[1]**(2*l) * Rho**k[1] * (Lambda0 + 2*Mu0), m[1]**l * m[2]**l * Rho**((k[1]+k[2])/2) * Lambda0, 0, 0, 0],
                  [m[0]**l * m[2]**l * Rho**((k[0]+k[2])/2) * Lambda0, m[1]**l * m[2]**l * Rho**((k[1]+k[2])/2) * Lambda0, m[2]**(2*l) * Rho**k[2] * (Lambda0 + 2*Mu0), 0, 0, 0],
                  [0, 0, 0, 2 * m[1]**l *m[2]**l * Rho**((k[1]+k[2])/2) * Mu0, 0, 0],
                  [0, 0, 0, 0, 2 * m[0]**l *m[2]**l * Rho**((k[0]+k[2])/2) * Mu0, 0],
                  [0, 0, 0, 0, 0, 2 * m[1]**l *m[0]**l * Rho**((k[1]+k[0])/2) * Mu0]])
    
    # Build linear regression
    X = np.matrix([[1, 0, 0,   np.log(Rho),             0,             0, np.log(m[0] * m[0])],
                   [0, 1, 0, np.log(Rho)/2, np.log(Rho)/2,             0, np.log(m[0] * m[1])],
                   [0, 1, 0, np.log(Rho)/2,             0, np.log(Rho)/2, np.log(m[0] * m[2])],
                   [0, 1, 0, np.log(Rho)/2, np.log(Rho)/2,             0, np.log(m[0] * m[1])],
                   [1, 0, 0,             0,   np.log(Rho),             0, np.log(m[1] * m[1])],
                   [0, 1, 0,             0, np.log(Rho)/2, np.log(Rho)/2, np.log(m[1] * m[2])],
                   [0, 1, 0, np.log(Rho)/2,             0, np.log(Rho)/2, np.log(m[0] * m[2])],
                   [0, 1, 0,             0, np.log(Rho)/2, np.log(Rho)/2, np.log(m[1] * m[2])],
                   [1, 0, 0,             0,             0,   np.log(Rho), np.log(m[2] * m[2])],
                   [0, 0, 1,             0, np.log(Rho)/2, np.log(Rho)/2, np.log(m[1] * m[2])],
                   [0, 0, 1, np.log(Rho)/2,             0, np.log(Rho)/2, np.log(m[0] * m[2])],
                   [0, 0, 1, np.log(Rho)/2, np.log(Rho)/2,             0, np.log(m[0] * m[1])]])
    B = np.matrix([[np.log(Lambda0 + 2*Mu0)],
                   [np.log(Lambda0)],
                   [np.log(2*Mu0)],
                   [k[0]],
                   [k[1]],
                   [k[2]],
                   [l]])
    Y = np.exp(X*B)
    S = np.array([[Y[0,0], Y[1,0], Y[2,0],      0,       0,      0],
                  [Y[3,0], Y[4,0], Y[5,0],      0,       0,      0],
                  [Y[6,0], Y[7,0], Y[8,0],      0,       0,      0],
                  [     0,      0,      0, Y[9,0],       0,      0],
                  [     0,      0,      0,      0, Y[10,0],      0],
                  [     0,      0,      0,      0,       0, Y[11,0]]])
    print(S)

    # Fit homogenization with new theorical model
    Rho = CortRho
    m1 = CortVal[:,0]
    m2 = CortVal[:,1]
    m3 = CortVal[:,2]
    S = CortStiff

    X = np.matrix(np.zeros((len(Rho)*12, 7)))
    Y = np.matrix(np.zeros((len(Rho)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0,   np.log(Rho[i]),             0,             0, np.log(m1[i] * m1[i])],
                                  [0, 1, 0, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])],
                                  [0, 1, 0, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])],
                                  [1, 0, 0,             0,   np.log(Rho[i]),             0, np.log(m2[i] * m2[i])],
                                  [0, 1, 0,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [0, 1, 0, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 1, 0,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [1, 0, 0,             0,             0,   np.log(Rho[i]), np.log(m3[i] * m3[i])],
                                  [0, 0, 1,             0, np.log(Rho[i])/2, np.log(Rho[i])/2, np.log(m2[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i])/2,             0, np.log(Rho[i])/2, np.log(m1[i] * m3[i])],
                                  [0, 0, 1, np.log(Rho[i])/2, np.log(Rho[i])/2,             0, np.log(m1[i] * m2[i])]])

        Y[Start:Stop] = np.log([[Si[0,0]],
                                [Si[0,1]],
                                [Si[0,2]],
                                [Si[1,0]],
                                [Si[1,1]],
                                [Si[1,2]],
                                [Si[2,0]],
                                [Si[2,1]],
                                [Si[2,2]],
                                [Si[3,3]],
                                [Si[4,4]],
                                [Si[5,5]]])
    
    FName = str(Path(__file__).parent / 'RegressionFabricK_Isotropic3.png')
    Parameters, R2adj, NE = FitNewModel(X, Y, FName=FName)



    # Read trabecular ROIs data
    TrabPath = Path(__file__).parents[1] / '06_PorosityEffect'
    Data = pd.read_csv(TrabPath/'ROIsData.csv')

    # Keep ROIs with CV < 0.263
    F = Data['Variation Coefficient'] < 0.263
    Files = []
    for I, Row in Data[F].iterrows():
        N = Row['ROI Number']
        S = Row['Scan']
        Files.append(f'{N}_{S}')

    # Define paths to fabric and homogenisation data
    FabPath = TrabPath / 'Fabric'
    ElaPath = TrabPath / 'Elasticity'

    # Read fabric data
    TrabVal = np.zeros((len(Files),3))
    TrabVec = np.zeros((len(Files),3,3))
    TrabRho = np.zeros((len(Files)))
    for i, F in enumerate(Files):
        Fabric = Read.Fabric(FabPath / (F+'.fab'))
        TrabVal[i] = Fabric[0]
        TrabVec[i] = Fabric[1]
        TrabRho[i] = Fabric[2]

    # Read homogenisation results
    TrabStiff = np.zeros((len(Files),6,6))
    for i, F in enumerate(Files):

        # Read stiffness matrix
        ComplianceMatrix = Read.ComplianceDat(ElaPath / (F+'.dat'))
        ComplianceMatrix = (ComplianceMatrix+ComplianceMatrix.T)/2
        StiffnessMatrix = np.linalg.inv(ComplianceMatrix)
        StiffnessMatrix = (StiffnessMatrix+StiffnessMatrix.T)/2

        # Transform stiffness into fabric coordinate system
        Mandel = Tensor.Engineering2MandelNotation(StiffnessMatrix)
        S4 = Tensor.IsoMorphism66_3333(Mandel)
        St = Tensor.TransformTensor(S4, np.eye(3), TrabVec[i])
        StiffnessMatrix = Tensor.IsoMorphism3333_66(St)
        Se = Tensor.Mandel2EngineeringNotation(StiffnessMatrix)
        TrabStiff[i] = 1/2 * (Se+Se.T)


    # Investigate the effect of porosity on the stiffness tensor
    Figure, Axis = plt.subplots(1,3, dpi=192, sharex=True, figsize=(12,4))
    for i in range(6):
        for j in range(6):
            if i < 3 and j < 3:
                if i == j:
                    if i == 0:
                        Cort = (0,1,1)
                        Trab = (0.6,0,1)
                        El = r'$\lambda_{11}$'
                    elif i == 1:
                        Cort = (0,0.6,1)
                        Trab = (1,0,1)
                        El = r'$\lambda_{22}$'
                    else:
                        Cort = (0,0,1)
                        Trab = (1,0,0)
                        El = r'$\lambda_{33}$'

                    Axis[0].plot(1-CortRho, CortStiff[:,i,j]/1E3,
                                 color=Cort, marker='o', linestyle='none', fillstyle='none')
                    Axis[0].plot(1-TrabRho, TrabStiff[:,i,j]/1E3,
                                 color=Trab, marker='o', linestyle='none', fillstyle='none')
                    Axis[0].plot([], label=El + ' Cortical',
                                 color=Cort, marker='o', linestyle='none', fillstyle='none')
                    Axis[0].plot([], label=El + ' Trabecular',
                                 color=Trab, marker='o', linestyle='none', fillstyle='none')
                else:
                    if (i == 0 and j == 2) or (i == 2 and j == 0):
                        Cort = (0,1,1)
                        Trab = (0.6,0,1)
                        El = r'$\lambda_{13}$'
                    elif (i == 1 and j == 2) or (i == 2 and j == 1):
                        Cort = (0,0.6,1)
                        Trab = (1,0,1)
                        El = r'$\lambda_{23}$'
                    else:
                        Cort = (0,0,1)
                        Trab = (1,0,0)
                        El = r'$\lambda_{12}$'
                    Axis[1].plot(1-CortRho, CortStiff[:,i,j]/1E3,
                                 color=Cort, marker='o', linestyle='none', fillstyle='none')
                    Axis[1].plot(1-TrabRho, TrabStiff[:,i,j]/1E3,
                                 color=Trab, marker='o', linestyle='none', fillstyle='none')
                    if (i == 0 and j == 1) or (i == 1 and j == 2) or (i == 2 and j == 0):
                        Axis[1].plot([], label=El + ' Cortical',
                                    color=Cort, marker='o', linestyle='none', fillstyle='none')
                        Axis[1].plot([], label=El + ' Trabecular',
                                    color=Trab, marker='o', linestyle='none', fillstyle='none')
            elif i == j:
                if i == 3:
                    Cort = (0,1,1)
                    Trab = (0.6,0,1)
                    El = r'$\mu_{23}$'
                elif i == 4:
                    Cort = (0,0.6,1)
                    Trab = (1,0,1)
                    El = r'$\mu_{31}$'
                else:
                    Cort = (0,0,1)
                    Trab = (1,0,0)
                    El = r'$\mu_{12}$'
                Axis[2].plot(1-CortRho, CortStiff[:,i,j]/1E3,
                                color=Cort, marker='o', linestyle='none', fillstyle='none')
                Axis[2].plot(1-TrabRho, TrabStiff[:,i,j]/1E3,
                                color=Trab, marker='o', linestyle='none', fillstyle='none')
                Axis[2].plot([], label=El + ' Cortical',
                                 color=Cort, marker='o', linestyle='none', fillstyle='none')
                Axis[2].plot([], label=El + ' Trabecular',
                                 color=Trab, marker='o', linestyle='none', fillstyle='none')
            else:
                continue
    
    Axis[0].legend(loc='upper right',  fontsize=9)
    Axis[1].legend(loc='upper right',  fontsize=9)
    Axis[2].legend(loc='upper right',  fontsize=9)
    Axis[1].set_xlabel(r'1-$\rho$ (-)')
    Axis[0].set_ylabel('Stiffness (GPa)')
    plt.show(Figure)

    # Fit individual curves to each stiffness component
    def Curve(Rho, a, b):
        return a*Rho**b
    BVTV = np.hstack((CortRho, TrabRho))
    Args = np.argsort(BVTV)
    BVTV = BVTV[Args]
    S = np.vstack((CortStiff, TrabStiff))[Args]
    P11 = curve_fit(Curve, BVTV, S[:,0,0])
    P12 = curve_fit(Curve, BVTV, S[:,0,1])
    P13 = curve_fit(Curve, BVTV, S[:,0,2])
    P22 = curve_fit(Curve, BVTV, S[:,1,1])
    P23 = curve_fit(Curve, BVTV, S[:,1,2])
    P33 = curve_fit(Curve, BVTV, S[:,2,2])
    P44 = curve_fit(Curve, BVTV, S[:,3,3])
    P55 = curve_fit(Curve, BVTV, S[:,4,4])
    P66 = curve_fit(Curve, BVTV, S[:,5,5])

    # Plot individual curves with data
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,0,0], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P11[0]), color=(1,0,0))
    Axis.plot(BVTV, S[:,1,1], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P22[0]), color=(0,1,0))
    Axis.plot(BVTV, S[:,2,2], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P33[0]), color=(0,0,1))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,0,1], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P12[0]), color=(1,0,0))
    Axis.plot(BVTV, S[:,0,2], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P13[0]), color=(0,1,0))
    Axis.plot(BVTV, S[:,1,2], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P23[0]), color=(0,0,1))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,3,3], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P44[0]), color=(1,0,0))
    Axis.plot(BVTV, S[:,4,4], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P55[0]), color=(0,1,0))
    Axis.plot(BVTV, S[:,5,5], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Curve(BVTV, *P66[0]), color=(0,0,1))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    # Define general parameters
    Lambda0 = np.mean([P11[0][0], P22[0][0], P33[0][0]])
    Lambda0p = np.mean([P12[0][0], P23[0][0], P13[0][0]])
    Mu0 = np.mean([P44[0][0], P55[0][0], P66[0][0]])
    k12 = np.mean([P11[0][1], P22[0][1]])
    k3 = P33[0][1]
    l = np.mean([P12[0][1], P23[0][1], P13[0][1]])
    m = np.mean([P44[0][1], P55[0][1], P66[0][1]])

    # Plot with average constants
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,0,0], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0*BVTV**k12, color=(1,0,0))
    Axis.plot(BVTV, S[:,1,1], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0*BVTV**k12, color=(0,1,0))
    Axis.plot(BVTV, S[:,2,2], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0*BVTV**k3, color=(0,0,1))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,0,1], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0p*BVTV**l, color=(1,0,0))
    Axis.plot(BVTV, S[:,0,2], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0p*BVTV**l, color=(1,0,0))
    Axis.plot(BVTV, S[:,1,2], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Lambda0p*BVTV**l, color=(1,0,0))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(BVTV, S[:,3,3], color=(1,0,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Mu0*BVTV**m, color=(1,0,0))
    Axis.plot(BVTV, S[:,4,4], color=(0,1,0,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Mu0*BVTV**m, color=(1,0,0))
    Axis.plot(BVTV, S[:,5,5], color=(0,0,1,0.2), linestyle='none', marker='o', fillstyle='none')
    Axis.plot(BVTV, Mu0*BVTV**m, color=(1,0,0))
    plt.xlabel(r'$\rho$')
    plt.ylabel('Stiffness (MPa)')
    plt.show(Figure)

    
    # Fit homogenization with Zysset-Curnier theorical model
    BVTV = np.hstack((CortRho, TrabRho))
    m1 = np.hstack((CortVal[:,0], TrabVal[:,0]))
    m2 = np.hstack((CortVal[:,1], TrabVal[:,1]))
    m3 = np.hstack((CortVal[:,2], TrabVal[:,2]))
    S = np.vstack((CortStiff, TrabStiff))

    X = np.matrix(np.zeros((len(BVTV)*12, 5)))
    Y = np.matrix(np.zeros((len(BVTV)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(BVTV[i]), np.log(m1[i] ** 2)],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m1[i] * m2[i])],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m1[i] * m3[i])],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m2[i] * m1[i])],
                                  [1, 0, 0, np.log(BVTV[i]), np.log(m2[i] ** 2)],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m2[i] * m3[i])],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m3[i] * m1[i])],
                                  [0, 1, 0, np.log(BVTV[i]), np.log(m3[i] * m2[i])],
                                  [1, 0, 0, np.log(BVTV[i]), np.log(m3[i] ** 2)],
                                  [0, 0, 1, np.log(BVTV[i]), np.log(m2[i] * m3[i])],
                                  [0, 0, 1, np.log(BVTV[i]), np.log(m3[i] * m1[i])],
                                  [0, 0, 1, np.log(BVTV[i]), np.log(m1[i] * m2[i])]])
        
        Y[Start:Stop] = np.log([[Si[0,0]],
                         [Si[0,1]],
                         [Si[0,2]],
                         [Si[1,0]],
                         [Si[1,1]],
                         [Si[1,2]],
                         [Si[2,0]],
                         [Si[2,1]],
                         [Si[2,2]],
                         [Si[3,3]],
                         [Si[4,4]],
                         [Si[5,5]]])
    
    FName = Path(__file__).parent / 'RegressionZysset.png'
    Parameters, R2adj, NE = FitZysset(X, Y, FName=str(FName))

    # Fit homogenization with alternative theorical model
    BVTV = np.hstack((CortRho, TrabRho))
    S = np.vstack((CortStiff, TrabStiff))

    X = np.matrix(np.zeros((len(BVTV)*12, 12)))
    Y = np.matrix(np.zeros((len(BVTV)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0],
                                  [1, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0],
                                  [1, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0, 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0, 0],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, np.log(BVTV[i]), 0],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, np.log(BVTV[i])]])
        
        Y[Start:Stop] = np.log([[Si[0,0]],
                                [Si[0,1]],
                                [Si[0,2]],
                                [Si[1,0]],
                                [Si[1,1]],
                                [Si[1,2]],
                                [Si[2,0]],
                                [Si[2,1]],
                                [Si[2,2]],
                                [Si[3,3]],
                                [Si[4,4]],
                                [Si[5,5]]])
    
    FName = Path(__file__).parent / 'RegressionNewFull.png'
    Parameters, R2adj, NE = FitRhoModel(X, Y, FName=str(FName))

    # Fit homogenization with alternative theorical model
    BVTV = CortRho #np.hstack((CortRho, TrabRho))
    m1 = CortVal[:,0] #np.hstack((CortVal[:,0], TrabVal[:,0]))
    m2 = CortVal[:,1] #np.hstack((CortVal[:,1], TrabVal[:,1]))
    m3 = CortVal[:,2] #np.hstack((CortVal[:,2], TrabVal[:,2]))
    S =  CortStiff #np.vstack((CortStiff, TrabStiff))

    X = np.matrix(np.zeros((len(BVTV)*12, 5)))
    Y = np.matrix(np.zeros((len(BVTV)*12, 1)))
    M = np.matrix(np.zeros((len(BVTV)*12, 1)))
    for i, Si in enumerate(S):
        
        Start, Stop = 12*i, 12*(i+1)

        # Build system
        X[Start:Stop] = np.array([[1, 0, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [1, 0, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [0, 1, 0, np.log(BVTV[i]), 0],
                                  [1, 0, 0, 0, np.log(BVTV[i])],
                                  [0, 0, 1, 0, np.log(BVTV[i])],
                                  [0, 0, 1, np.log(BVTV[i]), 0],
                                  [0, 0, 1, np.log(BVTV[i]), 0]])
        
        Y[Start:Stop] = np.log([[Si[0,0]],
                         [Si[0,1]],
                         [Si[0,2]],
                         [Si[1,0]],
                         [Si[1,1]],
                         [Si[1,2]],
                         [Si[2,0]],
                         [Si[2,1]],
                         [Si[2,2]],
                         [Si[3,3]],
                         [Si[4,4]],
                         [Si[5,5]]])
        
        M[Start:Stop] = np.log([[m1[i]],
                                 [1],
                                 [1],
                                 [1],
                                 [m2[i]],
                                 [1],
                                 [1],
                                 [1],
                                 [m3[i]],
                                 [1],
                                 [1],
                                 [1]]) * (1-BVTV[i])
    
    FName = Path(__file__).parent / 'RegressionNew.png'
    Parameters, R2adj, NE = FitkModel(X, Y, FName=str(FName))
    L0 = Parameters.loc['Value', 'Lambda0']
    L0p = Parameters.loc['Value', 'Lambda0p']
    Mu0 = Parameters.loc['Value', 'Mu0']
    k12 = Parameters.loc['Value', 'k12']
    k3 = Parameters.loc['Value', 'k3']

    # Model parameters
    BVTV = np.hstack((CortRho, TrabRho))
    e1Val = np.hstack((CortVal[:,0], TrabVal[:,0]))
    e2Val = np.hstack((CortVal[:,1], TrabVal[:,1]))
    e3Val = np.hstack((CortVal[:,2], TrabVal[:,2]))
    S = np.vstack((CortStiff, TrabStiff))

    # Build model predicted values
    MStiff = np.zeros((len(BVTV),6,6))
    MStiff[:,0,0] = (L0+2*Mu0) * BVTV**k12
    MStiff[:,2,2] = (L0+2*Mu0) * BVTV**k3

    # Add architecture anisotropy
    m1 = e1Val * (1-BVTV)
    m2 = e2Val * (1-BVTV)
    m3 = e3Val * (1-BVTV)
    MArchStiff = np.zeros((len(BVTV),6,6))
    MArchStiff[:,0,0] = (L0+2*Mu0) * BVTV**k12 * m1
    MArchStiff[:,2,2] = (L0+2*Mu0) * BVTV**k3 * m3

    # Plot anisotropy vs BVTV
    Figure, Axis = plt.subplots(1,1, dpi=192)
    Colors = [(1,0,0), (0,0,1), (0,0,0)]
    Labels = ['Simulation $S_{33} / S_{11}$', 'Model $S_{33} / S_{11}$', 'Model+Fabric $S_{33} / S_{11}$']
    X = np.ones((len(BVTV),2), float)
    X[:,1] = 1-BVTV
    for i, V in enumerate([S, MStiff, MArchStiff]):
            Y = np.array([v[2,2] / v[0,0] for v in V])
            B = OLS(np.matrix(X), np.matrix(Y).T)
            # xLine = 1-np.linspace(BVTV.min(), BVTV.max(), 10)
            # yLine = B[0,0] + B[1,0] * xLine
            Axis.plot(X[:,1], Y, label=Labels[i], linestyle='none', marker='o', color=Colors[i])
            # Axis.plot(xLine, yLine, linestyle='--', color=Colors[i], linewidth=1)
            # Axis.annotate(f'y = {round(B[1,0],2)}x + {round(B[0,0],2)}',(0.05,2.3-0.1*i), color=Colors[i])
    Axis.set_xlabel(r'1-$\rho$ (-)')
    # Axis.set_ylim([0.95, 2.85])
    Axis.set_ylabel('Degree of Anisotropy (-)')
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.15))
    # plt.savefig(Path(__file__).parent / 'Plots/AnisotropyBVTV.png')
    plt.show(Figure)

    # Read homogenisation results with tranvserse isotropic material
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    CortStiffT = np.zeros((len(Folders), 16, 6, 6))
    for f, Folder in enumerate(Folders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get homogenization stress results
            File = open(ElaPath / Folder / (ROI + f'_Transverse.out'), 'r').readlines()
            Stress = np.zeros((6,6))
            for i in range(6):
                for j in range(6):
                    Stress[i,j] = float(File[i+4].split()[j+1])

            # Compute stiffness
            for i in range(6):
                for j in range(6):
                    CortStiffT[f,r,i,j] = Stress[i,j] / Strain[i]

        # Transform stiffness into fabric coordinate system
        StiffnessMatrix = 1/2 * (CortStiffT[f,r] + CortStiffT[f,r].T)
        Mandel = Tensor.Engineering2MandelNotation(StiffnessMatrix)
        S4 = Tensor.IsoMorphism66_3333(Mandel)
        St = Tensor.TransformTensor(S4, np.eye(3), CortVec[f+16*r])
        StiffnessMatrix = Tensor.IsoMorphism3333_66(St)
        Se = Tensor.Mandel2EngineeringNotation(StiffnessMatrix)
        CortStiffT[f,r] = 1/2 * (Se+Se.T)

    CortStiffT = CortStiffT.reshape(-1,6,6)

    # Add tissue anisotropy
    m1t, m2t, m3t = 0.841, 0.841, 1.318 # Virtual tissue fabric
    m1f = CortVal[:,0]
    m2f = CortVal[:,1]
    m3f = CortVal[:,2]
    m1 = m1t*CortRho + m1f*(1-CortRho)
    m2 = m2t*CortRho + m2f*(1-CortRho)
    m3 = m3t*CortRho + m3f*(1-CortRho)
    m1 = CortVal[:,0] ** ((1-CortRho)/2) * m1t ** CortRho
    m2 = CortVal[:,1] ** ((1-CortRho)/2) * m2t ** CortRho
    m3 = CortVal[:,2] ** ((1-CortRho)/2) * m3t ** CortRho

    # Build model predicted values for cortical bone
    MStiff = np.zeros((len(CortRho),6,6))
    MStiff[:,0,0] = (L0+2*Mu0) * CortRho**k12 * m1
    MStiff[:,2,2] = (L0+2*Mu0) * CortRho**k3 * m3

    # Plot anisotropy vs BVTV
    Figure, Axis = plt.subplots(1,1, dpi=192)
    Colors = [(1,0,0), (0,0,1), (1,0,1)]
    Labels = ['Simulation $S_{33} / S_{11}$', 'Model $S_{33} / S_{11}$']
    X = np.ones((len(CortRho),2), float)
    X[:,1] = 1-CortRho
    for i, V in enumerate([CortStiffT, MStiff]):
            Y = np.array([v[2,2] / v[0,0] for v in V])
            B = OLS(np.matrix(X), np.matrix(Y).T)
            xLine = 1-np.linspace(CortRho.min(), CortRho.max(), 10)
            yLine = B[0,0] + B[1,0] * xLine
            Axis.plot(X[:,1], Y, label=Labels[i], linestyle='none', marker='o', color=Colors[i])
            Axis.plot(xLine, yLine, linestyle='--', color=Colors[i], linewidth=1)
            Axis.annotate(f'y = {round(B[1,0],2)}x + {round(B[0,0],2)}',(0.05,2.3-0.1*i), color=Colors[i])
    Axis.set_xlabel(r'1-$\rho$ (-)')
    # Axis.set_ylim([0.95, 2.85])
    Axis.set_ylabel('Degree of Anisotropy (-)')
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.15))
    # plt.savefig(Path(__file__).parent / 'Plots/AnisotropyBVTV.png')
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
