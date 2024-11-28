#%% !/usr/bin/env python3

Description = """
Read ISQ files and plot them in 3D using pyvista
"""

__author__ = ['Mathieu Simon']
__date_created__ = '12-11-2024'
__date__ = '15-11-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats.distributions import norm, t

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

#%% Functions

def DyadicProduct(A,B):

    if A.size == 3:
        C = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                C[i,j] = A[i]*B[j]

    elif A.size == 9:
        C = np.zeros((3,3,3,3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                            C[i,j,k,l] = A[i,j] * B[k,l]

    else:
        print('Matrices sizes mismatch')

    return C

def SymmetricProduct(A,B):

    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i,j,k,l] = (1/2)*(A[i,k]*B[j,l]+A[i,l]*B[j,k])

    return C

def ComplianceTensor(E1, E2, E3, Mu23, Mu31, Mu12, Nu12, Nu13, Nu23, EigenVectors=np.eye(3)):

    # Define constants
    Mu32, Mu13, Mu21 = Mu23, Mu31, Mu12
    Nu21 = Nu12 * E2 / E1
    Nu31 = Nu13 * E3 / E1
    Nu32 = Nu23 * E3 / E2

    # Group into list for loop computation
    E = [E1, E2, E3]
    Nu = np.array([[Nu13, Nu12], [Nu21, Nu23], [Nu32, Nu31]])
    Mu = np.array([[Mu13, Mu12], [Mu21, Mu23], [Mu32, Mu31]])

    # Build compliance tensor
    ComplianceTensor = np.zeros((3, 3, 3, 3))
    for i in range(3):
        Mi = DyadicProduct(EigenVectors[i], EigenVectors[i])
        Part1 = 1 / E[i] * DyadicProduct(Mi, Mi)
        ComplianceTensor += Part1

        for ii in range(3 - 1):
            j = i - ii - 1
            Mj = DyadicProduct(EigenVectors[j], EigenVectors[j])
            Part2 = -Nu[i, ii] / E[i] * DyadicProduct(Mi, Mj)
            Part3 = 1 / (2 * Mu[i, ii]) * SymmetricProduct(Mi, Mj)
            ComplianceTensor += Part2 + Part3

    ComplianceTensor = IsoMorphism3333_66(ComplianceTensor)

    return ComplianceTensor

def Engineering2MandelNotation(A):

    B = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            if i < 3 and j >= 3:
                B[i,j] = A[i,j] * np.sqrt(2)
            elif i >= 3 and j < 3:
                B[i,j] = A[i,j] * np.sqrt(2)
            elif i >= 3 and j >= 3:
                B[i,j] = A[i,j] * 2
            else:
                B[i, j] = A[i, j]

    return B

def IsoMorphism66_3333(A):

    # Check symmetry
    Symmetry = True
    for i in range(6):
        for j in range(6):
            if not A[i,j] == A[j,i]:
                Symmetry = False
                break
    if Symmetry == False:
        print('Matrix is not symmetric!')
        return

    B = np.zeros((3,3,3,3))

    # Build 4th tensor
    B[0, 0, 0, 0] = A[0, 0]
    B[1, 1, 0, 0] = A[1, 0]
    B[2, 2, 0, 0] = A[2, 0]
    B[1, 2, 0, 0] = A[3, 0] / np.sqrt(2)
    B[2, 0, 0, 0] = A[4, 0] / np.sqrt(2)
    B[0, 1, 0, 0] = A[5, 0] / np.sqrt(2)

    B[0, 0, 1, 1] = A[0, 1]
    B[1, 1, 1, 1] = A[1, 1]
    B[2, 2, 1, 1] = A[2, 1]
    B[1, 2, 1, 1] = A[3, 1] / np.sqrt(2)
    B[2, 0, 2, 1] = A[4, 1] / np.sqrt(2)
    B[0, 1, 2, 1] = A[5, 1] / np.sqrt(2)

    B[0, 0, 2, 2] = A[0, 2]
    B[1, 1, 2, 2] = A[1, 2]
    B[2, 2, 2, 2] = A[2, 2]
    B[1, 2, 2, 2] = A[3, 2] / np.sqrt(2)
    B[2, 0, 2, 2] = A[4, 2] / np.sqrt(2)
    B[0, 1, 2, 2] = A[5, 2] / np.sqrt(2)

    B[0, 0, 1, 2] = A[0, 3] / np.sqrt(2)
    B[1, 1, 1, 2] = A[1, 3] / np.sqrt(2)
    B[2, 2, 1, 2] = A[2, 3] / np.sqrt(2)
    B[1, 2, 1, 2] = A[3, 3] / 2
    B[2, 0, 1, 2] = A[4, 3] / 2
    B[0, 1, 1, 2] = A[5, 3] / 2

    B[0, 0, 2, 0] = A[0, 4] / np.sqrt(2)
    B[1, 1, 2, 0] = A[1, 4] / np.sqrt(2)
    B[2, 2, 2, 0] = A[2, 4] / np.sqrt(2)
    B[1, 2, 2, 0] = A[3, 4] / 2
    B[2, 0, 2, 0] = A[4, 4] / 2
    B[0, 1, 2, 0] = A[5, 4] / 2

    B[0, 0, 0, 1] = A[0, 5] / np.sqrt(2)
    B[1, 1, 0, 1] = A[1, 5] / np.sqrt(2)
    B[2, 2, 0, 1] = A[2, 5] / np.sqrt(2)
    B[1, 2, 0, 1] = A[3, 5] / 2
    B[2, 0, 0, 1] = A[4, 5] / 2
    B[0, 1, 0, 1] = A[5, 5] / 2



    # Add minor symmetries ijkl = ijlk and ijkl = jikl

    B[0, 0, 0, 0] = B[0, 0, 0, 0]
    B[0, 0, 0, 0] = B[0, 0, 0, 0]

    B[0, 0, 1, 0] = B[0, 0, 0, 1]
    B[0, 0, 0, 1] = B[0, 0, 0, 1]

    B[0, 0, 1, 1] = B[0, 0, 1, 1]
    B[0, 0, 1, 1] = B[0, 0, 1, 1]

    B[0, 0, 2, 1] = B[0, 0, 1, 2]
    B[0, 0, 1, 2] = B[0, 0, 1, 2]

    B[0, 0, 2, 2] = B[0, 0, 2, 2]
    B[0, 0, 2, 2] = B[0, 0, 2, 2]

    B[0, 0, 0, 2] = B[0, 0, 2, 0]
    B[0, 0, 2, 0] = B[0, 0, 2, 0]



    B[0, 1, 0, 0] = B[0, 1, 0, 0]
    B[1, 0, 0, 0] = B[0, 1, 0, 0]

    B[0, 1, 1, 0] = B[0, 1, 0, 1]
    B[1, 0, 0, 1] = B[0, 1, 0, 1]

    B[0, 1, 1, 1] = B[0, 1, 1, 1]
    B[1, 0, 1, 1] = B[0, 1, 1, 1]

    B[0, 1, 2, 1] = B[0, 1, 1, 2]
    B[1, 0, 1, 2] = B[0, 1, 1, 2]

    B[0, 1, 2, 2] = B[0, 1, 2, 2]
    B[1, 0, 2, 2] = B[0, 1, 2, 2]

    B[0, 1, 0, 2] = B[0, 1, 2, 0]
    B[1, 0, 2, 0] = B[0, 1, 2, 0]



    B[1, 1, 0, 0] = B[1, 1, 0, 0]
    B[1, 1, 0, 0] = B[1, 1, 0, 0]

    B[1, 1, 1, 0] = B[1, 1, 0, 1]
    B[1, 1, 0, 1] = B[1, 1, 0, 1]

    B[1, 1, 1, 1] = B[1, 1, 1, 1]
    B[1, 1, 1, 1] = B[1, 1, 1, 1]

    B[1, 1, 2, 1] = B[1, 1, 1, 2]
    B[1, 1, 1, 2] = B[1, 1, 1, 2]

    B[1, 1, 2, 2] = B[1, 1, 2, 2]
    B[1, 1, 2, 2] = B[1, 1, 2, 2]

    B[1, 1, 0, 2] = B[1, 1, 2, 0]
    B[1, 1, 2, 0] = B[1, 1, 2, 0]



    B[1, 2, 0, 0] = B[1, 2, 0, 0]
    B[2, 1, 0, 0] = B[1, 2, 0, 0]

    B[1, 2, 1, 0] = B[1, 2, 0, 1]
    B[2, 1, 0, 1] = B[1, 2, 0, 1]

    B[1, 2, 1, 1] = B[1, 2, 1, 1]
    B[2, 1, 1, 1] = B[1, 2, 1, 1]

    B[1, 2, 2, 1] = B[1, 2, 1, 2]
    B[2, 1, 1, 2] = B[1, 2, 1, 2]

    B[1, 2, 2, 2] = B[1, 2, 2, 2]
    B[2, 1, 2, 2] = B[1, 2, 2, 2]

    B[1, 2, 0, 2] = B[1, 2, 2, 0]
    B[2, 1, 2, 0] = B[1, 2, 2, 0]



    B[2, 2, 0, 0] = B[2, 2, 0, 0]
    B[2, 2, 0, 0] = B[2, 2, 0, 0]

    B[2, 2, 1, 0] = B[2, 2, 0, 1]
    B[2, 2, 0, 1] = B[2, 2, 0, 1]

    B[2, 2, 1, 1] = B[2, 2, 1, 1]
    B[2, 2, 1, 1] = B[2, 2, 1, 1]

    B[2, 2, 2, 1] = B[2, 2, 1, 2]
    B[2, 2, 1, 2] = B[2, 2, 1, 2]

    B[2, 2, 2, 2] = B[2, 2, 2, 2]
    B[2, 2, 2, 2] = B[2, 2, 2, 2]

    B[2, 2, 0, 2] = B[2, 2, 2, 0]
    B[2, 2, 2, 0] = B[2, 2, 2, 0]



    B[2, 0, 0, 0] = B[2, 0, 0, 0]
    B[0, 2, 0, 0] = B[2, 0, 0, 0]

    B[2, 0, 1, 0] = B[2, 0, 0, 1]
    B[0, 2, 0, 1] = B[2, 0, 0, 1]

    B[2, 0, 1, 1] = B[2, 0, 1, 1]
    B[0, 2, 1, 1] = B[2, 0, 1, 1]

    B[2, 0, 2, 1] = B[2, 0, 1, 2]
    B[0, 2, 1, 2] = B[2, 0, 1, 2]

    B[2, 0, 2, 2] = B[2, 0, 2, 2]
    B[0, 2, 2, 2] = B[2, 0, 2, 2]

    B[2, 0, 0, 2] = B[2, 0, 2, 0]
    B[0, 2, 2, 0] = B[2, 0, 2, 0]


    # Complete minor symmetries
    B[0, 2, 1, 0] = B[0, 2, 0, 1]
    B[0, 2, 0, 2] = B[0, 2, 2, 0]
    B[0, 2, 2, 1] = B[0, 2, 1, 2]

    B[1, 0, 1, 0] = B[1, 0, 0, 1]
    B[1, 0, 0, 2] = B[1, 0, 2, 0]
    B[1, 0, 2, 1] = B[1, 0, 1, 2]

    B[2, 1, 1, 0] = B[2, 1, 0, 1]
    B[2, 1, 0, 2] = B[2, 1, 2, 0]
    B[2, 1, 2, 1] = B[2, 1, 1, 2]


    # Add major symmetries ijkl = klij
    B[0, 1, 1, 1] = B[1, 1, 0, 1]
    B[1, 0, 1, 1] = B[1, 1, 1, 0]

    B[0, 2, 1, 1] = B[1, 1, 0, 2]
    B[2, 0, 1, 1] = B[1, 1, 2, 0]


    return B

def CheckMinorSymmetry(A):
    MinorSymmetry = True
    for i in range(3):
        for j in range(3):
            PartialTensor = A[:,:, i, j]
            if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                MinorSymmetry = True
            else:
                MinorSymmetry = False
                break

    if MinorSymmetry == True:
        for i in range(3):
            for j in range(3):
                PartialTensor = np.squeeze(A[i, j,:,:])
                if PartialTensor[1, 0] == PartialTensor[0, 1] and PartialTensor[2, 0] == PartialTensor[0, 2] and PartialTensor[1, 2] == PartialTensor[2, 1]:
                    MinorSymmetry = True
                else:
                    MinorSymmetry = False
                    break

    return MinorSymmetry

def IsoMorphism3333_66(A):

    if CheckMinorSymmetry == False:
        print('Tensor does not present minor symmetry')
    else:

        B = np.zeros((6,6))

        B[0, 0] = A[0, 0, 0, 0]
        B[0, 1] = A[0, 0, 1, 1]
        B[0, 2] = A[0, 0, 2, 2]
        B[0, 3] = np.sqrt(2) * A[0, 0, 1, 2]
        B[0, 4] = np.sqrt(2) * A[0, 0, 2, 0]
        B[0, 5] = np.sqrt(2) * A[0, 0, 0, 1]

        B[1, 0] = A[1, 1, 0, 0]
        B[1, 1] = A[1, 1, 1, 1]
        B[1, 2] = A[1, 1, 2, 2]
        B[1, 3] = np.sqrt(2) * A[1, 1, 1, 2]
        B[1, 4] = np.sqrt(2) * A[1, 1, 2, 0]
        B[1, 5] = np.sqrt(2) * A[1, 1, 0, 1]

        B[2, 0] = A[2, 2, 0, 0]
        B[2, 1] = A[2, 2, 1, 1]
        B[2, 2] = A[2, 2, 2, 2]
        B[2, 3] = np.sqrt(2) * A[2, 2, 1, 2]
        B[2, 4] = np.sqrt(2) * A[2, 2, 2, 0]
        B[2, 5] = np.sqrt(2) * A[2, 2, 0, 1]

        B[3, 0] = np.sqrt(2) * A[1, 2, 0, 0]
        B[3, 1] = np.sqrt(2) * A[1, 2, 1, 1]
        B[3, 2] = np.sqrt(2) * A[1, 2, 2, 2]
        B[3, 3] = 2 * A[1, 2, 1, 2]
        B[3, 4] = 2 * A[1, 2, 2, 0]
        B[3, 5] = 2 * A[1, 2, 0, 1]

        B[4, 0] = np.sqrt(2) * A[2, 0, 0, 0]
        B[4, 1] = np.sqrt(2) * A[2, 0, 1, 1]
        B[4, 2] = np.sqrt(2) * A[2, 0, 2, 2]
        B[4, 3] = 2 * A[2, 0, 1, 2]
        B[4, 4] = 2 * A[2, 0, 2, 0]
        B[4, 5] = 2 * A[2, 0, 0, 1]

        B[5, 0] = np.sqrt(2) * A[0, 1, 0, 0]
        B[5, 1] = np.sqrt(2) * A[0, 1, 1, 1]
        B[5, 2] = np.sqrt(2) * A[0, 1, 2, 2]
        B[5, 3] = 2 * A[0, 1, 1, 2]
        B[5, 4] = 2 * A[0, 1, 2, 0]
        B[5, 5] = 2 * A[0, 1, 0, 1]

        return B
    
def TransformTensor(A,OriginalBasis,NewBasis):

    # Build change of coordinate matrix
    O = OriginalBasis
    N = NewBasis

    Q = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            Q[i,j] = np.dot(O[i,:],N[j,:])

    if A.size == 36:
        A4 = IsoMorphism66_3333(A)

    elif A.size == 81 and A.shape == (3,3,3,3):
        A4 = A

    TransformedA = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    TransformedA[i, j, k, l] += Q[m,i]*Q[n,j]*Q[o,k]*Q[p,l] * A4[m, n, o, p]
    if A.size == 36:
        TransformedA = IsoMorphism3333_66(TransformedA)

    return TransformedA

def Mandel2EngineeringNotation(A):

    B = np.zeros((6,6))

    for i in range(6):
        for j in range(6):

            if i < 3 and j >= 3:
                B[i,j] = A[i,j] / np.sqrt(2)

            elif i >= 3 and j < 3:
                B[i,j] = A[i,j] / np.sqrt(2)

            elif i >= 3 and j >= 3:
                B[i,j] = A[i,j] / 2

            else:
                B[i, j] = A[i, j]

    return B

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




def GetFabric(FileName):

    Text = open(FileName,'r').readlines()
    BVTV = float(Text[12].split('=')[1])
    eValues = np.array(Text[18].split(':')[1].split(),float)
    eVectors = np.zeros((3,3))
    for i in range(3):
        eVectors[i] = Text[19+i].split(':')[1].split()

    return eValues, eVectors, BVTV

def PlotFabricROI(ROI:np.array, eValues:np.array, eVectors:np.array, FileName:Path) -> None:

    """
    Plots a 3D ellipsoid representing a region of interest (ROI) with scaling based on the
    eigenvalues and eigenvectors provided. The ellipsoid is overlaid on a binary structure mesh,
    and the plot is generated with the ability to visualize the MIL (Mean Intercept Length) values.

    Parameters:
    -----------
    ROI (3D array): A 3D binary array representing the region of interest (ROI).
        
    eValues (1D array): A 1D array containing the eigenvalues of the fabric.
        
    eVectors (3D array) : A 2D array (shape: 3x3) containing the eigenvectors of the fabric.
        
    Returns:
    --------
    None
    """

    # Create a unit sphere and transform it to an ellipsoid
    Sphere = pv.Sphere(radius=ROI.shape[0]/2, theta_resolution=50, phi_resolution=50)

    # Scale the sphere by the square roots of the eigenvalues
    ScaleMatrix = np.diag(np.sqrt(eValues))
    TransformMatrix = np.matmul(eVectors, ScaleMatrix)

    # Transform the sphere points to ellipsoid points
    Points = np.matmul(Sphere.points, TransformMatrix.T)

    # Center the ellipsoid at the structure's midpoint
    Offset = np.array(ROI.shape) / 2
    EllispoidPoints = Points + Offset
    Ellispoid = pv.PolyData(EllispoidPoints, Sphere.faces)

    # Calculate the radius for each ellipsoid point to color by radius
    Radii = np.linalg.norm(Ellispoid.points - Offset, axis=1)
    Radii = (Radii - min(Radii)) / (max(Radii) - min(Radii))
    Radii = Radii * (max(eValues) - min(eValues)) + min(eValues)
    Ellispoid['MIL'] = Radii

    # Plotting
    sargs = dict(font_family='times', 
                    width=0.05,
                    height=0.75,
                    vertical=True,
                    position_x=0.9,
                    position_y=0.125,
                    title_font_size=30,
                    label_font_size=20
                    )
    
    pl = pv.Plotter(off_screen=True)
    actors = pl.add_volume(ROI,
                           cmap='bone',
                           show_scalar_bar=False, opacity=[0.005,0])
    actors.prop.interpolation_type = 'linear'
    pl.add_mesh(Ellispoid, scalars='MIL', cmap='jet', scalar_bar_args=sargs)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = 30
    pl.camera.zoom(1.0)
    pl.add_bounding_box(color=(0,0,0), line_width=1)
    # pl.add_axes(viewport=(0,0,0.25,0.25))
    pl.screenshot(FileName)
    # pl.show()

    return

def PlotStiffnessROI(ROI:np.array, StiffnessTensor:np.array, FileName:Path) -> None:

    """
    Plots a 3D ellipsoid representing a region of interest (ROI) with scaling based on the
    eigenvalues and eigenvectors provided. The ellipsoid is overlaid on a binary structure mesh,
    and the plot is generated with the ability to visualize the MIL (Mean Intercept Length) values.

    Parameters:
    -----------
    ROI (3D array): A 3D binary array representing the region of interest (ROI).
        
    eValues (1D array): A 1D array containing the eigenvalues of the fabric.
        
    eVectors (3D array) : A 2D array (shape: 3x3) containing the eigenvectors of the fabric.
        
    Returns:
    --------
    None
    """

    # Create a unit sphere and transform it to an ellipsoid
    Sphere = pv.Sphere(radius=1, theta_resolution=50, phi_resolution=50)

    # Compute elongation and bulk modulus
    I = np.eye(3)
    ElongationModulus = np.zeros(Sphere.points.shape)
    BulkModulus = np.zeros(len(Sphere.points))
    for p, Point in enumerate(Sphere.points):
        N = DyadicProduct(Point, Point)
        SN = Transform(StiffnessTensor, N)
        ElongationModulus[p] = FrobeniusProduct(N, SN) * Point
        BulkModulus[p] = FrobeniusProduct(I, SN)

    # # Scale the sphere by the square roots of the eigenvalues
    # Scale = ROI.shape[0]/2 / max(np.linalg.norm(ElongationModulus, axis=1))
    # Points = ElongationModulus * Scale

    # # Center the ellipsoid at the structure's midpoint
    # Offset = np.array(ROI.shape) / 2
    # EllispoidPoints = Points + Offset
    Ellispoid = pv.PolyData(ElongationModulus, Sphere.faces)
    Ellispoid['Bulk Modulus'] = BulkModulus

    # Plotting
    SArgs = dict(font_family='times', 
                 width=0.05,
                 height=0.75,
                 vertical=True,
                 position_x=0.85,
                 position_y=0.125,
                 title_font_size=30,
                 label_font_size=20
                 )
    
    BArgs = dict(font_family='times', 
                 font_size=30,
                 location='back',
                 n_xlabels=3,
                 n_ylabels=3,
                 n_zlabels=3,
                 all_edges=True,
                 fmt='%i',
                 xtitle='',
                 ytitle='',
                 ztitle='',
                 use_3d_text=False
                 )
    
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(Ellispoid, scalars='Bulk Modulus', cmap='jet', scalar_bar_args=SArgs)
    pl.camera_position = 'xz'
    pl.camera.roll = 0
    pl.camera.elevation = 30
    pl.camera.azimuth = 30
    pl.camera.zoom(0.8)
    pl.add_axes(viewport=(0,0,0.25,0.25),
                label_size=(0.065, 0.065),
                xlabel='e1',
                ylabel='e2',
                zlabel='e3')
    pl.show_bounds(**BArgs)
    pl.screenshot(FileName, return_img=False, scale=2)
    # pl.show()

    return


def BoxPlot(ArraysList:list, Labels=['', 'Y'], SetsLabels=None,
            Vertical=True, FigName=None, YLim=[], Ttest=None) -> None:
    
    """
    Save boxplot of a list of arrays
    Used for assessment on random effects and residuals
    """

    if Vertical == True:
        Width = 2.5 + len(ArraysList)
        Figure, Axis = plt.subplots(1,1, dpi=100, figsize=(Width,4.5))
    else:
        Height = len(ArraysList) - 0.5
        Figure, Axis = plt.subplots(1,1, dpi=100, figsize=(6.5,Height))

    for i, Array in enumerate(ArraysList):

        # Create random positions
        Array = np.sort(Array)
        Norm = norm.pdf(np.linspace(-3,3,len(Array)), scale=1.5)
        Norm = Norm / max(Norm)
        RandPos = np.random.normal(0,0.03,len(Array)) * Norm + i

        if Vertical == True:
            Axis.plot(RandPos - RandPos.mean() + i, Array, linestyle='none',
                        marker='o',fillstyle='none', color=(1,0,0), ms=5)
        else:
            Axis.plot(Array, RandPos - RandPos.mean() + i, linestyle='none',
                        marker='o',fillstyle='none', color=(1,0,0), ms=5)
            
        Axis.boxplot(Array, vert=Vertical, widths=0.35,
                    showmeans=True,meanline=True,
                    showfliers=False, positions=[i],
                    capprops=dict(color=(0,0,0)),
                    boxprops=dict(color=(0,0,0)),
                    whiskerprops=dict(color=(0,0,0),linestyle='--'),
                    medianprops=dict(color=(0,1,0)),
                    meanprops=dict(color=(0,0,1)))

    if Ttest:
        for i, A in enumerate(ArraysList[:-1]):

            # Perform t-test
            if Ttest == 'Rel':
                T_Tests = ttest_rel(np.array(ArraysList[i+1],float), np.array(A,float))
            else:
                T_Tests = ttest_ind(np.array(ArraysList[i+1],float), np.array(A,float))
            YLine = 1.05 * max(A.max(), ArraysList[i+1].max())
            Plot = Axis.plot([i+0.05, i+0.95], [YLine, YLine], color=(0,0,0), marker='|',linewidth=0.5)
            MarkerSize = Plot[0].get_markersize()
            
            # Mark significance level
            if T_Tests[1] < 0.001:
                Text = '***'
            elif T_Tests[1] < 0.01:
                Text = '**' 
            elif T_Tests[1] < 0.05:
                Text = '*'
            else:
                Text = 'n.s.'
            Axis.annotate(Text, xy=[i+0.5, YLine], ha='center',
                          xytext=(0, -1.5*MarkerSize), textcoords='offset points',)

            # Write confidence interveal
            CIl = round(T_Tests.confidence_interval()[0],1)
            CIu = round(T_Tests.confidence_interval()[1],1)
            Text = 'CI (' + str(CIl) + ',' + str(CIu) + ')'
            Axis.annotate(Text, xy=[i+0.5, YLine], ha='center',
                          xytext=(0, 1.2*MarkerSize), textcoords='offset points',)
            if i == 0:
                Max = YLine*1.05
            else:
                Max = max([Max, YLine*1.05])
            Axis.set_ylim([0.95*min([min(A)for A in ArraysList]), Max])
    
    Axis.plot([],linestyle='none',marker='o',fillstyle='none', color=(1,0,0), label='Data')
    Axis.plot([],color=(0,0,1), label='Mean', linestyle='--')
    Axis.plot([],color=(0,1,0), label='Median')
    Axis.set_xlabel(Labels[0])
    Axis.set_ylabel(Labels[1])

    if SetsLabels and Vertical==True:
        Axis.set_xticks(np.arange(len(SetsLabels)))
        Axis.set_xticklabels(SetsLabels, rotation=0)
    elif SetsLabels and Vertical==False:
        Axis.set_yticks(np.arange(len(SetsLabels)))
        Axis.set_yticklabels(SetsLabels, rotation=0)
    else:
        Axis.set_xticks([])

    if len(YLim) == 2:
        Axis.set_ylim(YLim)
    
    if Vertical == True:
        plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.125))
        plt.subplots_adjust(left=0.25, right=0.75)
    else:
        plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.25))
        plt.subplots_adjust(left=0.25, right=0.75,top=0.8)
    
    if FigName:
        plt.savefig(FigName, bbox_inches='tight', pad_inches=0.02, dpi=196)
    plt.show(Figure)

    return


#%% Main

def Main():

    DataPath = Path(__file__).parents[1] / 'Tests'
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])

    # Get fabric info
    eValues, eVectors, BVTV = GetFabric(DataPath / 'Fabric.fab')

    # Sort fabric
    Arg = np.argsort(eValues)
    eValues = eValues[Arg]
    eVectors = eVectors[Arg]
    m1, m2, m3 = eValues

    # Get homogenization stress results
    Isotropic = open(DataPath / 'Cube_Isotropic.out', 'r').readlines()
    Transverse = open(DataPath / 'Cube_Transverse.out', 'r').readlines()
    Fabric = open(DataPath / 'Cube_Fabric.out', 'r').readlines()

    IsoStress = np.zeros((6,6))
    TransStress = np.zeros((6,6))
    FabStress = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            IsoStress[i,j] = float(Isotropic[i+4].split()[j+1])
            TransStress[i,j] = float(Transverse[i+4].split()[j+1])
            FabStress[i,j] = float(Fabric[i+4].split()[j+1])

    IsoStiffness = np.zeros((6,6))
    TransStiffness = np.zeros((6,6))
    FabStiffness = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            IsoStiffness[i,j] = IsoStress[i,j] / Strain[i]
            TransStiffness[i,j] = TransStress[i,j] / Strain[i]
            FabStiffness[i,j] = FabStiffness[i,j] / Strain[i]

    


    

    # Symetrize matrix
    Stiffness = 1/2 * (Stiffness + Stiffness.T)

    # Write tensor into mandel notation
    Mandel = Engineering2MandelNotation(Stiffness)

    # Step 3: Transform tensor into fabric coordinate system
    I = np.eye(3)
    Q = np.array(eVectors)
    Transformed = TransformTensor(Mandel, I, Q)

    # Project onto orthotropy
    Orthotropic = np.zeros(Transformed.shape)
    for i in range(Orthotropic.shape[0]):
        for j in range(Orthotropic.shape[1]):
            if i < 3 and j < 3:
                Orthotropic[i, j] = Transformed[i, j]
            elif i == j:
                Orthotropic[i, j] = Transformed[i, j]

    # Get tensor back to engineering notation
    Stiffness = Mandel2EngineeringNotation(Orthotropic)

    # Build linear system
    Start, Stop = 12*s, 12*(s+1)
    Y[Start:Stop] = np.log([[Stiffness[0,0]],
                            [Stiffness[0,1]],
                            [Stiffness[0,2]],
                            [Stiffness[1,0]],
                            [Stiffness[1,1]],
                            [Stiffness[1,2]],
                            [Stiffness[2,0]],
                            [Stiffness[2,1]],
                            [Stiffness[2,2]],
                            [Stiffness[1,2]],
                            [Stiffness[2,0]],
                            [Stiffness[0,1]]])
    
    X[Start:Stop] = np.array([[1, 0, 0, np.log(BVTV), np.log(m1 ** 2)],
                            [0, 1, 0, np.log(BVTV), np.log(m1 * m2)],
                            [0, 1, 0, np.log(BVTV), np.log(m1 * m3)],
                            [0, 1, 0, np.log(BVTV), np.log(m2 * m1)],
                            [1, 0, 0, np.log(BVTV), np.log(m2 ** 2)],
                            [0, 1, 0, np.log(BVTV), np.log(m2 * m3)],
                            [0, 1, 0, np.log(BVTV), np.log(m3 * m1)],
                            [0, 1, 0, np.log(BVTV), np.log(m3 * m2)],
                            [1, 0, 0, np.log(BVTV), np.log(m3 ** 2)],
                            [0, 0, 1, np.log(BVTV), np.log(m2 * m3)],
                            [0, 0, 1, np.log(BVTV), np.log(m3 * m1)],
                            [0, 0, 1, np.log(BVTV), np.log(m1 * m2)]])
    

    # Solve linear system
    Parameters, R2adj, NE = OLS(X, Y)


if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--Sample', help='Sample main file name', type=str)
    Parser.add_argument('--OutputPath', help='Output path for the ROI and png image of the plot', default=Path(__file__).parents[1] / '02_Results/Scans')
    Parser.add_argument('--NROIs', help='Number of region of interests to extract', type=int, default=3)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main(Arguments)
