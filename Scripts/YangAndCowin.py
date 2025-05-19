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
from math import sec
from pathlib import Path
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def TransverseTensor(C11, C12, C13, C33, C44):

    if type(C11) == float or type(C11) == np.float64:
        C = np.array([[C11, C12, C13,     0,      0,       0],
                      [C12, C11, C13,     0,      0,       0],
                      [C13, C13, C33,     0,      0,       0],
                      [  0,   0,   0, 2*C44,      0,       0],
                      [  0,   0,   0,     0,  2*C44,       0],
                      [  0,   0,   0,     0,      0, C11-C12]])
        

    else:
        C = np.zeros((len(C11), 6, 6))
        C[:,0,0] = C11
        C[:,1,1] = C11
        C[:,0,1] = C12
        C[:,1,0] = C12
        C[:,0,2] = C13
        C[:,1,2] = C13
        C[:,2,0] = C13
        C[:,2,1] = C13
        C[:,2,2] = C33
        C[:,3,3] = 2*C44
        C[:,4,4] = 2*C44
        C[:,5,5] = C11 - C12

    return C

def ProjectTransverse(C):
    
    """
    Vectorized projection of m 6x6 stiffness matrices to transverse isotropy.
    
    Parameters:
        C : ndarray, shape (m, 6, 6) — symmetric stiffness matrices
    
    Returns:
        ti_params : ndarray, shape (m, 5) — [C11, C12, C13, C33, C44]
        errors : ndarray, shape (m,) — Frobenius norm error ||C - C_TI||
    """
        
    # Least-squares design matrix
    X = np.array([[1.0, 0.0],
                  [1.0, 0.0],
                  [0.0, 1.0],
                  [0.0, 1.0],
                  [1.0, -1.0]])

    # Targets for least-squares
    Y = np.stack([C[:, 0, 0],
                  C[:, 1, 1],
                  C[:, 0, 1],
                  C[:, 1, 0],
                  C[:, 5, 5]], axis=-1)

    # Solve least-squares in batch: X = (A^T A)^{-1} A^T y
    XTXi = np.linalg.inv(X.T @ X)
    B = (XTXi @ X.T @ Y.transpose(1, 0)).T  # shape (m, 2)

    # Extract TI constants
    C11, C12 = B[:, 0], B[:, 1]
    C13 = (C[:, 0, 2] + C[:, 1, 2]) / 2
    C33 = C[:, 2, 2]
    C44 = (C[:, 3, 3] + C[:, 4, 4]) / 2
    Parameters = np.stack([C11, C12, C13, C33, C44/2], axis=1)

    # Build transverse isotropic tensors
    T = TransverseTensor(C11, C12, C13, C33, C44/2)

    # Compute Frobenius norm error
    NE = np.sum((C - T)**2, axis=(1,2)) / np.sum(C**2, axis=(1,2))
    NE = np.sqrt(NE)
    
    return Parameters, NE

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

def RotateTensors(C, Q):

    A = np.zeros(C.shape)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    A[:,i,j] += Q[:,i,k] * Q[:,j,l] * C[:,k,l]

    return A

def EigenValuesFunction(Rho, a, b, c):
    return a*Rho**b

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

                # Use Mandel notation
                S[f,r] = Tensor.Engineering2MandelNotation(S[f,r])

                # Compute corresponding compliance
                Compliance[f,r] = np.linalg.inv(S[f,r])

    # Reshape arrays
    Isotropic = np.reshape(Isotropic, (-1, 6, 6))
    Transverse = np.reshape(Transverse, (-1, 6, 6))
    Compliance = np.reshape(Compliance, (-1, 6, 6))
    Fabric = np.reshape(Fabric, (-1, 3))
    Rho = np.reshape(Rho, (-1))
    CV = np.reshape(CV, (-1))

    # Plot density distrbution
    FName = Path(__file__).parents[1] / '_Results/DensityDistribution.png'
    PlotHistogram(Rho, str(FName))

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

    # Project tensors onto transverse isotropy
    Projected, NE = Tensor.TransProj(Transverse)

    # Compute eigenvalues
    Alpha = np.arctan(2**0.5 / (4*Projected[:,0,2]) * (Projected[:,0,0] + Projected[:,0,1] - Projected[:,2,2]))
    L1 = Projected[:,2,2] + 2**0.5 * Projected[:,0,2] * (np.tan(Alpha) + 1/np.cos(Alpha))
    L2 = Projected[:,2,2] + 2**0.5 * Projected[:,0,2] * (np.tan(Alpha) - 1/np.cos(Alpha))
    L34 = Projected[:,0,0] - Projected[:,0,1]
    L56 = Projected[:,3,3]

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

    # Test building tensor back
    Test = L1[0] * Tensor.Dyadic(E1[0], E1[0])
    Test += L2[0] * Tensor.Dyadic(E2[0], E2[0])
    Test += L34[0] * Tensor.Dyadic(E3[0], E3[0])
    Test += L34[0] * Tensor.Dyadic(E4[0], E4[0])
    Test += L56[0] * Tensor.Dyadic(E5[0], E5[0])
    Test += L56[0] * Tensor.Dyadic(E6[0], E6[0])
    Test = Tensor.IsoMorphism3333_66(Test)


    # Spectral decomposition of all tensors
    Moduli, Bases = np.linalg.eig(Projected[0])
    Test = np.zeros((6,6))
    for i in range(6):
        Test += Moduli[i] * np.outer(Bases.T[i], Bases.T[i])

    # Sort eigenvectors and eigenvalues
    Sort = [0, 2, 1, 3, 4, 5]
    Bases = np.transpose(Bases, (0,2,1))
    Bases = Bases[:,Sort]
    Bases = np.transpose(Bases, (0,2,1))
    Moduli = Moduli[:, Sort]

    Ref = np.array([[1,  1,  1,  0,  0,  0],
                    [1,  1, -1,  0,  0,  0],
                    [1, -1,  0,  0,  0,  0],
                    [0,  0,  0,  1,  0,  0],
                    [0,  0,  0,  0,  1,  0],
                    [0,  0,  0,  0,  0,  1]])    

    N_NA = np.sum(Bases, axis=0) / len(Bases)
    
    # J
    J_NA = np.zeros((6,6))
    for i in range(6):
        J_NA += np.outer(N_NA[i], N_NA[i])
    L, D = np.linalg.eig(J_NA)
    J_NA = D @ np.diag(L**(-1/2)) @ D.T

    # Average tensor
    N_AVG = np.zeros((6,6))
    for i in range(6):
        N_AVG[i] = J_NA @ N_NA[i]

    # Compute rotation matrices
    Q = np.einsum('mij,ik -> mjk', Bases, N_AVG)

    # Rotate tensors
    C = RotateTensors(Projected, Q)

    # Compute corresponding eigenvalues
    eValues, eVectors = np.linalg.eig(C)
    eVectors = np.transpose(eVectors, (0,2,1))
    eVectors = eVectors[:,Sort]
    eVectors = np.transpose(eVectors, (0,2,1))
    eValues = eValues[:, Sort]
    Lambdas = eValues[:,:-2]

    # Plot eigen values as function of BVTV
    X = np.linspace(Rho.min(), Rho.max(), 100)
    Figure, Axis = plt.subplots(2,2, sharex=True, figsize=(10,7))
    for i in range(2):
        for j in range(2):
            Y = Lambdas[:,2*i+j]/1E3
            P = curve_fit(EigenValuesFunction, Rho, Y)
            YFit = P[0][0] + P[0][1] * X ** P[0][2]
            Residuals = Y - (P[0][0] + P[0][1] * Rho ** P[0][2])
            RSS = np.sum([R**2 for R in Residuals])
            TSS = np.sum([R**2 for R in (Y - Y.mean())])
            RegSS = TSS - RSS
            R2 = RegSS / TSS
            Axis[i,j].plot(Rho, Y, marker='o', linestyle='none', color=(1,0,0))
            Axis[i,j].plot(X, YFit, color=(0,0,1))
            Axis[i,j].set_xlabel(r'$\rho$ (-)')
            Axis[i,j].set_ylabel(f'$\Lambda_{2*i+j+1}$ (GPa)')
            Axis[i,j].annotate(f'$\Lambda_{2*i+j+1}$ = {round(P[0][0])} + {round(P[0][1])}' + r'$\rho$' + '$^{%.2f}$'%(round(P[0][2],2)), xy=(0.1, 0.9), xycoords='axes fraction')
            Axis[i,j].annotate(r'$R^2$: ' + format(round(R2, 3),'.3f'), xy=(0.1, 0.8), xycoords='axes fraction')
    plt.tight_layout()
    plt.show(Figure)

    # Build again tensor



    # Put average basis into second rank tensors
    e1 = N_AVG[0] # Dilatational mode 1
    E1 = np.array([[[e1[0], e1[5], e1[4]],
                    [e1[5], e1[1], e1[3]],
                    [e1[4], e1[3], e1[2]]]])
    np.linalg.det(E1)
    
    e2 = N_AVG[1] # Dilatational mode 2
    E2 = np.array([[[e2[0], e2[5], e2[4]],
                    [e2[5], e2[1], e2[3]],
                    [e2[4], e2[3], e2[2]]]])
    
    e3 = N_AVG[2] # Isochoric mode 1
    E3 = np.array([[[e3[0], e3[5], e3[4]],
                    [e3[5], e3[1], e3[3]],
                    [e3[4], e3[3], e3[2]]]])
    
    e4 = N_AVG[3] # Isochoric mode 2
    E4 = np.array([[[e4[0], e4[5], e4[4]],
                    [e4[5], e4[1], e4[3]],
                    [e4[4], e4[3], e4[2]]]])
    
    e5 = N_AVG[4] # Isochoric mode 2 bis
    E5 = np.array([[[e5[0], e5[5], e5[4]],
                    [e5[5], e5[1], e5[3]],
                    [e5[4], e5[3], e5[2]]]])
    
    e6 = N_AVG[5] # Isochoric mode 1 bis
    E6 = np.array([[[e6[0], e6[5], e6[4]],
                    [e6[5], e6[1], e6[3]],
                    [e6[4], e6[3], e6[2]]]])

    

    # Test methods of Cowin
    Test = np.array([[9.86, 2.773, 3.147, -0.115, -0.581, -2.072],
                     [2.773, 5.952, 2.346, -0.926, 0.382, -1.267],
                     [3.147, 2.346, 9.35, -0.851, 0.315, -0.662],
                     [-0.115, -0.926, -0.851, 2.805, -1.08, 0.399],
                     [-0.581, 0.382, 0.315, -1.08, 3.712, -0.168],
                     [-2.072, -1.267, -0.662, 0.399, -0.168, 3.016]])
    
    Test = Tensor.Engineering2MandelNotation(Test)
    
    # Permute coordinates
    Q = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    Q = Q66(Q)
    Test = Q.T @ Test @ Q
    
    # Change of coordinates
    Test = np.array([[9.86, 3.147, 2.773, -0.162, -2.93, -0.821],
                     [3.147, 9.35, 2.346, -1.204, -0.936, 0.445],
                     [2.773, 2.346, 5.952, -1.309, -1.791, 0.54],
                     [-0.162, -1.204, -1.309, 5.61, 0.798, -2.159],
                     [-2.93, -0.936, -1.791, 0.798, 6.032, -0.335],
                     [-0.821, 0.445, 0.54, -2.159, -0.335, 7.425]])
    
    Moduli, Base = np.linalg.eig(Test)

    # Correct values to mimic Cowin and Yang calculations
    Base = np.array([[0.64, 0.542, 0.396, -0.162, -0.333, 0.047],
                     [-0.278, 0.773, -0.117, -0.11, 0.5, -0.223],
                     [-0.238, 0.0068, 0.602, 0.712, -0.039, 0.268],
                     [0.472, 0.036, -0.511, 0.398, 0.242, 0.546],
                     [0.313, -0.287, 0.441, -0.216, 0.76, -0.019],
                     [-0.368, 0.157, 0.106, -0.499, 0.026, 0.761]])

    # # comparer angle
    # # produit scalaire entre les deux vecteurs


    # Nominal averages
    N_NA = np.zeros((6,6))
    N_NA[0] = [0.791, 0.396, 0.237, 0.0017, 0.031, 0.0067]
    N_NA[1] = [-0.365, 0.703, 0.144, 0.026, -0.035, -0.032]
    N_NA[2] = [-0.137, -0.238, 0.759, 0.011, 0.033, -0.0006]
    N_NA[3] = [0.0018, -0.0064, -0.015, 0.797, 0.014, -0.0064]
    N_NA[4] = [-0.039, 0.028, -0.05, -0.029, 0.753, -0.003]
    N_NA[5] = [-0.032, 0.02, 0.0048, -0.026, 0.0008, 0.773]

    # J
    J_NA = np.zeros((6,6))
    for i in range(6):
        J_NA += np.outer(N_NA[i], N_NA[i])
    L, D = np.linalg.eig(J_NA)
    J_NA = D @ np.diag(L**(-1/2)) @ D.T

    # Average tensor
    N_AVG = np.zeros((6,6))
    for i in range(6):
        N_AVG[i] = J_NA @ N_NA[i]

    # Compute rotation matrix¨
    Q = np.zeros((6,6))
    for i in range(6):
        Q += np.outer(Base[i], N_AVG[i])

    Q = np.array([[0.725, 0.088, -0.131, 0.458, 0.355, -0.338],
                  [0.116, 0.884, 0.307, 0.065, -0.296, 0.136],
                  [0.272, -0.078, 0.643, -0.523, 0.471, 0.101],
                  [-0.275, -0.202, 0.648, 0.429, -0.176, -0.498],
                  [-0.552, 0.329, -0.079, 0.223, 0.73, 0.0006],
                  [0.083, -0.235, 0.219, 0.528, 0.031, 0.781]])

    # Rotate tensor
    C = Q.T @ Test @ Q
    print(C)

    # Other way, compute averaged moduli (not working?)
    M_AVG = np.zeros(Moduli.shape)
    for j in range(6):
        for i in range(6):
            M_AVG[j] += Moduli[i] * np.dot(Base[i], N_AVG[j])**2

    C = np.zeros((6,6))
    for i in range(6):
        C += M_AVG * np.outer(Base.T[i], Base.T[i])
    print(C)

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
