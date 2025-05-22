#%% !/usr/bin/env python3

Description = """
Compare numerical simularions to experimental results
obtained by resonant ultrasound spectroscopy (RUS)
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
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats.distributions import norm, t

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Read, Tensor

#%% Functions

def OLS(X,Y):
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y
    return B

def ExpVsSim_OLS(X, Y, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute residuals, variance, and covariance matrix
    Y_Obs = Y
    Y_Fit = np.array(X * B).ravel()
    ArgSort = np.argsort(Y_Fit)
    Residuals = Y - X*B
    DOFs = len(Y) - X.shape[1]
    Sigma = Residuals.T * Residuals / DOFs
    Cov = Sigma[0,0] * XTXi

    # Compute B confidence interval
    t_Alpha = t.interval(Alpha, DOFs)
    B_CI_Low = B.T + t_Alpha[0] * np.sqrt(np.diag(Cov))
    B_CI_Top = B.T + t_Alpha[1] * np.sqrt(np.diag(Cov))

    # Store parameters in data frame
    Parameters = pd.DataFrame(columns=[f'Parameter {i+1}' for i in range(len(B))])
    Parameters.loc['Value'] = [P[0] for P in np.array(B)]
    Parameters.loc['95% CI Low'] = [P for P in np.array(B_CI_Low)[0]]
    Parameters.loc['95% CI Top'] = [P for P in np.array(B_CI_Top)[0]]

    # Compute R2 and standard error of the estimate
    RSS = np.sum([R**2 for R in Residuals])
    SE = np.sqrt(RSS / DOFs)
    TSS = np.sum([R**2 for R in (Y - Y.mean())])
    RegSS = TSS - RSS
    R2 = RegSS / TSS

    # Compute R2adj and NE
    R2adj = 1 - RSS/TSS * (len(Y)-1)/(len(Y)-X.shape[1]-1)

    NE = []
    Step = 12
    for i in range(0,len(Y),Step):
        T_Obs = X[i:i+Step,1]
        T_Fit = Y[i:i+Step]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(np.min(np.concatenate([X[:,1],Y])),
                       np.max(np.concatenate([X[:,1],Y])), len(Y))
    # B_0 = np.sort(np.sqrt(np.diag(X * Cov * X.T)))
    # CI_Line_u = (Y_Fit + t_Alpha[0]) * B_0
    # CI_Line_o = (Y_Fit + t_Alpha[1]) * B_0

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    # Elements
    ii = np.tile([1,0,0,0,1,0,0,0,1,0,0,0],len(X)//Step).astype(bool)
    ij = np.tile([0,1,1,1,0,1,1,1,0,0,0,0],len(X)//Step).astype(bool)
    jj = np.tile([0,0,0,0,0,0,0,0,0,1,1,1],len(X)//Step).astype(bool)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.array(X[:,1]).ravel(), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(X[ii, 1], Y_Obs[ii], color=Colors[0], linestyle='none', marker='s')
    Axes.plot(X[ij, 1], Y_Obs[ij], color=Colors[1], linestyle='none', marker='o')
    Axes.plot(X[jj, 1], Y_Obs[jj], color=Colors[2], linestyle='none', marker='^')
    Axes.plot([], color=Colors[0], linestyle='none', marker='s', label=r'$\lambda_{ii}$')
    Axes.plot([], color=Colors[1], linestyle='none', marker='o', label=r'$\lambda_{ij}$')
    Axes.plot([], color=Colors[2], linestyle='none', marker='^', label=r'$\mu_{ij}$')
    Axes.plot(Line, Line, color=(0.8, 0.8, 0.8), linestyle='--', linewidth=1)
    Axes.plot(np.array(X[:,1]).ravel()[ArgSort], Y_Fit[ArgSort], color=(0, 0, 0), linewidth=1)
    Axes.annotate(r'N ROIs   : ' + str(len(Y)//Step), xy=(0.3, 0.1), xycoords='axes fraction')
    Axes.annotate(r'N Points : ' + str(len(Y)), xy=(0.3, 0.025), xycoords='axes fraction')
    Axes.annotate(r'$R^2_{ajd}$: ' + format(round(R2adj, 3),'.3f'), xy=(0.65, 0.1), xycoords='axes fraction')
    Axes.annotate(r'NE : ' + format(round(NE.mean(), 2), '.2f') + '$\pm$' + format(round(NE.std(), 2), '.2f'), xy=(0.65, 0.025), xycoords='axes fraction')
    Axes.annotate(r'y = ' + format(round(B[0,0],2), '.2f') + ' + ' + format(round(B[1,0], 2), '.2f') + 'x', xy=(0.3, 0.85), xycoords='axes fraction')
    Axes.set_ylabel('Simulation $\mathrm{\mathbb{S}}$ (MPa)')
    Axes.set_xlabel('Experiment $\mathrm{\mathbb{S}}$ (MPa)')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show(Figure)

    return Parameters, R2adj, NE

#%% Main

def Main():

    # Read CV data
    CVData = pd.read_csv(Path(__file__).parents[1] / 'Results' / 'CV.csv', index_col=[0,1])

    # Read RUS data
    DataPath = Path(__file__).parents[1] / 'Data'
    RUS = pd.read_csv(DataPath / '00_RUSData.csv')
    RUS['S12'] = RUS['S11'] - 2*RUS['S66']

    # Build RUS stiffness matrices
    mRUS = np.zeros((len(RUS),6,6))
    for i, Row in RUS.iterrows():
        mRUS[i,0,0] = Row['S11']
        mRUS[i,0,1] = Row['S12']
        mRUS[i,0,2] = Row['S13']
        mRUS[i,1,0] = Row['S12']
        mRUS[i,1,1] = Row['S11']
        mRUS[i,1,2] = Row['S13']
        mRUS[i,2,0] = Row['S13']
        mRUS[i,2,1] = Row['S13']
        mRUS[i,2,2] = Row['S33']
        mRUS[i,3,3] = Row['S44']
        mRUS[i,4,4] = Row['S44']
        mRUS[i,5,5] = Row['S66']

    # Define paths
    ResultsPath = Path(__file__).parents[1] / 'Results' / 'Homogenisation'
    FabricPath = Path(__file__).parents[1] / 'Results' /'Fabric'
    Folders = sorted([Folder.name for Folder in FabricPath.iterdir() if Folder.is_dir()])

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

                # Project onto orthotropy
                for i in range(6):
                    for j in range(6):
                        if i>2 or j>2:
                            if i != j:
                                S[f,r,i,j] = 0

    # Get average stiffness tensor per sample
    mRho = np.mean(Rho, axis=1)
    mFabric = np.mean(Fabric, axis=1)
    mIsotropic = np.mean(Isotropic, axis=1)
    mTransverse = np.mean(Transverse, axis=1)

    # Plot average ROI BV/TV vs RUS BV/TV
    X = np.ones((len(mRho),2))
    X[:,1] = 1-RUS['Rho']
    B = OLS(np.matrix(X), np.matrix(1-mRho).T)
    Min = 1-min([RUS['Rho'].min(), mRho.min()])
    Max = 1-max([RUS['Rho'].max(), mRho.max()])
    X = np.linspace(Min, Max)
    Y = B[0,0] + X * B[1,0]


    Figure, Axis = plt.subplots(1,1,dpi=200)
    Axis.plot(1-RUS['Rho'], 1-mRho, linestyle='none', color=(1,0,0),
              marker='o', label='Data')
    Axis.plot(X, Y, color=(0,0,0), linestyle='--', label='Fit')
    Axis.plot([Min,Max], [Min,Max], color=(0.4,0.4,0.4), linestyle='--', label='1:1 line')
    Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)} x', (0.25,0.05), color=(0,0,0))
    Axis.set_xlabel(r'RUS 1-$\rho$ (-)')
    Axis.set_ylabel(r'Average ROI 1-$\rho$ (-)')
    plt.legend()
    plt.savefig(Path(__file__).parents[1] / 'Results/ExpSim_Rho.png')
    plt.show(Figure)

    # Compare with RUS
    X = np.matrix(np.ones((len(mRho)*12, 2)))
    Y = np.matrix(np.zeros((len(mRho)*12, 1)))
    for f in range(len(mRho)):
        
        Start, Stop = 12*f, 12*(f+1)
        
        X[Start:Stop,1] = [[mRUS[f][0,0]/1E3],
                           [mRUS[f][0,1]/1E3],
                           [mRUS[f][0,2]/1E3],
                           [mRUS[f][1,0]/1E3],
                           [mRUS[f][1,1]/1E3],
                           [mRUS[f][1,2]/1E3],
                           [mRUS[f][2,0]/1E3],
                           [mRUS[f][2,1]/1E3],
                           [mRUS[f][2,2]/1E3],
                           [mRUS[f][3,3]/1E3],
                           [mRUS[f][4,4]/1E3],
                           [mRUS[f][5,5]/1E3]]
        
        Y[Start:Stop] = [[mTransverse[f][0,0]],
                         [mTransverse[f][0,1]],
                         [mTransverse[f][0,2]],
                         [mTransverse[f][1,0]],
                         [mTransverse[f][1,1]],
                         [mTransverse[f][1,2]],
                         [mTransverse[f][2,0]],
                         [mTransverse[f][2,1]],
                         [mTransverse[f][2,2]],
                         [mTransverse[f][3,3]],
                         [mTransverse[f][4,4]],
                         [mTransverse[f][5,5]/2]]
       
    FName = Path(__file__).parents[1] / 'Results/ExpSim_S_Raw.png'
    Parameters, R2adj, NE = ExpVsSim_OLS(X, Y/1E3, FName=str(FName))

    # Project tensors onto transverse isotropy
    mProjected, NE = Tensor.TransProj(mTransverse)

    # Compare with RUS
    X = np.matrix(np.ones((len(mRho)*12, 2)))
    Y = np.matrix(np.zeros((len(mRho)*12, 1)))
    for f in range(len(mRho)):
        
        Start, Stop = 12*f, 12*(f+1)
        
        X[Start:Stop,1] = [[mRUS[f][0,0]/1E3],
                           [mRUS[f][0,1]/1E3],
                           [mRUS[f][0,2]/1E3],
                           [mRUS[f][1,0]/1E3],
                           [mRUS[f][1,1]/1E3],
                           [mRUS[f][1,2]/1E3],
                           [mRUS[f][2,0]/1E3],
                           [mRUS[f][2,1]/1E3],
                           [mRUS[f][2,2]/1E3],
                           [mRUS[f][3,3]/1E3],
                           [mRUS[f][4,4]/1E3],
                           [mRUS[f][5,5]/1E3]]
        
        Y[Start:Stop] = [[mProjected[f][0,0]],
                         [mProjected[f][0,1]],
                         [mProjected[f][0,2]],
                         [mProjected[f][1,0]],
                         [mProjected[f][1,1]],
                         [mProjected[f][1,2]],
                         [mProjected[f][2,0]],
                         [mProjected[f][2,1]],
                         [mProjected[f][2,2]],
                         [mProjected[f][3,3]],
                         [mProjected[f][4,4]],
                         [mProjected[f][5,5]/2]]
       
    FName = Path(__file__).parents[1] / 'Results/ExpSim_S_Proj.png'
    Parameters, R2adj, NE = ExpVsSim_OLS(X, Y/1E3, FName=str(FName))

    # Filter out samples with high CV
    F = np.all(CV < 0.263, axis=1)
    F12 = np.repeat(F, 12)
    FName = Path(__file__).parents[1] / 'Results/ExpSim_S_Filt.png'
    Parameters, R2adj, NE = ExpVsSim_OLS(X[F12], Y[F12]/1E3, FName=str(FName))
    


    # Compare degree of anisotropy
    X = np.ones((len(mRho[F]),2))
    X[:,1] = 1-mRho[F]
    Args = np.argsort(1-mRho[F])
    Colors = [(0,0,1),(1,0,0)]
    Labels = ['Experiment','Simulation']

    Figure, Axis = plt.subplots(1,1,dpi=200)
    for i, (Y, C, L) in enumerate(zip([mRUS[F], mTransverse[F]], Colors, Labels)):
        Y = Y[:,2,2] / Y[:,0,0]
        B = OLS(np.matrix(X),np.matrix(Y).T)
        Axis.plot(1-mRho[F], Y, color=C, label=L,
                marker='o', fillstyle='none', linestyle='none')
        Axis.plot((1-mRho[F])[Args], B[0,0] + (1-mRho[F])[Args]*B[1,0], color=C, linestyle='--')
        Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x', color=C,
                      xy=(0.0475, 1.76-i*0.05))
    Axis.set_xlabel(r'1-$\rho$ (-)')
    Axis.set_ylabel(r'$S_{33}$ / $S_{11}$ (-)')
    plt.legend()
    plt.savefig(Path(__file__).parents[1] / 'Results/ExpSim_AniS.png')
    plt.show(Figure)

    # Same for compliance
    E_Sim = np.zeros((len(mTransverse),6,6))
    for i, S in enumerate(mTransverse):
        E_Sim[i] = np.linalg.inv(S)

    E_RUS = np.zeros((len(mRUS),6,6))
    for i, S in enumerate(mRUS):
        E_RUS[i] = np.linalg.inv(S)

    X = np.ones((len(mRho[F]),2))
    X[:,1] = 1-mRho[F]
    Args = np.argsort(1-mRho[F])
    Colors = [(0,0,1),(1,0,0)]
    Labels = ['Experiment','Simulation']

    Figure, Axis = plt.subplots(1,1,dpi=200)
    for i, (Y, C, L) in enumerate(zip([E_RUS[F], E_Sim[F]], Colors, Labels)):
        Y = Y[:,0,0] / Y[:,2,2]
        B = OLS(np.matrix(X),np.matrix(Y).T)
        Axis.plot(1-mRho[F], Y, color=C, label=L,
                marker='o', fillstyle='none', linestyle='none')
        Axis.plot((1-mRho[F])[Args], B[0,0] + (1-mRho[F])[Args]*B[1,0], color=C, linestyle='--')
        Axis.annotate(f'y = {round(B[0,0],2)} + {round(B[1,0],2)}x', color=C,
                      xy=(0.2, 1.45-i*0.05))
    Axis.set_xlabel(r'1-$\rho$ (-)')
    Axis.set_ylabel(r'$E_{33}$ / $E_{11}$ (-)')
    plt.legend()
    plt.savefig(Path(__file__).parents[1] / 'Results/ExpSim_AniE.png')
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
