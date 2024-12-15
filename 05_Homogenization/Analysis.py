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

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats.distributions import norm, t

sys.path.append(str(Path(__file__).parents[1]))
from Utils import Time, Read, Image, Tensor

#%% Functions

def OLS(X, Y, Alpha=0.95, FName=''):

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
    if len(FName) > 0:
        plt.savefig(FName)
    plt.show()

    return Parameters, R2adj, NE

def OLS2(X, Y, Alpha=0.95, FName=''):

    # Solve linear system
    XTXi = np.linalg.inv(X.T * X)
    B = XTXi * X.T * Y

    # Compute residuals, variance, and covariance matrix
    Y_Obs = Y
    Y_Fit = X * B
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
    for i in range(0,len(Y),12):
        T_Obs = Y_Obs[i:i+12]
        T_Fit = Y_Fit[i:i+12]
        Numerator = np.sum([T**2 for T in (T_Obs-T_Fit)])
        Denominator = np.sum([T**2 for T in T_Obs])
        NE.append(np.sqrt(Numerator/Denominator))
    NE = np.array(NE)


    # Prepare data for plot
    Line = np.linspace(np.min(np.concatenate([X,Y])),
                       np.max(np.concatenate([X,Y])), len(Y))
    # B_0 = np.sort(np.sqrt(np.diag(X * Cov * X.T)))
    # CI_Line_u = np.exp(Line + t_Alpha[0] * B_0)
    # CI_Line_o = np.exp(Line + t_Alpha[1] * B_0)

    # Plots
    DPI = 500
    Colors=[(0,0,1),(0,1,0),(1,0,0)]

    # Elements
    ii = np.tile([1,0,0,1,0,1,0,0,0],len(X)//9).astype(bool)
    ij = np.tile([0,1,1,0,1,0,0,0,0],len(X)//9).astype(bool)
    jj = np.tile([0,0,0,0,0,0,1,1,1],len(X)//9).astype(bool)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI)
    # Axes.fill_between(np.exp(Line), CI_Line_u, CI_Line_o, color=(0.8,0.8,0.8))
    Axes.plot(X[ii, 0], Y_Obs[ii],
              color=Colors[0], linestyle='none', marker='s')
    Axes.plot(X[ij, 0], Y_Obs[ij],
              color=Colors[1], linestyle='none', marker='o')
    Axes.plot(X[jj, 0], Y_Obs[jj],
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
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    if len(FName) > 0:
        plt.savefig(FName)
    plt.close(Figure)

    return Parameters, R2adj, NE

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

    # Read RUS data
    DataPath = Path(__file__).parents[1] / '00_Data'
    Data = pd.read_excel(DataPath / 'Cij.xls')

    # List fabric files
    FabricPath = Path(__file__).parents[1] / '01_Fabric/Results'
    FabricFiles = [F.name for F in FabricPath.iterdir() if F.name.endswith('.fab')]

    # List simulations output files
    ResultsPath =Path(__file__).parent / 'Elasticity'
    Folders = [F for F in ResultsPath.iterdir() if F.is_dir()]
    ROIs = []
    for x in range(2):
        for y in range(2):
            for z in range(4):
                    ROIs.append(f'ROI_{x+1}{y+1}{z+1}')

    # List complete folders
    cFolders = []
    for Folder in Folders:
        Outs = [F for F in Folder.iterdir() if F.name.endswith('.out')]
        if len(Outs) == 32:
            cFolders.append(Folder)

    # Collect data
    Strain = np.array([0.001, 0.001, 0.001, 0.002, 0.002, 0.002])
    Isotropic = np.zeros((len(cFolders), len(ROIs), 6, 6))
    Transverse = np.zeros((len(cFolders), len(ROIs), 6, 6))
    RUS = np.zeros((len(cFolders), 6, 6))
    Fabric = np.zeros((len(cFolders), 3))
    for f, Folder in enumerate(cFolders):

        # Iterate over each ROI
        for r, ROI in enumerate(ROIs):

            # Get simulation results
            for S, sType in zip([Isotropic, Transverse],['Isotropic','Transverse']):

                # Get homogenization stress results
                File = open(Folder / (ROI + f'_{sType}.out'), 'r').readlines()
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

        # Get sample fabric in original resolution
        F1 = np.array([Folder.name[:4] in File for File in FabricFiles])
        F2 = np.array([Folder.name[5:8] in File for File in FabricFiles])
        F3 = np.array([Folder.name[-1] in File for File in FabricFiles])
        Idx = np.argwhere(F1 & F2 & F3)[0][0]
        FabricFile = FabricFiles[Idx]
        FabricData = Read.Fabric(FabricPath / FabricFile)
        Fabric[f] = FabricData[0]

        # Get RUS data
        F1 = Data['sample'] == Folder.name[:-2]
        F2 = Data['Octant'] == Folder.name[-1]

        # Build stiffness tensor
        C11 = Data[F1&F2]['RUS:C11'].values[0]
        C33 = Data[F1&F2]['RUS:C33'].values[0]
        C44 = Data[F1&F2]['RUS:C44'].values[0]
        C66 = Data[F1&F2]['RUS:C66'].values[0]
        C13 = Data[F1&F2]['RUS:C13'].values[0]
        C12 = C11 - 2*C66

        C = np.array([[C11, C12, C13,   0,   0,   0],
                        [C12, C11, C13,   0,   0,   0],
                        [C13, C13, C33,   0,   0,   0],
                        [  0,   0,   0, C44,   0,   0],
                        [  0,   0,   0,   0, C44,   0],
                        [  0,   0,   0,   0,   0, C66]])

        RUS[f] = C

    # Get average stiffness tensor per sample
    mIsotropic = np.mean(Isotropic, axis=1)
    mTransverse = np.mean(Transverse, axis=1)

    # Compute anisotropy
    aRUS = np.zeros((len(cFolders),3))
    aFab = np.zeros((len(cFolders),3))
    aIso = np.zeros((len(cFolders),3))
    aTra = np.zeros((len(cFolders),3))
    for A, C in zip([aRUS, aIso, aTra],[RUS, mIsotropic, mTransverse]):
        A[:,0] = C[:,1,1] / C[:,0,0]
        A[:,1] = C[:,2,2] / C[:,1,1]
        A[:,2] = C[:,2,2] / C[:,0,0]
    aFab[:,0] = Fabric[:,1] / Fabric[:,0]
    aFab[:,1] = Fabric[:,2] / Fabric[:,1]
    aFab[:,2] = Fabric[:,2] / Fabric[:,0]

    # Plot anisotropies
    Colors = [(0,0,0),(1,0,0),(0,0,1),(1,0,1)]
    Labels = ['RUS','Fabric','Isotropic','Transverse']
    FName = Path(__file__).parent / 'Plots/Anisotropy.png'
    Figure, Axis = plt.subplots(1,1)
    for c, C in enumerate([aRUS, aFab, aIso, aTra]):
        Axis.plot(np.arange(3)+1, C[0], marker='o', color=Colors[c], label=Labels[c])
        for A in C[1:]:
            Axis.plot(np.arange(3)+1, A, marker='o', color=Colors[c])
    Axis.set_xticks(np.arange(3)+1)
    Axis.set_xticklabels(['$e_2$/$e_1$', '$e_3$/$e_2$', '$e_3$/$e_1$'])
    Axis.set_xlabel('Directions')
    Axis.set_ylabel('Anisotropy')
    plt.legend()
    plt.savefig(FName)
    plt.close(Figure)


    # Build linear system
    X = np.matrix(np.ones((len(cFolders)*9, 1)))
    Y = np.matrix(np.zeros((len(cFolders)*9, 1)))
    for f in range(len(cFolders)):
        
        Start, Stop = 9*f, 9*(f+1)
        X[Start:Stop] = [[RUS[f][0,0]],
                         [RUS[f][0,1]],
                         [RUS[f][0,2]],
                         [RUS[f][1,1]],
                         [RUS[f][1,2]],
                         [RUS[f][2,2]],
                         [RUS[f][3,3]],
                         [RUS[f][4,4]],
                         [RUS[f][5,5]]]
        
        Y[Start:Stop] = [[mIsotropic[f][0,0]],
                         [mIsotropic[f][0,1]],
                         [mIsotropic[f][0,2]],
                         [mIsotropic[f][1,1]],
                         [mIsotropic[f][1,2]],
                         [mIsotropic[f][2,2]],
                         [mIsotropic[f][3,3]],
                         [mIsotropic[f][4,4]],
                         [mIsotropic[f][5,5]]]
       
    FName = Path(__file__).parent / 'Plots/Elasticity_IsoRUS.png'
    Parameters, R2adj, NE = OLS2(np.log(X), np.log(Y/1E3), FName=str(FName))

    for f in range(len(cFolders)):
        Start, Stop = 9*f, 9*(f+1)
        Y[Start:Stop] = [[mTransverse[f][0,0]],
                         [mTransverse[f][0,1]],
                         [mTransverse[f][0,2]],
                         [mTransverse[f][1,1]],
                         [mTransverse[f][1,2]],
                         [mTransverse[f][2,2]],
                         [mTransverse[f][3,3]],
                         [mTransverse[f][4,4]],
                         [mTransverse[f][5,5]]]
       
    FName = Path(__file__).parent / 'Plots/Elasticity_TraRUS.png'
    Parameters, R2adj, NE = OLS2(np.log(X), np.log(Y/1E3), FName=str(FName))


    #     # Get fabric
    #     F1 = np.array([Folder.name[:4] in File for File in FabricFiles])
    #     F2 = np.array([Folder.name[5:8] in File for File in FabricFiles])
    #     F3 = np.array([Folder.name[-1] in File for File in FabricFiles])
    #     Idx = np.argwhere(F1 & F2 & F3)[0][0]
    #     FabricFile = FabricFiles[Idx]
    #     Fabric = Read.Fabric(FabricPath / FabricFile)
    #     m1, m2, m3 = Fabric[0]
    #     BVTV = Fabric[2]

    #     # Get corresponding stiffness
    #     Stiffness = mIsotropic[f]
        
    #     Start, Stop = 12*f, 12*(f+1)
    #     Y[Start:Stop] = np.log([[Stiffness[0,0]],
    #                             [Stiffness[0,1]],
    #                             [Stiffness[0,2]],
    #                             [Stiffness[1,0]],
    #                             [Stiffness[1,1]],
    #                             [Stiffness[1,2]],
    #                             [Stiffness[2,0]],
    #                             [Stiffness[2,1]],
    #                             [Stiffness[2,2]],
    #                             [Stiffness[1,2]],
    #                             [Stiffness[2,0]],
    #                             [Stiffness[0,1]]])
        
    #     X[Start:Stop] = np.array([[1, 0, 0, np.log(BVTV), np.log(m1 ** 2)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m1 * m2)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m1 * m3)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m2 * m1)],
    #                               [1, 0, 0, np.log(BVTV), np.log(m2 ** 2)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m2 * m3)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m3 * m1)],
    #                               [0, 1, 0, np.log(BVTV), np.log(m3 * m2)],
    #                               [1, 0, 0, np.log(BVTV), np.log(m3 ** 2)],
    #                               [0, 0, 1, np.log(BVTV), np.log(m2 * m3)],
    #                               [0, 0, 1, np.log(BVTV), np.log(m3 * m1)],
    #                               [0, 0, 1, np.log(BVTV), np.log(m1 * m2)]])

    # # Build linear system
    # X = np.matrix(np.ones((len(cFolders)*9, 1)))
    # Y = np.matrix(np.zeros((len(cFolders)*9, 1)))
    # for f, Folder in enumerate(cFolders):
        
    #     Start, Stop = 9*f, 9*(f+1)
    #     X[Start:Stop] = [[RUS[f][0,0]],
    #                      [RUS[f][0,1]],
    #                      [RUS[f][0,2]],
    #                      [RUS[f][1,1]],
    #                      [RUS[f][1,2]],
    #                      [RUS[f][2,2]],
    #                      [RUS[f][3,3]],
    #                      [RUS[f][4,4]],
    #                      [RUS[f][5,5]]]
        
    #     Y[Start:Stop] = [[mIsotropic[f][0,0]],
    #                      [mIsotropic[f][0,1]],
    #                      [mIsotropic[f][0,2]],
    #                      [mIsotropic[f][1,1]],
    #                      [mIsotropic[f][1,2]],
    #                      [mIsotropic[f][2,2]],
    #                      [mIsotropic[f][3,3]],
    #                      [mIsotropic[f][4,4]],
    #                      [mIsotropic[f][5,5]]]
       
    # FName = Path(__file__).parent / 'Plots/Elasticity_IsoRUS.png'
    # Parameters, R2adj, NE = OLS2(np.log(X), np.log(Y/1E3), FName=str(FName))

    # TransverseNorm = Tensor.Norm(mTransverse)
    # IsotropicNorm = Tensor.Norm(mIsotropic)
    # mIsotropic / IsotropicNorm * TransverseNorm

    # # Get Experiment RUS results
    # Octant = 'L' if '2L' in Sample else 'M'
    # Year = Sample[16:20]
    # Number = Sample[12:15]
    # F1 = RUSData['Octant'] == Octant
    # F2 = RUSData['sample'] == f'{Year}_{Number}'
    # Loc = RUSData[F1 & F2].index
    # C11 = RUSData.loc[Loc,'RUS:C11'].values[0]
    # C33 = RUSData.loc[Loc,'RUS:C33'].values[0]
    # C44 = RUSData.loc[Loc,'RUS:C44'].values[0]
    # C66 = RUSData.loc[Loc,'RUS:C66'].values[0]
    # C13 = RUSData.loc[Loc,'RUS:C13'].values[0]
    # C12 = C11 - 2*C66

    # RUSStiffness = np.array([[C11, C12, C13,   0,   0,   0],
    #                             [C12, C11, C13,   0,   0,   0],
    #                             [C13, C13, C33,   0,   0,   0],
    #                             [  0,   0,   0, C44,   0,   0],
    #                             [  0,   0,   0,   0, C44,   0],
    #                             [  0,   0,   0,   0,   0, C66]])

    # # Investigate anisotropy
    # F_DA = [eValues[1] / eValues[0], eValues[2] / eValues[0], eValues[2] / eValues[1]]
    # I_DA = [IsoStiffness[1,1] / IsoStiffness[0,0], IsoStiffness[2,2] / IsoStiffness[0,0], IsoStiffness[2,2] / IsoStiffness[1,1]]
    # T_DA = [TransStiffness[1,1] / TransStiffness[0,0], TransStiffness[2,2] / TransStiffness[0,0], TransStiffness[2,2] / TransStiffness[1,1]]
    # E_DA = [RUSStiffness[1,1] / RUSStiffness[0,0], RUSStiffness[2,2] / RUSStiffness[0,0], RUSStiffness[2,2] / RUSStiffness[1,1]]
    # DA = ['$e_2$ / $e_1$', '$e_3$ / $e_1$', '$e_3$ / $e_2$']

    # Figure, Axis = plt.subplots(1,1)
    # Axis.plot(F_DA, color=(1,0,0), linestyle='none', marker='o', label='Fabric')
    # Axis.plot(I_DA, color=(1,0,1), linestyle='none', marker='o', label='Isotropic')
    # Axis.plot(T_DA, color=(0,0,1), linestyle='none', marker='o', label='Transverse')
    # Axis.plot(E_DA, color=(0,0,0), linestyle='none', marker='o', label='Experiment')
    # Axis.set_xticks(np.arange(3))
    # Axis.set_xticklabels(DA)
    # Axis.set_xlabel('Plane')
    # Axis.set_ylabel('Degree of anisotropy')
    # plt.legend()
    # plt.show(Figure)


if __name__ == '__main__':
    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add optional argument
    ScriptVersion = Parser.prog + ' version ' + __version__
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)

    # Read arguments from the command line
    Arguments = Parser.parse_args()
    Main()
