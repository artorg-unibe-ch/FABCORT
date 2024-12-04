#%% !/usr/bin/env python3

Description = """
Investigates the relationships between resolution,
Otsu threshold for segmentation, and resulting BVTV
"""

__author__ = ['Mathieu Simon']
__date_created__ = '04-12-2024'
__date__ = '04-12-2024'
__license__ = 'GPL'
__version__ = '1.0'


#%% Imports

import argparse
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
from scipy.stats import t
from matplotlib.cm import jet
import matplotlib.pyplot as plt

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

#%% Functions

def Resample(Image, Factor=None, Size=[None], Spacing=[None], Order=1):

    """
    Resample a SimpleITK image by either a given factor, a new size, or
    a new voxel spacing. Order stands for interpolation order e.g.
    Order = 1: Linear interpolation 
    """

    Dimension = Image.GetDimension()
    OriginalSpacing = np.array(Image.GetSpacing())
    OriginalSize = np.array(Image.GetSize())
    PhysicalSize = OriginalSize * OriginalSpacing

    Origin = Image.GetOrigin()
    Direction = Image.GetDirection()

    if Factor:
        NewSize = [round(Size/Factor) for Size in Image.GetSize()] 
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Size[0]:
        NewSize = Size
        NewSpacing = [PSize/Size for Size, PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing[0]:
        NewSpacing = Spacing
        NewSize = [np.floor(Size/Spacing).astype('int') + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]

    NewArray = np.zeros(NewSize[::-1],'int')
    NewImage = sitk.GetImageFromArray(NewArray)
    NewImage.SetOrigin(Origin - OriginalSpacing/2)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    Resampled = sitk.Resample(Image, NewImage, Transform, Order)
    
    return Resampled

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
    # SMax = 1e4
    # SMin = 1e-3

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

    DataPath = Path(__file__).parents[1] / '00_Data'
    Files = sorted([F for F in Path.iterdir(DataPath) if F.name.endswith('.mhd')])


    # Investigate different resolution factors
    Factors = [1,2,4]
    O1, O2, O4 = [], [], []
    BVTV1, BVTV2, BVTV4 = [], [], []
    Otsus = [O1, O2, O4]
    BVTVs = [BVTV1, BVTV2, BVTV4]
    for File in Files:

        # Read file
        Image = sitk.ReadImage(str(File))

        for F, O, BVTV in zip(Factors, Otsus, BVTVs):
            
            # Resample
            Resampled = Resample(Image, Factor=F)

            # Compute Otsu's threshold
            Otsu = sitk.OtsuThresholdImageFilter()
            Otsu.SetInsideValue(0)
            Otsu.SetOutsideValue(1)
            Seg = Otsu.Execute(Resampled)
            Thresold = Otsu.GetThreshold()

            # Compute bone volume fraction
            Array = sitk.GetArrayFromImage(Seg)
            VolumeFraction = Array.sum() / Array.size

            # Store values
            O.append(Thresold)
            BVTV.append(VolumeFraction)

    # Group by file instead than by factors
    Otsus = np.array(Otsus).T
    BVTVs = np.array(BVTVs).T

    # Sort by bone volume fraction for plotting
    Args = np.argsort(BVTVs[:,0])

    # Plot
    Figure, Axis = plt.subplots(1,2, figsize=(11,5), dpi=200)
    for i, Arg in enumerate(Args):
        Color = jet(i/(len(Args)-1))
        Axis[0].plot(Factors, Otsus[Args[i]], color=Color)
        Axis[1].plot(Factors, BVTVs[Args[i]], color=Color)
    for i in range(2):
        Axis[i].set_xlabel('Downscale factor (-)')
        Axis[i].set_xticks(Factors)
    Axis[0].set_ylabel('Otsu\'s threshold (-)')
    Axis[1].set_ylabel('Bone volume fraction (-)')
    plt.show(Figure)

    # Iteratively change threshold to reach same BV/TV as initial resolution
    BVTV_N2 = []
    BVTV_N4 = []
    Otsu_N2 = []
    Otsu_N4 = []
    BVTV_N = [BVTV_N2, BVTV_N4]
    Otsu_N = [Otsu_N2, Otsu_N4]
    for f, File in enumerate(Files):
        
        # Read file
        Image = sitk.ReadImage(str(File))

        for i, (F, B, O) in enumerate(zip(Factors[1:], BVTV_N, Otsu_N)):

            # Resample
            Resampled = Resample(Image, Factor=F)

            # Iteratively decrease the threshold to correct BV/TV
            Threshold = Otsus[f,i+1]
            Difference = abs(BVTVs[f,0] - BVTVs[f,i+1])

            # Store initial values
            BestThreshold = Threshold
            BestVolumeFraction = BVTVs[f,i+1]

            for i in range(1,int(Threshold)):

                # Threshold image
                Binarize = sitk.BinaryThresholdImageFilter()
                Binarize.SetUpperThreshold(Threshold-i)
                Binarize.SetOutsideValue(1)
                Binarize.SetInsideValue(0)
                Bin = Binarize.Execute(Resampled)

                # Compute new bone volume fraction
                Array = sitk.GetArrayFromImage(Bin)
                VolumeFraction = Array.sum() / Array.size

                # Compute new difference
                NewDifference = abs(BVTVs[f,0] - VolumeFraction)

                # If difference decreases, store it otherwise stop
                if NewDifference < Difference:
                    Difference = NewDifference
                    BestThreshold = Threshold-i
                    BestVolumeFraction = VolumeFraction
                else:
                    # Store new threshold and bone volume fraction
                    O.append(BestThreshold)
                    B.append(BestVolumeFraction)
                    break
    
    # Group by file instead than by factors
    Otsu_N = np.array(Otsu_N).T
    BVTV_N = np.array(BVTV_N).T

    # Concatenate with original resolution values
    Otsu_N = np.vstack([Otsus[:,0], Otsu_N.T]).T
    BVTV_N = np.vstack([BVTVs[:,0], BVTV_N.T]).T

    # Plot
    Figure, Axis = plt.subplots(1,2, figsize=(11,5), dpi=200)
    for i, Arg in enumerate(Args):
        Color = jet(i/(len(Args)-1))
        Axis[0].plot(Factors, Otsu_N[Arg], color=Color)
        Axis[1].plot(Factors, BVTV_N[Arg], color=Color)
    for i in range(2):
        Axis[i].set_xlabel('Downscale factor (-)')
        Axis[i].set_xticks(Factors)
    Axis[0].set_ylabel('Otsu\'s threshold (-)')
    Axis[1].set_ylabel('Bone volume fraction (-)')
    plt.show(Figure)

    # Save thresholds for fabric analysis with Medtool
    with open(Path(__file__).parent / 'Thresholds.csv','w') as F:
        F.write('$Sample,$Factor,$Threshold,\n')
        for i, File in enumerate(Files):
            for Factor, Threshold in zip(Factors[1:], Otsu_N[i,1:]):
                Line = File.name[:-4] + f',{Factor},' + f'{int(Threshold)},\n'
                F.write(Line)



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
