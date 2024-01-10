# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 15:14:21 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is used to compile the moment arm and line of action data calculated
    from individual specimen to summarise the group data. The settings are used
    to extract the data from the analysed specimen and conditions. This should be
    run from the 'Code' folder.     
    
"""

# %% Import packages

from matplotlib import pyplot as plt
import seaborn as sns
import os
from glob import glob
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
# import math
from skimage import io
# import re
# from shutil import copy2, rmtree

# %% Define functions

#Circle fitting function
def fitCircle(xPoints, yPoints):
    
    #Inputs:
        #xPoints - numpy array of x coordinates
        #yPoints - numpy array of y coordinates
    
    #Calculate means of XY points
    x_m = np.mean(xPoints)
    y_m = np.mean(yPoints)
    
    #Calculation of the reduced coordinates
    u = xPoints - x_m
    v = yPoints - y_m
    
    #Linear system defining the center (uc, vc) in reduced coordinates:
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)
    
    #Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = np.linalg.solve(A, B)
    
    #Get centres
    centreX = x_m + uc
    centreY = y_m + vc
    
    #Calculate distance to centre and radius
    Ri = np.sqrt((xPoints-centreX)**2 + (yPoints-centreY)**2)
    radius = np.mean(Ri)
    
    return centreX, centreY, radius

#Identifying 2D line intersections
def lineIntersect(Ax1, Ay1, Ax2, Ay2, Bx1, By1, Bx2, By2):
    #returns a (x, y) tuple or None if there is no intersection """
    d = (By2 - By1) * (Ax2 - Ax1) - (Bx2 - Bx1) * (Ay2 - Ay1)
    uA = ((Bx2 - Bx1) * (Ay1 - By1) - (By2 - By1) * (Ax1 - Bx1)) / d
    ###uB = ((Ax2 - Ax1) * (Ay1 - By1) - (Ay2 - Ay1) * (Ax1 - Bx1)) / d
    x = Ax1 + uA * (Ax2 - Ax1)
    y = Ay1 + uA * (Ay2 - Ay1)
 
    return x, y

#Procrustes function for identifying transform
def procrustes(X, Y, scaling = False, reflection = 'best'):
    
    """
    A port of MATLAB's `procrustes` function to Numpy.
    
    See: https://stackoverflow.com/questions/18925181/procrustes-analysis-with-numpy

    Procrustes analysis determines a linear transformation (translation,
    reflection, orthogonal rotation and scaling) of the points in Y to best
    conform them to the points in matrix X, using the sum of squared errors
    as the goodness of fit criterion.

        d, Z, [tform] = procrustes(X, Y)

    Inputs:
    ------------
    X, Y    
        matrices of target and input coordinates. they must have equal
        numbers of  points (rows), but Y may have fewer dimensions
        (columns) than X.

    scaling 
        if False, the scaling component of the transformation is forced
        to 1

    reflection
        if 'best' (default), the transformation solution may or may not
        include a reflection component, depending on which fits the data
        best. setting reflection to True or False forces a solution with
        reflection or no reflection respectively.

    Outputs
    ------------
    d       
        the residual sum of squared errors, normalized according to a
        measure of the scale of X, ((X - X.mean(0))**2).sum()

    Z
        the matrix of transformed Y-values

    tform   
        a dict specifying the rotation, translation and scaling that
        maps X --> Y

    """

    n,m = X.shape
    ny,my = Y.shape

    muX = X.mean(0)
    muY = Y.mean(0)

    X0 = X - muX
    Y0 = Y - muY

    ssX = (X0**2.).sum()
    ssY = (Y0**2.).sum()

    # centred Frobenius norm
    normX = np.sqrt(ssX)
    normY = np.sqrt(ssY)

    # scale to equal (unit) norm
    X0 /= normX
    Y0 /= normY

    if my < m:
        Y0 = np.concatenate((Y0, np.zeros(n, m-my)),0)

    # optimum rotation matrix of Y
    A = np.dot(X0.T, Y0)
    U,s,Vt = np.linalg.svd(A,full_matrices=False)
    V = Vt.T
    T = np.dot(V, U.T)

    if reflection != 'best':

        # does the current solution use a reflection?
        have_reflection = np.linalg.det(T) < 0

        # if that's not what was specified, force another reflection
        if reflection != have_reflection:
            V[:,-1] *= -1
            s[-1] *= -1
            T = np.dot(V, U.T)

    traceTA = s.sum()

    if scaling:

        # optimum scaling of Y
        b = traceTA * normX / normY

        # standarised distance between X and b*Y*T + c
        d = 1 - traceTA**2

        # transformed coords
        Z = normX*traceTA*np.dot(Y0, T) + muX

    else:
        b = 1
        d = 1 + ssY/ssX - 2 * traceTA * normY / normX
        Z = normY*np.dot(Y0, T) + muX

    # transformation matrix
    if my < m:
        T = T[:my,:]
    c = muX - b*np.dot(muY, T)
    
    #transformation values 
    tform = {'rotation':T, 'scale':b, 'translation':c}
   
    return d, Z, tform


# %% Parameters to change

#Provide list of specimen to extract data from
specimenNames = ['S1_r', 'S2_r', 'S3_r', 'S4_l', 'S5_l', 'S6_r', 'S7_r', 'S9_l']

#Provide list of conditions to extract data from
conditionNames = ['split50']
conditionLabels = ['50% Split']

#Provide list of conditions to extract data from
posNames = ['abd0', 'abd90', 'ABER', 'APP']
posLabels = ['0'+u'\u00b0'+' Abd.', 
             '90'+u'\u00b0'+' Abd.',
             'ABER', 'Apprehension']

#Provide list of loads to extract data from
loadNames = ['_0N', '_20N', '_40N']
loadLabels = ['0N', '20N', '40N']

#Provide list of planes to extract data from
planeNames = ['_SP', '_TP']
planeLabels = ['Scapular Plane', 'Transverse Plane']

#Provide list of muscle lines to extract data from
subscapNames = ['ss1', 'ss2', 'ss3', 'ss4']
subscapLabels = ['Superior', 'Middle-Superior', 'Middle-Inferior', 'Inferior']

# %% Set-up

#Set matplotlib parameters
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.weight'] = 'bold'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelweight'] = 'bold'
rcParams['legend.fontsize'] = 10
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['legend.framealpha'] = 0.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'

#Set grey color palette for three loads
greyPal = ['#e2e2e2', '#6a6a6a', '#000000']

#Set code directory to return to
homeDir = os.getcwd()

#Navigate to data directory
os.chdir('..\\Data\\')

#Set code to return to data directory
datDir = os.getcwd()

#Set dictionary to store data in
dataDict = {'specimen': [], 'condition': [], 'region': [],
            'position': [], 'load': [], 'plane': [],
            'lineOfAction': [], 'momentArm': [], 'stabilityRatio': []}

# %% Extract data

#Loop through specimen
for specimen in specimenNames:
    
    #Loop through conditions
    for condition in conditionNames:
        
        #Navigate to the condition directory
        os.chdir(f'{specimen}\\{condition}')
        
        #Get directory list that has processed data
        procDir = []
        for currDir in glob(os.getcwd()+'\\*\\'):
            #Navigate to directory
            os.chdir(currDir)
            #Check for processed data based on processed data files
            if os.path.isfile('lineOfActionData.csv') and os.path.isfile('momentArmData.csv'):
                #Append to list
                procDir.append(os.getcwd())
            #Return to upper directory
            os.chdir('..')

        #Loop through positions
        for pos in posNames:
            
            #Loop through loads
            for load in loadNames:
                
                #Loop through planes
                for plane in planeNames:
                    
                    #Try to find processed directory if available. Some specimen
                    #may have missing data for certain analysis conditions
                    try:
                        
                        #Identify the relevant directory to navigate into
                        currDir = [ii for ii in procDir if pos in ii and load in ii and plane in ii][0]
                                            
                        #Load into directory
                        os.chdir(currDir)
                        
                        #Read line of action data
                        loa = pd.read_csv('lineOfActionData.csv')
                        
                        #Read moment arm data
                        ma = pd.read_csv('momentArmData.csv')
                        
                        #Read in stability ratio data
                        sr = pd.read_csv('stabilityRatioData.csv')
                        
                        #Loop through muscle regions
                        for subscap in subscapNames:
                            
                            #Get loa value
                            loaVal = loa.loc[loa['subscapLine'] == subscap]['lineOfAction'].to_numpy()[0]
                            
                            #Get ma value
                            maVal = ma.loc[ma['subscapLine'] == subscap]['momentArm'].to_numpy()[0]
                            
                            #Get sr value
                            srVal = sr.loc[ma['subscapLine'] == subscap]['stabilityRatio'].to_numpy()[0]
                            
                            #Append values and info to data dictionary
                            dataDict['specimen'].append(specimen.split('_')[0])
                            dataDict['condition'].append(condition)
                            dataDict['region'].append(subscap)
                            dataDict['position'].append(pos)
                            dataDict['load'].append(load.split('_')[-1])
                            dataDict['plane'].append(plane.split('_')[-1])
                            dataDict['lineOfAction'].append(loaVal)
                            dataDict['momentArm'].append(maVal)
                            dataDict['stabilityRatio'].append(srVal)
                            
                    except:
                        
                        #Place nan values for current condition                        
                        print(f'NaNs inserted for condition {pos}{load}{plane} for {specimen}')
                        
                        #Loop through muscle regions
                        for subscap in subscapNames:
                            
                            #Append values and info to data dictionary
                            dataDict['specimen'].append(specimen.split('_')[0])
                            dataDict['condition'].append(condition)
                            dataDict['region'].append(subscap)
                            dataDict['position'].append(pos)
                            dataDict['load'].append(load.split('_')[-1])
                            dataDict['plane'].append(plane.split('_')[-1])
                            dataDict['lineOfAction'].append(np.nan)
                            dataDict['momentArm'].append(np.nan)
                            dataDict['stabilityRatio'].append(np.nan)
                            
                    #Return to specimen directory
                    os.chdir('..')
                    
        #Return to data directory
        os.chdir(datDir)
                        
#Convert to dataframe
groupData = pd.DataFrame.from_dict(dataDict)

#Export group data
groupData.to_csv(homeDir+'\\..\\Results\\Processed\\groupData.csv',
                 index = False)
                    
# %% Summarise data

#Group data by region, position, load and plane and describe
describedData = groupData.groupby(['region','position','load','plane']).describe().loc[['mean','sd']]

#Group and describe just loa data
describedData_loa = groupData.drop(['momentArm', 'stabilityRatio'],
                                   axis = 1).groupby(['region',
                                                      'position',
                                                      'load',
                                                      'plane']).describe(
                                                          ).loc[:,(slice(None),['mean','std'])]

#Export loa summary data
describedData_loa.to_csv('..\\Results\\Processed\\summaryLoA_groupData.csv',
                         index = True)

# %% Optimise curve fitting for line of action to stability ratio

#Function to calculate the cubic curve with constants
def cubic(x, x3, x2, x1, c):
    return x3*x**3 + x2*x**2 + x1*x + c

#Scapular plane

#Extract data for curve fitting in scapular plane
xFit = groupData.loc[groupData['plane'] == 'SP']['lineOfAction']
yFit = groupData.loc[groupData['plane'] == 'SP']['stabilityRatio']

#Fit the cubic function
pars_SP, cov_SP = curve_fit(f = cubic, xdata = xFit.dropna(), ydata = yFit.dropna(), #function & data
                            p0 = [0, 0, 0, 0], bounds = (-np.inf, np.inf)) #initial guess & bounds

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevsPars_SP = np.sqrt(np.diag(cov_SP))

#Calculate the residuals
fitRes_SP = yFit - cubic(xFit, *pars_SP)

# #Plot the fit data as an overlay on the scatter data
# plt.scatter(xFit, yFit, s = 20, color = '#00b3b3', label = 'Data')
# plt.plot(np.linspace(100, 220, 100),
#           cubic(np.linspace(100, 220, 100), *pars_SP),
#           linestyle = '--', linewidth = 2, color = 'black')

#Transverse plane

#Extract data for curve fitting in scapular plane
xFit = groupData.loc[groupData['plane'] == 'TP']['lineOfAction']
yFit = groupData.loc[groupData['plane'] == 'TP']['stabilityRatio']

#Fit the cubic function
pars_TP, cov_TP = curve_fit(f = cubic, xdata = xFit.dropna(), ydata = yFit.dropna(), #function & data
                            p0 = [0, 0, 0, 0], bounds = (-np.inf, np.inf)) #initial guess & bounds

# Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
stdevsPars_TP = np.sqrt(np.diag(cov_TP))

#Calculate the residuals
fitRes_TP = yFit - cubic(xFit, *pars_TP)

# #Plot the fit data as an overlay on the scatter data
# plt.scatter(xFit, yFit, s = 20, color = '#00b3b3', label = 'Data')
# plt.plot(np.linspace(130, 250, 100),
#           cubic(np.linspace(130, 250, 100), *pars_TP),
#           linestyle = '--', linewidth = 2, color = 'black')

# %% Visualise data - LoA

#Get max and min values in each plane to set y-axis limits
#Scapular plane
loaMax_SP = np.max(groupData.loc[groupData['plane'] == 'SP']['lineOfAction'])
loaMin_SP = np.min(groupData.loc[groupData['plane'] == 'SP']['lineOfAction'])
#Transverse plane
loaMax_TP = np.max(groupData.loc[groupData['plane'] == 'TP']['lineOfAction'])
loaMin_TP = np.min(groupData.loc[groupData['plane'] == 'TP']['lineOfAction'])

#Create figure
fig, ax = plt.subplots(nrows = 4, ncols = 2,
                       figsize = (10,10))

#Loop through regions for scapula plane
for subscap in subscapNames:
    
    #Create mean/SD point plot on relevant axis
    mn = sns.pointplot(data = groupData.loc[(groupData['plane'] == 'SP') &
                                            (groupData['region'] == subscap)],
                       x = 'position', y = 'lineOfAction', errorbar = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's', #size = 1,
                       join = False, dodge = 0.5, #zorder = 3,
                       ax = ax[subscapNames.index(subscap),0])
    
    # #Create boxplot on relevant axis
    # #Line of action data
    # bx = sns.boxplot(data = groupData.loc[(groupData['plane'] == 'SP') &
    #                                       (groupData['region'] == subscap)],
    #                  x = 'position', y = 'lineOfAction',
    #                  order = posNames, hue = 'load',
    #                  palette = 'Greys', width = 0.3, whis = [0,100],
    #                  ax = ax[subscapNames.index(subscap),0])
    
    # #Adjust colours of boxplot lines and fill
    # for ii in range(len(bx.artists)):
        
    #     #Get the current artist
    #     artist = bx.artists[ii]
        
    #     #Set the linecolor on the artist to the facecolor, and set the facecolor to None
    #     col = artist.get_facecolor()
    #     artist.set_edgecolor(col)
    #     artist.set_facecolor('None')
    
    #     #Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
    #     #Loop over them here, and use the same colour as above
    #     for jj in range(ii*6,ii*6+6):
    #         line = ax[subscapNames.index(subscap),0].lines[jj]
    #         line.set_color(col)
    #         line.set_mfc(col)
    #         line.set_mec(col)
    
    #Add point plot
    sp = sns.stripplot(data = groupData.loc[(groupData['plane'] == 'SP') &
                                          (groupData['region'] == subscap)],
                       x = 'position', y = 'lineOfAction',
                       order = posNames, hue = 'load',
                       palette = greyPal, 
                       edgecolor = 'black', #facecolor = 'white',
                       size = 3, linewidth = 1,
                       dodge = 0.5, alpha = 1, zorder = 4,
                       ax = ax[subscapNames.index(subscap),0])
    
    #Remove x-axis label
    ax[subscapNames.index(subscap),0].set_xlabel('')
    
    #Set x-axes tick labels if bottom row
    if subscapNames.index(subscap) == len(subscapNames)-1:
        #Set axis tick labels
        ax[subscapNames.index(subscap),0].set_xticklabels(posLabels,
                                                          fontsize = 10,
                                                          fontweight = 'bold',
                                                          rotation = 45,
                                                          ha = 'right')
    else:
        #Remove tick labels
        ax[subscapNames.index(subscap),0].set_xticklabels([])
        
    #Set title
    ax[subscapNames.index(subscap),0].set_title(subscapLabels[subscapNames.index(subscap)],
                                                fontsize = 12, fontweight = 'bold')
    
    #Set y-axis label
    ax[subscapNames.index(subscap),0].set_ylabel('Line of Action (\u00b0)\n(\u2190Inferior / Superior\u2192)',
                                                 fontsize = 10, fontweight = 'bold')
    
    #Add 'zero' line at 180 degrees
    ax[subscapNames.index(subscap),0].axhline(y = 180,
                                              color = 'lightgrey',
                                              linewidth = 1, linestyle = '--',
                                              zorder = 0)
    
    #Remove legend
    ax[subscapNames.index(subscap),0].get_legend().remove()
    
    #Set y-axis limit to max and minimum values in dataset for plane
    #Add a buffer of 2.5%
    ax[subscapNames.index(subscap),0].set_ylim([loaMin_SP - (loaMax_SP*0.025),
                                                loaMax_SP + (loaMax_SP*0.025)])
    
    # #Add twin axis for stability ratio values
    # ax2 = ax[subscapNames.index(subscap),0].twinx()
    
    # #Get y-axis limits of line of action axis and map these to stability ratio
    # #axis using fit cubic function
    # #Get the limits and set to the same on secondary axis
    # loaLimits = ax[subscapNames.index(subscap),0].get_ylim()
    # ax2.set_ylim(loaLimits)
    
    # #Calculate values for stability ratio y-ticks using function
    # srTicks = np.round(cubic(ax2.get_yticks(), *pars_SP), 1)
    # #Fix up -0 value if present
    # for ii, element in enumerate(srTicks):
    #     if element == -0.0:
    #         srTicks[ii] = 0.0
    
    # #Change tick labels on stability ratio axis
    # ax2.set_yticklabels(srTicks)
    
    # #Add y-label for stability ratio
    # ax2.set_ylabel('Stability Ratio',
    #                fontsize = 10, fontweight = 'bold')
    
    ### TODO: add images and labels on axes...
    
#Loop through regions for transverse plane
for subscap in subscapNames:
    
    #Create mean/SD point plot on relevant axis
    mn = sns.pointplot(data = groupData.loc[(groupData['plane'] == 'TP') &
                                            (groupData['region'] == subscap)],
                       x = 'position', y = 'lineOfAction', errorbar = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's',# markersize = 1,
                       join = False, dodge = 0.5, #zorder = 3,
                       ax = ax[subscapNames.index(subscap),1])
    
    # #Create boxplot on relevant axis
    # #Line of action data
    # bx2 = sns.boxplot(data = groupData.loc[(groupData['plane'] == 'TP') &
    #                                       (groupData['region'] == subscap)],
    #                  x = 'position', y = 'lineOfAction',
    #                  order = posNames, hue = 'load',
    #                  palette = ['black', 'red'], width = 0.4, whis = [0,100],
    #                  ax = ax[subscapNames.index(subscap),1])
    
    # #Adjust colours of boxplot lines and fill
    # for ii in range(len(bx2.artists)):
        
    #     #Get the current artist
    #     artist = bx2.artists[ii]
        
    #     #Set the linecolor on the artist to the facecolor, and set the facecolor to None
    #     col = artist.get_facecolor()
    #     artist.set_edgecolor(col)
    #     artist.set_facecolor('None')
    
    #     #Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
    #     #Loop over them here, and use the same colour as above
    #     for jj in range(ii*6,ii*6+6):
    #         line = ax[subscapNames.index(subscap),1].lines[jj]
    #         line.set_color(col)
    #         line.set_mfc(col)
    #         line.set_mec(col)
    
    #Add point plot
    sp = sns.stripplot(data = groupData.loc[(groupData['plane'] == 'TP') &
                                          (groupData['region'] == subscap)],
                       x = 'position', y = 'lineOfAction',
                       order = posNames, hue = 'load',
                       palette = greyPal, 
                       edgecolor = 'black', #facecolor = 'white',
                       size = 3, linewidth = 1,
                       dodge = 0.5, alpha = 1, zorder = 4,
                       ax = ax[subscapNames.index(subscap),1])
    
    #Remove x-axis label
    ax[subscapNames.index(subscap),1].set_xlabel('')
    
    #Set x-axes tick labels if bottom row
    if subscapNames.index(subscap) == len(subscapNames)-1:
        #Set axis tick labels
        ax[subscapNames.index(subscap),1].set_xticklabels(posLabels,
                                                          fontsize = 10,
                                                          fontweight = 'bold',
                                                          rotation = 45,
                                                          ha = 'right')
    else:
        #Remove tick labels
        ax[subscapNames.index(subscap),1].set_xticklabels([])
        
    #Set title
    ax[subscapNames.index(subscap),1].set_title(subscapLabels[subscapNames.index(subscap)],
                                                fontsize = 12, fontweight = 'bold')
    
    #Set y-axis label
    ax[subscapNames.index(subscap),1].set_ylabel('Line of Action (\u00b0)\n(\u2190Anterior / Posterior\u2192)',
                                                 fontsize = 10, fontweight = 'bold')
    
    #Add 'zero' line at 180 degrees
    ax[subscapNames.index(subscap),1].axhline(y = 180,
                                              color = 'lightgrey',
                                              linewidth = 1, linestyle = '--',
                                              zorder = 0)
    
    #Remove legend
    ax[subscapNames.index(subscap),1].get_legend().remove()
    
    #Set y-axis limit to max and minimum values in dataset for plane
    #Add a buffer of 2.5%
    ax[subscapNames.index(subscap),1].set_ylim([loaMin_TP - (loaMax_TP*0.025),
                                                loaMax_TP + (loaMax_TP*0.025)])
    
    #Set y-axis labels to act around 180
    ##### TODO: sort this out better...
    ax[subscapNames.index(subscap),1].set_yticks(np.array([150, 180, 210, 240]))
    
    # #Add twin axis for stability ratio values
    # ax2 = ax[subscapNames.index(subscap),1].twinx()
    
    # #Get y-axis limits of line of action axis and map these to stability ratio
    # #axis using fit cubic function
    # #Get the limits and set to the same on secondary axis
    # loaLimits = ax[subscapNames.index(subscap),1].get_ylim()
    # ax2.set_ylim(loaLimits)
    
    # #Match y-axis ticks
    # #### TODO: fix up less manual
    # ax2.set_yticks(np.array([150, 180, 210, 240]))
    
    # #Calculate values for stability ratio y-ticks using function
    # srTicks = np.round(cubic(ax2.get_yticks(), *pars_TP), 1)
    # #Fix up -0 value if present
    # for ii, element in enumerate(srTicks):
    #     if element == -0.0:
    #         srTicks[ii] = 0.0
    
    # #Change tick labels on stability ratio axis
    # ax2.set_yticklabels(srTicks)
    
    # #Add y-label for stability ratio
    # ax2.set_ylabel('Stability Ratio',
    #                fontsize = 10, fontweight = 'bold')
    
    ### TODO: add images and labels on axes...
    
#Tight layout
plt.tight_layout()

#Save figure
fig.savefig('..\\Results\\Figures\\loa_Fig_noStabilityRatio.png', dpi = 300, format = 'png')
fig.savefig('..\\Results\\Figures\\loa_Fig_noStabilityRatio.eps', dpi = 600, format = 'eps')

#Close figure
plt.close()

# %% Visualise data - moment arms

#Get max and min values in each plane to set y-axis limits
#Scapular plane
maMax_SP = np.max(groupData.loc[groupData['plane'] == 'SP']['momentArm'])
maMin_SP = np.min(groupData.loc[groupData['plane'] == 'SP']['momentArm'])
#Transverse plane
maMax_TP = np.max(groupData.loc[groupData['plane'] == 'TP']['momentArm'])
maMin_TP = np.min(groupData.loc[groupData['plane'] == 'TP']['momentArm'])

#Create figure
fig, ax = plt.subplots(nrows = 4, ncols = 2,
                       figsize = (10,10))

#Loop through regions for scapula plane
for subscap in subscapNames:
    
    #Create mean/SD point plot on relevant axis
    mn = sns.pointplot(data = groupData.loc[(groupData['plane'] == 'SP') &
                                            (groupData['region'] == subscap)],
                       x = 'position', y = 'momentArm', errorbar = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's', #markersize = 1,
                       join = False, dodge = 0.5, #zorder = 3,
                       ax = ax[subscapNames.index(subscap),0])
    
    #Add point plot
    sp = sns.stripplot(data = groupData.loc[(groupData['plane'] == 'SP') &
                                          (groupData['region'] == subscap)],
                       x = 'position', y = 'momentArm',
                       order = posNames, hue = 'load',
                       palette = greyPal, 
                       edgecolor = 'black', #facecolor = 'white',
                       size = 3, linewidth = 1,
                       dodge = 0.5, alpha = 1, zorder = 4,
                       ax = ax[subscapNames.index(subscap),0])
    
    #Remove x-axis label
    ax[subscapNames.index(subscap),0].set_xlabel('')
    
    #Set x-axes tick labels if bottom row
    if subscapNames.index(subscap) == len(subscapNames)-1:
        #Set axis tick labels
        ax[subscapNames.index(subscap),0].set_xticklabels(posLabels,
                                                          fontsize = 10,
                                                          fontweight = 'bold',
                                                          rotation = 45,
                                                          ha = 'right')
    else:
        #Remove tick labels
        ax[subscapNames.index(subscap),0].set_xticklabels([])
        
    #Set title
    ax[subscapNames.index(subscap),0].set_title(subscapLabels[subscapNames.index(subscap)],
                                                fontsize = 12, fontweight = 'bold')
    
    #Set y-axis label
    ax[subscapNames.index(subscap),0].set_ylabel('Moment Arm (mm)\n(\u2190Abduction / Adduction\u2192)',
                                                 fontsize = 10, fontweight = 'bold')
    
    #Add 'zero' line
    ax[subscapNames.index(subscap),0].axhline(y = 0,
                                              color = 'lightgrey',
                                              linewidth = 1, linestyle = '--',
                                              zorder = 0)
    
    #Remove legend
    ax[subscapNames.index(subscap),0].get_legend().remove()
    
    #Set y-axis limit to max and minimum values in dataset for plane
    #Add a buffer of 10%
    ax[subscapNames.index(subscap),0].set_ylim([maMin_SP - (maMax_SP*0.1),
                                                maMax_SP + (maMax_SP*0.1)])
    
    #Set y-axis labels appropriately
    ax[subscapNames.index(subscap),0].set_yticks(np.array([-30,-20,-10,0,10,20,30]))
    
#Loop through regions for transverse plane
for subscap in subscapNames:
    
    #Create mean/SD point plot on relevant axis
    mn = sns.pointplot(data = groupData.loc[(groupData['plane'] == 'TP') &
                                            (groupData['region'] == subscap)],
                       x = 'position', y = 'momentArm', errorbar = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's', #markersize = 1,
                       join = False, dodge = 0.5, #zorder = 3,
                       ax = ax[subscapNames.index(subscap),1])
    
    #Add point plot
    sp = sns.stripplot(data = groupData.loc[(groupData['plane'] == 'TP') &
                                          (groupData['region'] == subscap)],
                       x = 'position', y = 'momentArm',
                       order = posNames, hue = 'load',
                       palette = greyPal, 
                       edgecolor = 'black', #facecolor = 'white',
                       size = 3, linewidth = 1,
                       dodge = 0.5, alpha = 1, zorder = 4,
                       ax = ax[subscapNames.index(subscap),1])
    
    #Remove x-axis label
    ax[subscapNames.index(subscap),1].set_xlabel('')
    
    #Set x-axes tick labels if bottom row
    if subscapNames.index(subscap) == len(subscapNames)-1:
        #Set axis tick labels
        ax[subscapNames.index(subscap),1].set_xticklabels(posLabels,
                                                          fontsize = 10,
                                                          fontweight = 'bold',
                                                          rotation = 45,
                                                          ha = 'right')
    else:
        #Remove tick labels
        ax[subscapNames.index(subscap),1].set_xticklabels([])
        
    #Set title
    ax[subscapNames.index(subscap),1].set_title(subscapLabels[subscapNames.index(subscap)],
                                                fontsize = 12, fontweight = 'bold')
    
    #Set y-axis label
    ax[subscapNames.index(subscap),1].set_ylabel('Moment Arm (mm)\n(\u2190Int. Rot. / Ext. Rot.\u2192)',
                                                 fontsize = 10, fontweight = 'bold')
    
    #Add 'zero' line
    ax[subscapNames.index(subscap),1].axhline(y = 0,
                                              color = 'lightgrey',
                                              linewidth = 1, linestyle = '--',
                                              zorder = 0)
    
    #Remove legend
    ax[subscapNames.index(subscap),1].get_legend().remove()
    
    #Set y-axis limit to max and minimum values in dataset for plane
    #Add a buffer of 10%
    ax[subscapNames.index(subscap),1].set_ylim([maMin_TP - (maMax_TP*0.1),
                                                maMax_TP + (maMax_TP*0.1)])
    
#Tight layout
plt.tight_layout()

#Save figure
fig.savefig('..\\Results\\Figures\\ma_Fig.png', dpi = 300, format = 'png')
fig.savefig('..\\Results\\Figures\\ma_Fig.eps', dpi = 600, format = 'eps')

#Close figure
plt.close()

# %% Visualise methodological processes on radiograph

#Uses S4 90 abduction, 40N, scapular plane image, ss4

#Create figure
fig, ax = plt.subplots(figsize = (20,6), nrows = 1, ncols = 4)

#Load image
img = io.imread('S4_l\\split50\\abd90_40N_SP_D0302_T1333\\abd90_40N_SP_D0302_T1333_1_proc.tif')

#Turn off axes and display image in each
for currAx in ax.flatten():
    #Turn off axes
    currAx.axis('off')
    #Display image
    currAx.imshow(img, cmap = 'gray', origin = 'upper', vmin = 0, vmax = 255)
    
#Set tight layout
plt.tight_layout()

#Get the digitised points data for this reference image
#Subscap points
ssPts = pd.read_csv('S4_l\\split50\\abd90_40N_SP_D0302_T1333\\ss4.csv')
ssPts_x = ssPts['X'].to_numpy()
ssPts_y = ssPts['Y'].to_numpy()
#Humeral head points
hhPts = pd.read_csv('S4_l\\split50\\abd90_40N_SP_D0302_T1333\\hh.csv')
hhPts_x = hhPts['X'].to_numpy()
hhPts_y = hhPts['Y'].to_numpy()

#Glenoid plane points
#Need to calculate via reference points
#Estimate glenoid plane based on image transformation

#Import the digitised reference points from the current image
transPts = pd.read_csv('S4_l\\split50\\abd90_40N_SP_D0302_T1333\\pts.csv')
transPts = np.transpose(np.array([transPts['X'].to_numpy(), transPts['Y'].to_numpy()]))

#Get original digitised points and glenoid pts
refPts = pd.read_csv('S4_l\\split50\\abd0_0N_SP_D0302_T1308\\pts.csv')
refPts = np.transpose(np.array([refPts['X'].to_numpy(), refPts['Y'].to_numpy()]))

#Calculate procrustes transformation between points
resErr, newPts, tform = procrustes(refPts, transPts)

#Get the original gp points
gpPts = pd.read_csv('S4_l\\split50\\abd0_0N_SP_D0302_T1308\\gp.csv')
gpPts = np.transpose(np.array([gpPts['X'].to_numpy(), gpPts['Y'].to_numpy()]))

#Apply the transformation to the reference glenoid plane points
#Transform reference glenoid points to new array
gpPts_trans = np.zeros(gpPts.shape)
for pp in range(gpPts.shape[0]):
    gpPts_trans[pp] = tform['rotation'].dot(gpPts[pp]) - tform['translation']
    
#Convert to same format
gpPts_x = gpPts_trans[:,0]
gpPts_y = gpPts_trans[:,1]
 
#Display digitised points on first image
#Subscap points
ax.flatten()[0].scatter(ssPts_x, ssPts_y, s = 20, c = 'red', zorder = 5)
#Glenoid plane points
ax.flatten()[0].scatter(gpPts_x, gpPts_y, s = 20, c = 'yellow', zorder = 4)
#Humeral head points
ax.flatten()[0].scatter(hhPts_x, hhPts_y, s = 20, c = 'cyan', zorder = 3)
#Add title
ax.flatten()[0].set_title('a) Digitised Points', fontsize = 16, fontweight = 'bold')

#Display fits to these points on second image
#Subscap
#Fit line
m,c = np.polyfit(ssPts_x, ssPts_y, 1)
#Plot fitted line
ax.flatten()[1].plot(ssPts_x, (m * ssPts_x + c), c = 'red', lw = 2, zorder = 5)
#Glenoid plane
#Fit line
m,c = np.polyfit(gpPts_x, gpPts_y, 1)
#Plot fitted line
ax.flatten()[1].plot(gpPts_x, (m * gpPts_x + c), c = 'yellow', lw = 2, zorder = 4)
#Humeral head
#Fit circle
hhCentreX, hhCentreY, hhRadius = fitCircle(hhPts_x, hhPts_y)
#Plot cirlcle
ax.flatten()[1].add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                                      edgecolor = 'cyan', facecolor = 'none', lw = 2,
                                      zorder = 3))
ax.flatten()[1].scatter(hhCentreX, hhCentreY, s = 20, c = 'cyan', zorder = 3)
#Add title
ax.flatten()[1].set_title('b) Geometrical Fits', fontsize = 16, fontweight = 'bold')

#Display line of action calculation on third image
#Subscap
#Fit line
subscapM, subscapC = np.polyfit(ssPts_x, ssPts_y, 1)
#Plot fitted line
ax.flatten()[2].plot(ssPts_x, (subscapM * ssPts_x + subscapC), c = 'red', lw = 2, zorder = 5)
#Calculate intersection of subscap to glenoid plane
glenoidM, glenoidC = np.polyfit(gpPts_x, gpPts_y, 1)
#Extract two points on a line fit to the glenoid plane
gx1 = 0
gx2 = 100
gy1 = glenoidM * gx1 + glenoidC #y-intercept; x = 0
gy2 = glenoidM * gx2 + glenoidC #x = 100
#Extract two points on the subscapularis linefit
sx1 = 0
sx2 = 100
sy1 = subscapM * sx1 + subscapC #y-intercept; x = 0
sy2 = subscapM * sx2 + subscapC #x = 100
#Calculate line intersection
intX, intY = lineIntersect(gx1, gy1, gx2, gy2, sx1, sy1, sx2, sy2)
ax.flatten()[2].scatter(intX, intY, s = 20, c = 'yellow', zorder = 6)
#Extend subscapularis line
#Line
ax.flatten()[2].plot(np.array(((ax.flatten()[2].get_xlim()[1] / 3), np.min(ssPts_x))),
          subscapM * np.array(((ax.flatten()[2].get_xlim()[1] / 3), np.min(ssPts_x))) + subscapC,
          c = 'red', lw = 2, ls = '--', zorder = 5)
#Arrow head
ax.flatten()[2].arrow(np.max(ssPts_x), subscapM * np.max(ssPts_x) + subscapC,
                      (ax.flatten()[2].get_xlim()[1] / 3) - np.max(ssPts_x),
                      (subscapM * (ax.flatten()[2].get_xlim()[1] / 3) + subscapC) - (subscapM * np.max(ssPts_x) + subscapC),
                      head_width = 30, head_length = 30, 
                      lw = 0, color = 'red', ls = ':', edgecolor = 'none',
                      length_includes_head = False)
#Plot the glenoid plane
#Line
ax.flatten()[2].plot(gpPts_x, (glenoidM * gpPts_x + glenoidC), c = 'yellow', lw = 2, zorder = 4)
#Arrow head
ax.flatten()[2].arrow(gpPts_x[0], glenoidM * gpPts_x[0] + glenoidC,
                      gpPts_x[1] - gpPts_x[0],
                      gpPts_x[1] - (glenoidM * gpPts_x[0] + glenoidC),
                      head_width = 30, head_length = 30, 
                      lw = 0, color = 'yellow', ls = ':', edgecolor = 'none',
                      length_includes_head = True)
#Extend the glenoid plane to the intersection point
ax.flatten()[2].plot(np.array(((intX), np.min(gpPts_x))),
          glenoidM * np.array(((intX), np.min(gpPts_x))) + glenoidC,
          c = 'yellow', lw = 2, ls = '--', zorder = 4)
#Plot the perpendicular line
#Get gradient
mPerp = 1 / (-glenoidM)
#Work out y-intercept
cPerp = ((mPerp*intX) - intY) / -1
#Plot line
ax.flatten()[2].plot(np.array(((intX), intX-250)),
          mPerp * np.array(((intX), intX-250)) + cPerp,
          c = 'lime', lw = 2, zorder = 3)
#Arrow head
ax.flatten()[2].arrow(intX, intY,
                      -250, 
                      (mPerp * (intX-250) + cPerp) - intY,
                      head_width = 30, head_length = 30, 
                      lw = 0, color = 'lime', ls = ':', edgecolor = 'none',
                      length_includes_head = False)
#Add title
ax.flatten()[2].set_title('c) Line of Action Calculation', fontsize = 16, fontweight = 'bold')

#Display moment arm calculation on fourth image
#Subscap
#Fit line
subscapM, subscapC = np.polyfit(ssPts_x, ssPts_y, 1)
#Plot fitted line
ax.flatten()[3].plot(ssPts_x, (subscapM * ssPts_x + subscapC), c = 'red', lw = 2, zorder = 5)
#Extend subscap line
ax.flatten()[3].plot(np.array(((ax.flatten()[2].get_xlim()[0] + (np.max(ssPts_x)/4)), np.max(ssPts_x))),
          subscapM * np.array(((ax.flatten()[2].get_xlim()[0] + (np.max(ssPts_x)/4)), np.max(ssPts_x))) + subscapC,
          c = 'red', lw = 2, ls = '--', zorder = 5)
#Plot humeral head centre
ax.flatten()[3].scatter(hhCentreX, hhCentreY, s = 20, c = 'cyan', zorder = 4)
#Calculate the point on the extended subscapularis line that is intercepted
#by a perpendicular line from the humerus
#Set the points at the start and the end of the subscapularis line
x1 = ax.flatten()[2].get_xlim()[0] + (np.max(ssPts_x)/4)
y1 = subscapM * x1 + subscapC
x2 = np.max(ssPts_x)
y2 = subscapM * x2 + subscapC
#Set the humeral head point
x3 = hhCentreX
y3 = hhCentreY
#Calculate intersecting point with perpendicular line
#SEE: https://stackoverflow.com/questions/1811549/perpendicular-on-a-line-from-a-given-point
k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)**2 + (x2-x1)**2)
x4 = x3 - k * (y2-y1)
y4 = y3 + k * (x2-x1)
#Plot the perpendicular line between the points
#This is effectively the moment arm - but only use here as visualisation
ax.flatten()[3].plot(np.array((x3,x4)),np.array((y3,y4)),
                     c = 'cyan', lw = 2, zorder = 3)
#Add title
ax.flatten()[3].set_title('d) Moment Arm Calculation', fontsize = 16, fontweight = 'bold')

#Save figure
fig.savefig('..\\Results\\Figures\\methods_Fig.png', dpi = 300, format = 'png')

#Close figure
plt.close()

# %% OLD TEST CODE BELOW...





# %%

# %% Displaysome basic mean/SD parameters

# #Each region for X abduction in X plane with X load
# currPos = 'abd90'
# currPlane = 'TP'
# currLoad = '40N'
# for subscap in subscapNames:
#     mu = np.median(groupData.loc[(groupData['region'] == subscap) &
#                                (groupData['plane'] == currPlane) &
#                                (groupData['position'] == currPos) &
#                                (groupData['load'] == currLoad),['lineOfAction']].to_numpy())
#     sigma25,sigma75 = np.percentile(groupData.loc[(groupData['region'] == subscap) &
#                                                      (groupData['plane'] == currPlane) &
#                                                      (groupData['position'] == currPos) &
#                                                      (groupData['load'] == currLoad),['lineOfAction']].to_numpy(), [25,75])
#     print(f'{currPlane}.{subscap}_{currPos}_{currLoad} = [{mu},{sigma25},{sigma75}];')

