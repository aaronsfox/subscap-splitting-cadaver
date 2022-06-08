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
# from skimage import io
# import re
# from shutil import copy2, rmtree

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

# %% Visualise data

##### TODO: this is currently done for minimal dataset --- needs to be optimised
##### to include all data/better formatted

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
                       x = 'position', y = 'lineOfAction', ci = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's', markersize = 1,
                       join = False, dodge = 0.5, zorder = 3,
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
    
    #Add twin axis for stability ratio values
    ax2 = ax[subscapNames.index(subscap),0].twinx()
    
    #Get y-axis limits of line of action axis and map these to stability ratio
    #axis using fit cubic function
    #Get the limits and set to the same on secondary axis
    loaLimits = ax[subscapNames.index(subscap),0].get_ylim()
    ax2.set_ylim(loaLimits)
    
    #Calculate values for stability ratio y-ticks using function
    srTicks = np.round(cubic(ax2.get_yticks(), *pars_SP), 1)
    #Fix up -0 value if present
    for ii, element in enumerate(srTicks):
        if element == -0.0:
            srTicks[ii] = 0.0
    
    #Change tick labels on stability ratio axis
    ax2.set_yticklabels(srTicks)
    
    #Add y-label for stability ratio
    ax2.set_ylabel('Stability Ratio',
                   fontsize = 10, fontweight = 'bold')
    
    ### TODO: add images and labels on axes...
    
#Loop through regions for transverse plane
for subscap in subscapNames:
    
    #Create mean/SD point plot on relevant axis
    mn = sns.pointplot(data = groupData.loc[(groupData['plane'] == 'TP') &
                                            (groupData['region'] == subscap)],
                       x = 'position', y = 'lineOfAction', ci = 'sd',
                       order = posNames, hue = 'load',
                       palette = greyPal, markers = 's', markersize = 1,
                       join = False, dodge = 0.5, zorder = 3,
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
    
    #Add twin axis for stability ratio values
    ax2 = ax[subscapNames.index(subscap),1].twinx()
    
    #Get y-axis limits of line of action axis and map these to stability ratio
    #axis using fit cubic function
    #Get the limits and set to the same on secondary axis
    loaLimits = ax[subscapNames.index(subscap),1].get_ylim()
    ax2.set_ylim(loaLimits)
    
    #Match y-axis ticks
    #### TODO: fix up less manual
    ax2.set_yticks(np.array([150, 180, 210, 240]))
    
    #Calculate values for stability ratio y-ticks using function
    srTicks = np.round(cubic(ax2.get_yticks(), *pars_TP), 1)
    #Fix up -0 value if present
    for ii, element in enumerate(srTicks):
        if element == -0.0:
            srTicks[ii] = 0.0
    
    #Change tick labels on stability ratio axis
    ax2.set_yticklabels(srTicks)
    
    #Add y-label for stability ratio
    ax2.set_ylabel('Stability Ratio',
                   fontsize = 10, fontweight = 'bold')
    
    ### TODO: add images and labels on axes...
    
#Tight layout
plt.tight_layout()

#Save figure
fig.savefig('..\\Results\\Figures\\loa_Fig.png', dpi = 300, format = 'png')

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

