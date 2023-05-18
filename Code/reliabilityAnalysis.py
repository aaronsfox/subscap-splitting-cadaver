# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 15:14:21 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Code used for allocating reliability conditions and analysing these.
    
"""

# %% Import packages

import random
import pandas as pd
from glob import glob
import os
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import math
import pingouin as pg

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

#Set specimen
specimenNames = ['S1_r', 'S2_r', 'S3_r', 'S4_l', 'S5_l', 'S6_r', 'S7_r', 'S8_l', 'S9_l']

#Set conditions
conditionNames = ['split50', 'split25_upper', 'split25_lower']

#Set loads
loadNames = ['0N', '20N', '40N']

#Set muscle lines
muscleNames = ['ss1', 'ss2', 'ss3', 'ss4']

#Set arm positions
posNames = ['abd0', 'abd90', 'ABER', 'APP']

#Set phantom diameter for scaling
phantomDiameter_SP = 6.5
phantomDiameter_TP = 6.5

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

# %% Identify the conditions to assess

"""

Here we allocate a load and muscle for each specimen and condition randomly. The 
selected option will be assessed in both the transverse and scapula plane.

"""

#Set dicitonary to store data in
reliabilityConditions = {'specimen': [], 'condition': [],
                         'load': [], 'muscle': [], 'position': []}

#Loop through specimen
for specimen in specimenNames:
    
    #Loop through conditions
    for condition in conditionNames:
        
        #Set seed for random selection of load
        random.seed(specimenNames.index(specimen) + conditionNames.index(condition) + 12345)
        
        #Select load
        loadChoice = random.choice(loadNames)
        
        #Set seed for random selection of muscle
        random.seed(specimenNames.index(specimen) + conditionNames.index(condition) + 54321)
        
        #Select muscle
        muscleChoice = random.choice(muscleNames)
        
        #Set seed for random selection of position
        random.seed(specimenNames.index(specimen) + conditionNames.index(condition) + 111999)
        
        #Select muscle
        posChoice = random.choice(posNames)
        
        #Append to dictionary
        reliabilityConditions['specimen'].append(specimen)
        reliabilityConditions['condition'].append(condition)
        reliabilityConditions['load'].append(loadChoice)
        reliabilityConditions['muscle'].append(muscleChoice)
        reliabilityConditions['position'].append(posChoice)
        
#Convert to dataframe
reliabilityConditions_df = pd.DataFrame.from_dict(reliabilityConditions)

#Export to file
reliabilityConditions_df.to_excel('..\\Results\\Reliability\\reliabilitySelectConditions.xlsx',
                                  index = False)

# %% Extract data for reliability calculations

#Set lists to store data to
lineOfAction1 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}
momentArm1 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}
stabilityRatio1 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}
lineOfAction2 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}
momentArm2 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}
stabilityRatio2 = {'specimen': [], 'condition': [], 'load': [], 'muscle': [], 'position': [], 'plane': [], 'val': []}

#Remove some specific conditions
#Specimen 1 split 25 upper has no phantom (index 1)
#Specimen 2 split 50 phantom not appropriate (index 3)
#Specimen 6 split 50 has no phantom (index 15)
#Specimen 6 split 25 upper has no phantom (index 16)
#Specimen 8 not collected at all (index 21, 22, 23)
reliabilityConditions_df.drop([reliabilityConditions_df.index[1],
                               reliabilityConditions_df.index[3],
                               reliabilityConditions_df.index[15],
                               reliabilityConditions_df.index[16],
                               reliabilityConditions_df.index[21],
                               reliabilityConditions_df.index[22],
                               reliabilityConditions_df.index[23]],
                              inplace = True)
reliabilityConditions_df = reliabilityConditions_df.reset_index(drop = True)

#Loop through condition indices
for ii in reliabilityConditions_df.index:
    
    #Get characteristics for current condition
    specimen = reliabilityConditions_df.iloc[ii]['specimen']
    condition = reliabilityConditions_df.iloc[ii]['condition']
    position = reliabilityConditions_df.iloc[ii]['position']
    load = reliabilityConditions_df.iloc[ii]['load']
    muscle = reliabilityConditions_df.iloc[ii]['muscle']
    
    #List all folders in specimen and condition that have the load and position
    dirList = glob(f'..\\Data\\{specimen}\\{condition}\\{position}_{load}_*\\')
    
    #Get scapula and transverse plane directories
    boolDirSP = ['_SP_' in currDir for currDir in dirList]
    boolDirTP = ['_TP_' in currDir for currDir in dirList]
    dirSP = [xx for ii, xx in enumerate(dirList) if boolDirSP[ii]]
    dirTP = [xx for ii, xx in enumerate(dirList) if boolDirTP[ii]]
    allDir = dirSP + dirTP
    
    #Loop through directory list and run calculations if data present
    for currDir in allDir:
        
        #Check for muscle file
        if os.path.exists(f'{currDir}{muscle}_2.csv'):
            
            #Identify plane for current data
            if '_SP' in currDir:
                currPlane = 'SP'
            else:
                currPlane = 'TP'
            
            #Extract original data
            loaOrig = pd.read_csv(f'{currDir}lineOfActionData.csv')
            maOrig = pd.read_csv(f'{currDir}momentArmData.csv')
            srOrig = pd.read_csv(f'{currDir}stabilityRatioData.csv')
        
            #Read in the reliability data
            ssPoints = pd.read_csv(f'{currDir}{muscle}_2.csv')
            gpPoints = pd.read_csv(f'{currDir}gp_2.csv')
            hhPoints = pd.read_csv(f'{currDir}hh_2.csv')
            phantomPoints = pd.read_csv(f'{currDir}phantom.csv')
            
            #Identify the processed .tif file
            try:
                imageFile = glob(f'{currDir}*_proc.tif')[0]
            except:
                imageFile = glob(f'{currDir}*.tif')[0]
                
            #Load image
            img = io.imread(imageFile)
            
            #Fit circle to new humeral head
            hhCentreX, hhCentreY, hhRadius = fitCircle(hhPoints['X'].to_numpy(),
                                                       hhPoints['Y'].to_numpy())
            
            #Fit circle to new phantom bead
            _, _, phantomRadius = fitCircle(phantomPoints['X'].to_numpy(),
                                            phantomPoints['Y'].to_numpy())
            
            # %% Moment arm calculation
            
            #Set point to extend subscapularis line to based on humeral head
            #position (i.e. halfway between humeral head and edge of image)
            if hhCentreX < np.mean(gpPoints['X']):
                extendX = hhCentreX / 2
            else:
                extendX = hhCentreX + ((img.shape[1] - hhCentreX) / 2)
                
            #Set humeral head point as p3 for moment arm calculations
            p3 = np.asarray((hhCentreX, hhCentreY))
            
            #Create axes to visualise on
            fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)

            #Display image on current axes
            ax.imshow(img, cmap = 'gray', origin = 'upper')
                
            #Turn off axes
            ax.axis('off')
                    
            #Display humeral head points
            ax.scatter(hhPoints['X'], hhPoints['Y'], s = 10, c = 'green')
            ax.add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                                     edgecolor = 'green', facecolor = 'none'))
            ax.scatter(hhCentreX, hhCentreY, s = 500, c = 'green', marker = '+')
            
            #Display subscapularis line
            #Plot points
            ax.scatter(ssPoints['X'],ssPoints['Y'], s = 10, c = 'red')
            #Fit line
            m,c = np.polyfit(ssPoints['X'], ssPoints['Y'], 1)
            #Plot fitted line
            ax.plot(ssPoints['X'], m * ssPoints['X'] + c, c = 'red', lw = 1)
            #Extend line of action
            #Check which way to extend and whether to use min or max
            if extendX < hhCentreX:
                ax.plot(np.array((extendX,np.min(ssPoints['X']))), m * np.array((extendX,np.min(ssPoints['X']))) + c,
                        c = 'red', lw = 1, ls = '--')
            else:
                ax.plot(np.array((np.max(ssPoints['X']),extendX)), m * np.array((np.max(ssPoints['X']),extendX)) + c,
                        c = 'red', lw = 1, ls = '--')
                
            #Extract two points on the line for later moment arm calculation
            #Take the end points of the line on the x-axes
            p1 = np.asarray((ax.get_xlim()[0],
                             m * ax.get_xlim()[0] + c))
            p2 = np.asarray((ax.get_xlim()[1],
                             m * ax.get_xlim()[1] + c))
            
            #Calculate the point on the extended subscapularis line that is intercepted
            #by a perpendicular line from the humerus
            
            #Set the points at the start and the end of the subscapularis line
            x1 = hhCentreX/2
            y1 = m * (hhCentreX/2)+ c
            x2 = max(ssPoints['X'])
            y2 = m * x2 + c
            
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
            ax.plot(np.array((x3,x4)),np.array((y3,y4)),
                    c = 'yellow', lw = 1, ls = '--')
            
            #Calculate moment arm distance
            #Includes scale by bead size
            #Use cross product to calculate moment arm
            if currPlane == 'SP':
                ma = np.cross(p2-p1,p3-p1) / np.linalg.norm(p2-p1) / (phantomRadius / (phantomDiameter_SP/2))
            else:
                ma = np.cross(p2-p1,p3-p1) / np.linalg.norm(p2-p1) / (phantomRadius / (phantomDiameter_TP/2))
            
            #Add title with label and moment arm
            ax.set_title(f'Moment Arm Reliability Data ({muscle} / {round(ma,3)}mm)',
                         fontsize = 12, fontweight = 'bold')
            
            #Set tight layout
            plt.tight_layout()
            
            #Save figure
            plt.savefig(f'{currDir}reliabilityMomentArm.png', format = 'png', dpi = 300)
            
            #Close figure
            plt.close()
            
            #Append data
            #Original data
            momentArm1['specimen'].append(specimen)
            momentArm1['condition'].append(condition)
            momentArm1['load'].append(load)
            momentArm1['position'].append(position)
            momentArm1['muscle'].append(muscle)
            momentArm1['plane'].append(currPlane)
            momentArm1['val'].append(maOrig.loc[maOrig['subscapLine'] == muscle,]['momentArm'].values[0])
            #Reliability data
            momentArm2['specimen'].append(specimen)
            momentArm2['condition'].append(condition)
            momentArm2['load'].append(load)
            momentArm2['position'].append(position)
            momentArm2['muscle'].append(muscle)
            momentArm2['plane'].append(currPlane)
            momentArm2['val'].append(ma)
            
            # %% Line of aciton calculation
            
            #Create axes to visualise on
            fig, ax = plt.subplots(figsize = (12,6), nrows = 1, ncols = 2)

            #Show image
            ax[0].imshow(img, cmap = 'gray', origin = 'upper')
            
            #Turn off axes labels
            ax[0].axis('off')
            
            #Display subscapularis line
            #Plot points
            ax[0].scatter(ssPoints['X'], ssPoints['Y'], s = 5, c = 'red')
            #Fit line
            m, c = np.polyfit(ssPoints['X'], ssPoints['Y'], 1)
            #Plot fitted line
            ax[0].plot(ssPoints['X'], m * ssPoints['X'] + c, c = 'red', lw = 1)
            #Extend line of action to min and max x-values
            #Get current axes limits
            retXLim = ax[0].get_xlim()
            retYLim = ax[0].get_ylim()
            #Plot lines
            ax[0].plot(np.array((retXLim[0],np.min(ssPoints['X']))),
                      m * np.array((retXLim[0],np.min(ssPoints['X']))) + c,
                      c = 'red', lw = 1, ls = '--')
            ax[0].plot(np.array((retXLim[1],np.max(ssPoints['X']))),
                      m * np.array((retXLim[1],np.max(ssPoints['X']))) + c,
                      c = 'red', lw = 1, ls = '--')
            #Reset axes (note that axes don't seem to be changed, but here just in case)
            ax[0].set_xlim(retXLim)
            ax[0].set_ylim(retYLim)
            
            #Display glenoid plane
            ax[0].scatter(gpPoints['X'], gpPoints['Y'], s = 5, c = 'yellow')
            glenoidM,glenoidC = np.polyfit(gpPoints['X'], gpPoints['Y'], 1)
            ax[0].plot(gpPoints['X'], glenoidM * gpPoints['X'] + glenoidC, c = 'yellow', lw = 1)
            #Extend plane to max and min of x-axis
            ax[0].plot(np.array((retXLim[0],np.min(gpPoints['X']))),
                      glenoidM * np.array((retXLim[0],np.min(gpPoints['X']))) + glenoidC,
                      c = 'yellow', lw = 1, ls = '--')
            ax[0].plot(np.array((retXLim[1],np.max(gpPoints['X']))),
                      glenoidM * np.array((retXLim[1],np.max(gpPoints['X']))) + glenoidC,
                      c = 'yellow', lw = 1, ls = '--')
            #Reset axes (note that axes don't seem to be changed, but here just in case)
            ax[0].set_xlim(retXLim)
            ax[0].set_ylim(retYLim)
            
            #Calculate intersect of glenoid plane and lowest subscapularis
            #Set glenoid points as the top and bottom points
            #Extract two points on a line fit to the glenoid plane
            gx1 = 0
            gx2 = 100
            gy1 = glenoidM * gx1 + glenoidC #y-intercept; x = 0
            gy2 = glenoidM * gx2 + glenoidC #x = 100
            
            #Extract two points on the subscapularis linefit
            sx1 = 0
            sx2 = 100
            sy1 = m * sx1 + c #y-intercept; x = 0
            sy2 = m * sx2 + c #x = 100
            
            #Calculate line intersection
            intX, intY = lineIntersect(gx1, gy1, gx2, gy2, sx1, sy1, sx2, sy2)
            ax[0].scatter(intX, intY, c = 'blue', s = 5, zorder = 4)
            
            #Put a check in place to add a point to the subscap line if the points
            #are all on the wrong side of the intercept
            if currPlane == 'SP':
                
                #Check if all points are to the left of the intercept
                if np.all((intX - ssPoints['X'].to_numpy()) > 0):
                    #Calculate a point on the subscap line to the right of the intercept
                    #The makeshift X point is 5% across from the intercept relative
                    #to the maximum of the X axis                
                    makeshiftX = intX + (ax[0].get_xlim()[1] - intX) * 0.05
                    makeshiftY = m * makeshiftX + c
                    #Plot the point
                    ax[0].scatter(makeshiftX, makeshiftY, c = 'red', s = 5, zorder = 4)
                    #Add this point to the current subscapularis line for calculations
                    ssPoints = ssPoints.append({'X': makeshiftX, 'Y': makeshiftY}, ignore_index = True)  
                    
            elif currPlane == 'TP':
                
                if specimen.split('_')[1] == 'l':
                    
                    #Check if all points are to the left of the intercept
                    if np.all((intX - ssPoints['X'].to_numpy()) > 0):
                        #Calculate a point on the subscap line to the right of the intercept
                        #The makeshift X point is 5% across from the intercept relative
                        #to the maximum of the X axis                
                        makeshiftX = intX + (ax[0].get_xlim()[1] - intX) * 0.05
                        makeshiftY = m * makeshiftX + c
                        #Plot the point
                        ax[0].scatter(makeshiftX, makeshiftY, c = 'red', s = 5, zorder = 4)
                        #Add this point to the current subscapularis line for calculations
                        ssPoints = ssPoints.append({'X': makeshiftX, 'Y': makeshiftY}, ignore_index = True)
                        
                elif specimen.split('_')[1] == 'r':
                    
                    #Check if all points are to the left of the intercept
                    if np.all((intX - ssPoints['X'].to_numpy()) < 0):
                        #Calculate a point on the subscap line to the right of the intercept
                        #The makeshift X point is 5% across from the intercept relative
                        #to the maximum of the X axis                
                        makeshiftX = intX - (intX - ax[0].get_xlim()[0]) * 0.05
                        makeshiftY = m * makeshiftX + c
                        #Plot the point
                        ax[0].scatter(makeshiftX, makeshiftY, c = 'red', s = 5, zorder = 4)
                        #Add this point to the current subscapularis line for calculations
                        ssPoints = ssPoints.append({'X': makeshiftX, 'Y': makeshiftY}, ignore_index = True)
            
            #Determine slope (i.e. negative reciprocal) of perpendicular line to glenoid plane
            mPerp = 1 / (-glenoidM)
            
            #Solve for y-intercept if scapular plane
            if currPlane == 'SP':            
                #Solve for y-intercept of new line
                #Use intercept points for line
                cPerp = (mPerp * intX - intY) * -1
                #Visualise
                ax[0].plot([0,intX], [cPerp, intY], c = 'blue', lw = 1, ls = '--')
            #Solve for y at maximum or minimum X if transverse plane
            #Depends on if left or right limb for specimen
            elif currPlane == 'TP':
                if specimen.split('_')[1] == 'r':
                    #Solver for y at maximum x-value
                    cPerp = (mPerp * intX - intY) * -1
                    cMax = mPerp * img.shape[1] + cPerp
                    #Visualise
                    ax[0].plot([intX, img.shape[1]], [intY, cMax], c = 'blue', lw = 1, ls = '--')
                elif specimen.split('_')[1] == 'l':
                    #Solve for y-intercept of new line
                    #Use intercept points for line
                    cPerp = (mPerp * intX - intY) * -1
                    #Visualise
                    ax[0].plot([0,intX], [cPerp, intY], c = 'blue', lw = 1, ls = '--')
            
            #To calculate the line of action we revert to a standard lower left origin
            #compared to the upper left origin of the images --- therefore we can work
            #more easily by converting the data points to this coordinate system
            #Note that only the y coordinates need to change
            #Glenoid points
            gpPointsFlipX = gpPoints['X']
            gpPointsFlipY = img.shape[0] - gpPoints['Y']
            #Intersect point
            intFlipX = intX
            intFlipY = img.shape[0] - intY
            #Subscapularis points
            ssPointsFlipX = ssPoints['X']
            ssPointsFlipY = img.shape[0] - ssPoints['Y']
            
            #Reorient the data points so that the glenoid plane = y-axis and the perpendicular
            #line = x-axis
            
            #Recalculate glenoid plane gradient with flipped coordinates
            glenoidFlipM = np.polyfit(gpPointsFlipX, gpPointsFlipY, 1)[0]
                
            #Determine angle to rotate clockwise about intersection point based on
            #glenoid plane gradient
            #Negative gradient needs to be rotated clockwise by 90 - degrees
            if glenoidFlipM < 0:
                #Convert gradient to radians
                rotAng = np.radians(90) - np.arctan(glenoidM)
            elif glenoidFlipM > 0:
                #Convert gradient to radians
                rotAng = (np.radians(90) + np.arctan(glenoidM)) * -1
            elif glenoidFlipM == 0:
                #90 degree rotation required
                rotAng = np.radians(90)
            
            #Rotate the current subscapularis points clockwise by the prescribed degrees
            #Origin is the intersection point
            rx = []
            ry = []
            for nn in range(len(ssPointsFlipX)):
                #Set current point
                #Convert to 0,0 origin here too
                px = ssPointsFlipX[nn] - intFlipX
                py = ssPointsFlipY[nn] - intFlipY
                #Rotate about 0,0 origin by specified angle
                #Append to list here too
                rx.append(px * math.cos(rotAng) + py * math.sin(rotAng))
                ry.append(-px * math.sin(rotAng) + py * math.cos(rotAng))
                
            #Plot subscapularis points and origin
            ax[1].scatter(rx, ry, s = 15, c = 'red')
            ax[1].scatter(0, 0, s = 30, c = 'black')
            
            #Square up axes based on maximum value across the current axes limits
            #Get axes limits
            lnXLim = ax[1].get_xlim()
            lnYLim = ax[1].get_ylim()
            #Find absolute maximum of all values
            maxLim = np.max(np.array([np.abs(lnXLim),np.abs(lnYLim)]))
            #Set axes to square using maximum value
            ax[1].set_xlim([maxLim*-1,maxLim])
            ax[1].set_ylim([maxLim*-1,maxLim])
        
            #Create axes lines scaled to axes lengths
            if currPlane == 'SP':
                ax[1].plot([0,maxLim/2*-1], [0,0], lw = 2, c = 'black') #x-axis
            elif currPlane == 'TP':
                #Direction depends on specimen limb
                if specimen.split('_')[1] == 'r':
                    ax[1].plot([0,maxLim/2], [0,0], lw = 2, c = 'black') #x-axis
                elif specimen.split('_')[1] == 'l':
                    ax[1].plot([0,maxLim/2*-1], [0,0], lw = 2, c = 'black') #x-axis
            ax[1].plot([0,0], [0,maxLim/2], lw = 2, c = 'black') #y-axis
                
            #Fit muscle line
            rotM, rotC = np.polyfit(rx, ry, 1)
            ax[1].plot(np.array(rx), rotM * np.array(rx) + rotC, c = 'red', lw = 1, zorder = 0)
            #Extend to origin
            ax[1].plot(np.array((0,np.min(rx))),
                       rotM * np.array((0,np.min(rx))) + rotC,
                       c = 'red', lw = 1, ls = '--', zorder = 0)
                    
            #Turn off axes labels
            ax[1].axis('off')
        
            #Convert gradient of line to angle in degrees
            #Scapular plane calculations
            if currPlane == 'SP':
                if rotM > 0:
                    lineOfAction = 180 + np.degrees(np.arctan(rotM))
                elif rotM < 0:
                    #Inferiorly directed line of action (i.e. < 180 degrees)
                    lineOfAction = 180 - (np.degrees(np.arctan(rotM))*-1)
                elif rotM == 0:
                    #180 degree (i.e. straight compression) line of action
                    lineOfAction = 180
            #Transverse plane calculations
            elif currPlane == 'TP':
                #Calculations depend on specimen limb
                if specimen.split('_')[1] == 'r':
                    if rotM > 0:
                        lineOfAction = 180 + np.degrees(np.arctan(rotM))
                    elif rotM < 0:
                        #Inferiorly directed line of action (i.e. < 180 degrees)
                        lineOfAction = 180 - (np.degrees(np.arctan(rotM))*-1)
                    elif rotM == 0:
                        #180 degree (i.e. straight compression) line of action
                        lineOfAction = 180
                        
                elif specimen.split('_')[1] == 'l':
                    if rotM > 0:
                        lineOfAction = 180 - np.degrees(np.arctan(rotM))
                    elif rotM < 0:
                        #Inferiorly directed line of action (i.e. < 180 degrees)
                        lineOfAction = 180 + (np.degrees(np.arctan(rotM))*-1)
                    elif rotM == 0:
                        #180 degree (i.e. straight compression) line of action
                        lineOfAction = 180
            
            #Calculate stability ratio
            
            #Set a point at the end of the line of action
            #If this is in the scapular plane, we'll use the max x-limit
            if currPlane == 'SP':
                
                #Get the point at the edge
                srCalcPt = np.array([maxLim, rotM * maxLim + rotC])
                
            #If it is in the transverse plane, we'll use the max or min x-limit
            #This depends on the limb used
            if currPlane == 'TP':
                if '_r' in specimen:
                    #Get the point at the edge
                    srCalcPt = np.array([maxLim*-1, rotM * (maxLim*-1) + rotC])
                elif '_l' in specimen:
                    #Get the point at the edge
                    srCalcPt = np.array([maxLim, rotM * maxLim + rotC])
    
            #Calculate the directional cosines for the axes
            #This notes that the vector starts at the axes origin
            
            #Calculate the length of the vector between the origin and end point
            vecLen = np.sqrt(srCalcPt[0]**2 + srCalcPt[1]**2)
            
            #Calculate the directional cosines for the X and Y axes
            cosX = srCalcPt[0] / vecLen
            cosY = srCalcPt[1] / vecLen
            
            #Calculate stability ratio
            stabRatio = cosY / np.abs(cosX)
            
            #Print stability ratio in bootom right corner of axes
            ax[1].text(1-0.025, 0.1,
                      'SR = '+str(np.round(stabRatio,3)),
                      ha = 'right', va = 'center', transform = ax[1].transAxes)

            
        #Add title with label and line of action/stability ratio
        ax[0].set_title('Line of Action Reliability Data ('+muscle+' / '+str(np.round(lineOfAction,2))+u'\u00b0'+')',
                        fontsize = 12, fontweight = 'bold')
        ax[1].set_title('Stability Ratio Reliability Data ('+muscle+' / '+str(np.round(stabRatio,3))+')',
                        fontsize = 12, fontweight = 'bold')
        
        #Tight layout
        plt.tight_layout()
        
        #Save figure
        plt.savefig(f'{currDir}reliabilityLineOfAction.png', format = 'png', dpi = 300)
        
        #Close figure
        plt.close()
        
        #Append data
        
        #Line of action
        #Original data
        lineOfAction1['specimen'].append(specimen)
        lineOfAction1['condition'].append(condition)
        lineOfAction1['load'].append(load)
        lineOfAction1['position'].append(position)
        lineOfAction1['muscle'].append(muscle)
        lineOfAction1['plane'].append(currPlane)
        lineOfAction1['val'].append(loaOrig.loc[loaOrig['subscapLine'] == muscle,]['lineOfAction'].values[0])
        #Reliability data
        lineOfAction2['specimen'].append(specimen)
        lineOfAction2['condition'].append(condition)
        lineOfAction2['load'].append(load)
        lineOfAction2['position'].append(position)
        lineOfAction2['muscle'].append(muscle)
        lineOfAction2['plane'].append(currPlane)
        lineOfAction2['val'].append(lineOfAction)
        
        #Stability ratio
        stabilityRatio1['specimen'].append(specimen)
        stabilityRatio1['condition'].append(condition)
        stabilityRatio1['load'].append(load)
        stabilityRatio1['position'].append(position)
        stabilityRatio1['muscle'].append(muscle)
        stabilityRatio1['plane'].append(currPlane)
        stabilityRatio1['val'].append(srOrig.loc[srOrig['subscapLine'] == muscle,]['stabilityRatio'].values[0])
        #Reliability data
        stabilityRatio2['specimen'].append(specimen)
        stabilityRatio2['condition'].append(condition)
        stabilityRatio2['load'].append(load)
        stabilityRatio2['position'].append(position)
        stabilityRatio2['muscle'].append(muscle)
        stabilityRatio2['plane'].append(currPlane)
        stabilityRatio2['val'].append(stabRatio)

# %% Convert reliability results to dataframes

#Merge dataframes of different metrics
#Moment arms
maAgreement = pd.merge(pd.DataFrame.from_dict(momentArm1),
                       pd.DataFrame.from_dict(momentArm2),
                       how = 'left',
                       left_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'],
                       right_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'])
#Line of action
loaAgreement = pd.merge(pd.DataFrame.from_dict(lineOfAction1),
                        pd.DataFrame.from_dict(lineOfAction2),
                        how = 'left',
                        left_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'],
                        right_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'])
#Stability ratio
srAgreement = pd.merge(pd.DataFrame.from_dict(stabilityRatio1),
                       pd.DataFrame.from_dict(stabilityRatio2),
                       how = 'left',
                       left_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'],
                       right_on = ['specimen', 'condition', 'load', 'muscle', 'position', 'plane'])

#Calculate difference
maAgreement['difference'] = maAgreement['val_x'] - maAgreement['val_y']
loaAgreement['difference'] = loaAgreement['val_x'] - loaAgreement['val_y']
srAgreement['difference'] = srAgreement['val_x'] - srAgreement['val_y']

#Split into the intra and inter datasets
#Moment arms
maAgreementInter = maAgreement.loc[maAgreement['condition'] == 'split50',]
maAgreementIntra = maAgreement.loc[maAgreement['condition'] != 'split50',]
#Line of action
loaAgreementInter = loaAgreement.loc[loaAgreement['condition'] == 'split50',]
loaAgreementIntra = loaAgreement.loc[loaAgreement['condition'] != 'split50',]
#Stability ratio
srAgreementInter = srAgreement.loc[srAgreement['condition'] == 'split50',]
srAgreementIntra = srAgreement.loc[srAgreement['condition'] != 'split50',]

#Plot and save 95% LoA's

#Moment arms

#Inter-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(maAgreementInter['val_x'], maAgreementInter['val_y'])
ax.set_ylabel('Rater 1 - Rater 2')
ax.set_xlabel('Mean of Rater 1 & Rater 2')
ax.set_title('Inter-Rater Limits of Agreement: Moment Arms',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\InterRater\\Figures\\InterRater_MomentArm_LimitsOfAgreement.png', format = 'png', dpi = 300)
plt.close()

#Intra-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(maAgreementIntra['val_x'], maAgreementIntra['val_y'])
ax.set_ylabel('Rater 1A - Rater 1B')
ax.set_xlabel('Mean of Rater 1A & Rater 1B')
ax.set_title('Intra-Rater Limits of Agreement: Moment Arms',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\IntraRater\\Figures\\IntraRater_MomentArm_LimitsOfAgreement.png', format = 'png', dpi = 300)

#Line of Action

#Inter-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(loaAgreementInter['val_x'], loaAgreementInter['val_y'])
ax.set_ylabel('Rater 1 - Rater 2')
ax.set_xlabel('Mean of Rater 1 & Rater 2')
ax.set_title('Inter-Rater Limits of Agreement: Line of Action',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\InterRater\\Figures\\InterRater_LineOfAction_LimitsOfAgreement.png', format = 'png', dpi = 300)
plt.close()

#Intra-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(loaAgreementIntra['val_x'], loaAgreementIntra['val_y'])
ax.set_ylabel('Rater 1A - Rater 1B')
ax.set_xlabel('Mean of Rater 1A & Rater 1B')
ax.set_title('Intra-Rater Limits of Agreement: Line of Action',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\IntraRater\\Figures\\IntraRater_LineOfAction_LimitsOfAgreement.png', format = 'png', dpi = 300)
plt.close()

#Stability Ratio

#Inter-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(srAgreementInter['val_x'], srAgreementInter['val_y'])
ax.set_ylabel('Rater 1 - Rater 2')
ax.set_xlabel('Mean of Rater 1 & Rater 2')
ax.set_title('Inter-Rater Limits of Agreement: Stability Ratio',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\InterRater\\Figures\\InterRater_StabilityRatio_LimitsOfAgreement.png', format = 'png', dpi = 300)
plt.close()

#Intra-rater
fig, ax = plt.subplots(figsize = (6,6), nrows = 1, ncols = 1)
pg.plot_blandaltman(srAgreementIntra['val_x'], srAgreementIntra['val_y'])
ax.set_ylabel('Rater 1A - Rater 1B')
ax.set_xlabel('Mean of Rater 1A & Rater 1B')
ax.set_title('Intra-Rater Limits of Agreement: Stability Ratio',
             fontsize = 14, fontweight = 'bold')
plt.tight_layout()
plt.savefig('..\\Results\\Reliability\\IntraRater\\Figures\\IntraRater_StabilityRatio_LimitsOfAgreement.png', format = 'png', dpi = 300)
plt.close()

# %% ----- End of reliabilityAnalysis.py ----- %% #