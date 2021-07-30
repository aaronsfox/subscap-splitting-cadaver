# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:46:38 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is used to analyse data from the listed specimen across the Latarjet
    condition listed. It calculates subscapularis moment arms and lines of action
    for the measurement conditions that contain digitised data.
    
    This script should be run from the 'Code' folder, and with the appropriate
    specimen and Latarjet condition labels will be able to identify the data that
    needs to be processed.
	
	TODO: QUALITY CONTROL --- humeral head size consistency, phantom size consistency
	
	TODO: QUALITY CONTROL --- image series for matching across glenoid planes

	TODO: lack of phantom on certain images means we should take average where available (i.e.
	print out what the average and error is to confirm similar image scales across series)
    
"""

# %% Import packages

from matplotlib import pyplot as plt
import os
from glob import glob
import numpy as np
import pandas as pd
import math
from skimage import io
import re

# %% Parameters to change

##### ----- CHECK AND CHANGE THESE PARAMETERS EACH RUN ----- #####

#Specimen to analyse
specimen = 'S3_r' #specimen number and label for limb tested

#Latarjet condition to analyse
latarjetCond = 'split50' #split50, split25_upper, split25_lower

#Set phantom diameter for scaling
phantomDiameter = 6.5

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

#Plotting digitised points function
def visDigitised():
    
    #Plot image
    fig, ax = plt.subplots()
    ax.imshow(img, cmap = 'gray', origin = 'upper')
    
    #Turn off axes
    ax.axis('off')
    
    #Plot humeral head points
    ax.scatter(hhPoints['X'], hhPoints['Y'], s = 10, c = 'green')
    
    #Plot humeral head fitted circle
    ax.add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                             edgecolor = 'green', facecolor = 'none'))
    ax.scatter(hhCentreX, hhCentreY, s = 500, c = 'green', marker = '+')
    
    #Plot subscap points and fitted lines
    #Set point to extend subscapularis lines to
    #Currently half-way between image edge and humeral head centre
    #Need to set whether this is on the right or left of the image
    #Can do this by checking humeral head centre relative to glenoid plane
    if hhCentreX < np.mean(gpPoints['X']):
        extendX = hhCentreX / 2
    else:
        extendX = hhCentreX + ((img.shape[1] - hhCentreX) / 2)
    #Loop through points
    for ss in range(len(subscapNames)):
        #Plot points
        plt.scatter(subscapPoints[subscapNames[ss]]['X'],
                    subscapPoints[subscapNames[ss]]['Y'],
                    s = 10, c = 'red')
        #Fit line
        m,c = np.polyfit(subscapPoints[subscapNames[ss]]['X'],
                         subscapPoints[subscapNames[ss]]['Y'],
                         1)
        #Plot fitted line
        ax.plot(subscapPoints[subscapNames[ss]]['X'],
                m * subscapPoints[subscapNames[ss]]['X'] + c,
                c = 'red', lw = 1)
        #Extend line of action
        #Need to check which direction to extend in, and this will dictate whether
        #to use min or max
        if hhCentreX < np.mean(gpPoints['X']):
            ax.plot(np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))),
                    m * np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))) + c,
                    c = 'red', lw = 1, ls = '--')
        else:
            ax.plot(np.array((np.max(subscapPoints[subscapNames[ss]]['X']),extendX)),
                    m * np.array((np.max(subscapPoints[subscapNames[ss]]['X']),extendX)) + c,
                    c = 'red', lw = 1, ls = '--')
        
    #Plot phantom points
    if os.path.isfile('phantom.csv'):
        
        #Load file
        phantomPoints = pd.read_csv('phantom.csv')
        
        #Plot points
        ax.scatter(phantomPoints['X'], phantomPoints['Y'], s = 10, c = 'magenta')
        
        #Fit circle
        phantomCentreX, phantomCentreY, phantomRadius = fitCircle(phantomPoints['X'].to_numpy(),
                                                                  phantomPoints['Y'].to_numpy())
        
        #Plot phantom fitted circle
        ax.add_artist(plt.Circle((phantomCentreX, phantomCentreY), phantomRadius,
                                  edgecolor = 'magenta', facecolor = 'none'))
        ax.scatter(phantomCentreX, phantomCentreY, s = 75, c = 'magenta', marker = '+')
    
    #Plot glenoid plane points and fit line
    ax.scatter(gpPoints['X'], gpPoints['Y'], s = 10, c = 'yellow')
    m,c = np.polyfit(gpPoints['X'], gpPoints['Y'], 1)
    ax.plot(gpPoints['X'], m * gpPoints['X'] + c, c = 'yellow', lw = 1)
    
    #Set title
    plt.title(f'{currPosLabel} / {currLoadLabel} / {currPlaneLabel}\nDigitised Points & Fits',
              fontsize = 12, fontweight = 'bold')
    
    #Tight layout
    plt.tight_layout()

#Visualise glenoid plane fit
def visGlenoidFit():
    
    #Create figure and axes
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (8,4))
    
    #Plot reference image and digitised points              
    #Plot image
    ax[0].imshow(refImg, cmap = 'gray', origin = 'upper')
    #Plot originally digitised reference points
    ax[0].scatter(refPts['X'], refPts['Y'], c = 'green', s = 7)   
    #Plot glenoid points and fit a line
    ax[0].scatter(gpPts[:,0], gpPts[:,1], c = 'yellow', s = 7)
    m,c = np.polyfit(gpPts[:,0], gpPts[:,1], 1)
    ax[0].plot(gpPts[:,0], m * gpPts[:,0] + c, c = 'yellow', lw = 1)
    
    #Turn off axes
    ax[0].axis('off')
    
    #Plot transformed estimates on to original image
    ax[0].scatter(newPts[:,0], newPts[:,1],
                  edgecolor = 'red', facecolor = None,
                  s = 7)
    
    #Plot current glenoid image
    ax[1].imshow(img, cmap = 'gray', origin = 'upper')
    
    #Turn off axes
    ax[1].axis('off')
    
    #Plot current image digitised points
    ax[1].scatter(refPts['X'], refPts['Y'], c = 'red', s = 7)
    
    #Plot estimated glenoid plane points and fit line
    ax[1].scatter(gpPts_trans[:,0], gpPts_trans[:,1],
                  edgecolor = 'yellow', facecolor = None,
                  s = 7)
    m,c = np.polyfit(gpPts_trans[:,0], gpPts_trans[:,1], 1)
    ax[1].plot(gpPts_trans[:,0], m * gpPts_trans[:,0] + c, c = 'yellow',
               lw = 1, ls = '--')
    
    #Set titles
    ax[0].set_title(f'Ref. Img. / Res. Err. on Points: {np.round(resErr,5)}',
                    fontsize = 10, fontweight = 'bold')
    ax[1].set_title('Current Img. Estimated Glenoid Plane',
                    fontsize = 10, fontweight = 'bold')

    #Tight layout
    plt.tight_layout()

#Calculate and export moment arms
def calcMomentArms():
    
    #Set variable to store moment arm calculations in
    maList = []
    
    #Create axes to visualise on
    fig, ax = plt.subplots(figsize = (6,6), nrows = 2, ncols = 2)
    
    #Set variable for axes
    whichAx = [[0,0], [0,1],
               [1,0], [1,1]]
    
    #Set point to extend subscapularis line to based on humeral head
    #position (i.e. halfway between humeral head and edge of image)
    if hhCentreX < np.mean(gpPoints['X']):
        extendX = hhCentreX / 2
    else:
        extendX = hhCentreX + ((img.shape[1] - hhCentreX) / 2)
        
    #Set humeral head point as p3 for moment arm calculations
    p3 = np.asarray((hhCentreX, hhCentreY))
    
    #Loop through four for subscap lines
    for subscap in subscapNames:
        
        #Set current index
        ss = subscapNames.index(subscap)
        
        #Display image on current axes
        ax[whichAx[ss][0],whichAx[ss][1]].imshow(img,
                                                 cmap = 'gray',
                                                 origin = 'upper')
        
        #Turn off axes
        ax[whichAx[ss][0],whichAx[ss][1]].axis('off')
        
        #Display humeral head points
        ax[whichAx[ss][0],whichAx[ss][1]].scatter(hhPoints['X'], hhPoints['Y'], s = 10, c = 'green')
        ax[whichAx[ss][0],whichAx[ss][1]].add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                                 edgecolor = 'green', facecolor = 'none'))
        ax[whichAx[ss][0],whichAx[ss][1]].scatter(hhCentreX, hhCentreY, s = 500, c = 'green', marker = '+')
        
        #Display subscapularis line
        #Plot points
        ax[whichAx[ss][0],whichAx[ss][1]].scatter(subscapPoints[subscap]['X'],
                                                  subscapPoints[subscap]['Y'],
                                                  s = 10, c = 'red')
        #Fit line
        m,c = np.polyfit(subscapPoints[subscap]['X'],
                         subscapPoints[subscap]['Y'],
                         1)
        #Plot fitted line
        ax[whichAx[ss][0],whichAx[ss][1]].plot(subscapPoints[subscap]['X'],
                                               m * subscapPoints[subscap]['X'] + c,
                                               c = 'red', lw = 1)
        #Extend line of action
        #Check which way to extend and whether to use min or max
        if extendX < hhCentreX:
            ax[whichAx[ss][0],whichAx[ss][1]].plot(np.array((extendX,np.min(subscapPoints[subscap]['X']))),
                                                   m * np.array((extendX,np.min(subscapPoints[subscap]['X']))) + c,
                                                   c = 'red', lw = 1, ls = '--')
        else:
            ax[whichAx[ss][0],whichAx[ss][1]].plot(np.array((np.max(subscapPoints[subscap]['X']),extendX)),
                                                   m * np.array((np.max(subscapPoints[subscap]['X']),extendX)) + c,
                                                   c = 'red', lw = 1, ls = '--')
            
        #Extract two points on the line for later moment arm calculation
        #Take the end points of the line on the x-axes
        p1 = np.asarray((ax[whichAx[ss][0],whichAx[ss][1]].get_xlim()[0],
                         m * ax[whichAx[ss][0],whichAx[ss][1]].get_xlim()[0] + c))
        p2 = np.asarray((ax[whichAx[ss][0],whichAx[ss][1]].get_xlim()[1],
                         m * ax[whichAx[ss][0],whichAx[ss][1]].get_xlim()[1] + c))
        
        #Calculate the point on the extended subscapularis line that is intercepted
        #by a perpendicular line from the humerus
        
        #Set the points at the start and the end of the subscapularis line
        x1 = hhCentreX/2
        y1 = m * (hhCentreX/2)+ c
        x2 = max(subscapPoints[subscap]['X'])
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
        ax[whichAx[ss][0],whichAx[ss][1]].plot(
            np.array((x3,x4)),np.array((y3,y4)),
            c = 'yellow', lw = 1, ls = '--')
        
        #Calculate moment arm distance
        #Includes scale by bead size
        
        #Old method using distance - prone to small errors
        # ma = (math.sqrt(((x3 - x4)**2) + ((y3 - y4)**2))) / phantomScale
        
        #Use cross product to calculate moment arm
        if 'SP' in currPlaneName:
            ma = np.cross(p2-p1,p3-p1) / np.linalg.norm(p2-p1) / phantomScale_SP
        elif 'TP' in currPlaneName:
            ma = np.cross(p2-p1,p3-p1) / np.linalg.norm(p2-p1) / phantomScale_TP
        
        # #Check if moment arm is positive (abduction/IR) vs. negative (adduction/ER)
        # #This considers proper orientation of image
        ### Not needed with cross product approach
        # if y3 < y4:
        #     #Humeral head centre is above moment arm intersection
        #     ma = ma*-1
            
        #Store moment arm value
        maList.append(ma)
        
        #Add title with label and moment arm
        ax[whichAx[ss][0],whichAx[ss][1]].set_title(f'{muscleTitles[ss]} / {round(ma,3)}mm',
                                                    fontsize = 12,
                                                    fontweight = 'bold')
        
    #Set tight layout on figure
    plt.tight_layout()
    
    #Save figure
    plt.savefig('momentArms.png', format = 'png', dpi = 300)
    
    #Close figure
    plt.close()
    
    #Convert moment arm data to dataframe and export
    pd.DataFrame(list(zip(subscapNames, maList)),
                 columns = ['subscapLine', 'momentArm']).to_csv(
                     'momentArmData.csv',
                     index = False)
                     
#Calculate lines of action and stability ratios
def calcLoA():

    #Set variable to store line of action and stability ratio calculations in
    loaList = []
    srList = []
    
    #Create figure to visualise LoA's on
    fig = plt.figure(figsize = (4,8))
    
    #Loop through subscapularis lines
    for subscap in subscapNames:
        
        #Set current index
        ss = subscapNames.index(subscap)
        
        #Add axes to show image on
        imAx = plt.subplot2grid((8,8), (ss*2,0), rowspan = 2, colspan = 4)
        
        #Show image
        imAx.imshow(img, cmap = 'gray', origin = 'upper')
        
        #Turn off axes labels
        imAx.axis('off')
        
        #Display subscapularis line
        #Plot points
        imAx.scatter(subscapPoints[subscap]['X'],
                     subscapPoints[subscap]['Y'],
                     s = 5, c = 'red')
        #Fit line
        subscapM,subscapC = np.polyfit(subscapPoints[subscap]['X'],
                                       subscapPoints[subscap]['Y'],
                                       1)
        #Plot fitted line
        imAx.plot(subscapPoints[subscap]['X'],
                  subscapM * subscapPoints[subscap]['X'] + subscapC,
                  c = 'red', lw = 1)
        #Extend line of action to min and max x-values
        #Get current axes limits
        retXLim = imAx.get_xlim()
        retYLim = imAx.get_ylim()
        #Plot lines
        imAx.plot(np.array((retXLim[0],np.min(subscapPoints[subscap]['X']))),
                  subscapM * np.array((retXLim[0],np.min(subscapPoints[subscap]['X']))) + subscapC,
                  c = 'red', lw = 1, ls = '--')
        imAx.plot(np.array((retXLim[1],np.max(subscapPoints[subscap]['X']))),
                  subscapM * np.array((retXLim[1],np.max(subscapPoints[subscap]['X']))) + subscapC,
                  c = 'red', lw = 1, ls = '--')
        #Reset axes (note that axes don't seem to be changed, but here just in case)
        imAx.set_xlim(retXLim)
        imAx.set_ylim(retYLim)
        
        #Add title
        imAx.set_title(muscleTitles[ss],
                       fontsize = 12, fontweight = 'bold')
        
        #Display glenoid plane
        imAx.scatter(gpPoints['X'], gpPoints['Y'], s = 5, c = 'yellow')
        glenoidM,glenoidC = np.polyfit(gpPoints['X'], gpPoints['Y'], 1)
        imAx.plot(gpPoints['X'], glenoidM * gpPoints['X'] + glenoidC, c = 'yellow', lw = 1)
        #Extend plane to max and min of x-axis
        imAx.plot(np.array((retXLim[0],np.min(gpPoints['X']))),
                  glenoidM * np.array((retXLim[0],np.min(gpPoints['X']))) + glenoidC,
                  c = 'yellow', lw = 1, ls = '--')
        imAx.plot(np.array((retXLim[1],np.max(gpPoints['X']))),
                  glenoidM * np.array((retXLim[1],np.max(gpPoints['X']))) + glenoidC,
                  c = 'yellow', lw = 1, ls = '--')
        
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
        sy1 = subscapM * sx1 + subscapC #y-intercept; x = 0
        sy2 = subscapM * sx2 + subscapC #x = 100
        
        #Calculate line intersection
        intX, intY = lineIntersect(gx1, gy1, gx2, gy2, sx1, sy1, sx2, sy2)
        imAx.scatter(intX, intY, c = 'blue', s = 5, zorder = 4)
        
        #Determine slope (i.e. negative reciprocal) of perpendicular line to glenoid plane
        mPerp = 1 / (-glenoidM)
        
        #Solve for y-intercept if scapular plane
        if 'SP' in currPlaneName:            
            #Solve for y-intercept of new line
            #Use intercept points for line
            cPerp = (mPerp * intX - intY) * -1
            #Visualise
            imAx.plot([0,intX], [cPerp, intY], c = 'blue', lw = 1, ls = '--')
        #Solve for y at maximum or minimum X if transverse plane
        #Depends on if left or right limb for specimen
        elif 'TP' in currPlaneName:
            if specimen.split('_')[1] == 'r':
                #Solver for y at maximum x-value
                cPerp = (mPerp * intX - intY) * -1
                cMax = mPerp * img.shape[1] + cPerp
                #Visualise
                imAx.plot([intX, img.shape[1]], [intY, cMax], c = 'blue', lw = 1, ls = '--')
            elif specimen.split('_')[1] == 'l':
                #Solve for y-intercept of new line
                #Use intercept points for line
                cPerp = (mPerp * intX - intY) * -1
                #Visualise
                imAx.plot([0,intX], [cPerp, intY], c = 'blue', lw = 1, ls = '--')
        
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
        ssPointsFlipX = subscapPoints[subscapNames[ss]]['X']
        ssPointsFlipY = img.shape[0] - subscapPoints[subscapNames[ss]]['Y']
        
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
                
        #Visualise rotated points on new subplots
        lnAx = plt.subplot2grid((8,8), (ss*2,4), rowspan = 2, colspan = 4)
    
        #Plot subscapularis points and origin
        lnAx.scatter(rx, ry, s = 15, c = 'red')
        lnAx.scatter(0, 0, s = 30, c = 'black')
        
        #Square up axes based on maximum value across the current axes limits
        #Get axes limits
        lnXLim = lnAx.get_xlim()
        lnYLim = lnAx.get_ylim()
        #Find absolute maximum of all values
        maxLim = np.max(np.array([np.abs(lnXLim),np.abs(lnYLim)]))
        #Set axes to square using maximum value
        lnAx.set_xlim([maxLim*-1,maxLim])
        lnAx.set_ylim([maxLim*-1,maxLim])
    
        #Create axes lines scaled to axes lengths
        if 'SP' in currPlaneName:
            lnAx.plot([0,maxLim/2*-1], [0,0], lw = 2, c = 'black') #x-axis
        elif 'TP' in currPlaneName:
            #Direction depends on specimen limb
            if specimen.split('_')[1] == 'r':
                lnAx.plot([0,maxLim/2], [0,0], lw = 2, c = 'black') #x-axis
            elif specimen.split('_')[1] == 'l':
                lnAx.plot([0,maxLim/2*-1], [0,0], lw = 2, c = 'black') #x-axis
        lnAx.plot([0,0], [0,maxLim/2], lw = 2, c = 'black') #y-axis
        
        #Fit muscle line
        rotM, rotC = np.polyfit(rx, ry, 1)
        lnAx.plot(np.array(rx), rotM * np.array(rx) + rotC, c = 'red', lw = 1, zorder = 0)
        #Extend to origin
        lnAx.plot(np.array((0,np.min(rx))),
                  rotM * np.array((0,np.min(rx))) + rotC,
                  c = 'red', lw = 1, ls = '--', zorder = 0)
        
        #### TODO: EXTEND TO APPROPRIATE EDGE OF PLOT FOR BEST VISUALISATION...
        
        #Turn off axes labels
        lnAx.axis('off')
    
        #Convert gradient of line to angle in degrees
        #Scapular plane calculations
        if 'SP' in currPlaneName:
            if rotM > 0:
                lineOfAction = 180 + np.degrees(np.arctan(rotM))
            elif rotM < 0:
                #Inferiorly directed line of action (i.e. < 180 degrees)
                lineOfAction = 180 - (np.degrees(np.arctan(rotM))*-1)
            elif rotM == 0:
                #180 degree (i.e. straight compression) line of action
                lineOfAction = 180
        #Transverse plane calculations
        elif 'TP' in currPlaneName:
            #Calculations depend on specimen limb
            if specimen.split('_')[1] == 'r':
                if rotM > 0:
                    lineOfAction = 180 + np.degrees(np.arctan(rotM))
                elif rotM < 0:
                    #Inferiorly directed line of action (i.e. < 180 degrees)
                    lineOfAction = 360 - (np.degrees(np.arctan(rotM))*-1)
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
                    
        #Print line of action in top left corner of axes
        lnAx.text(0.025, 0.9,
                  'LoA = '+str(np.round(lineOfAction,2))+u'\u00b0',
                  ha = 'left', va = 'center', transform = lnAx.transAxes)
        
        #Store current LoA calculation in list
        loaList.append(lineOfAction)
        
        #Calculate stability ratio
        
        #Set a point at the end of the line of action
        #If this is in the scapular plane, we'll use the max x-limit
        if 'SP' in currPlaneName:
            
            #Get the point at the edge
            srCalcPt = np.array([maxLim, rotM * maxLim + rotC])
            
        #If it is in the transverse plane, we'll use the min x-limit
        if 'TP' in currPlaneName:
            
            #Get the point at the edge
            srCalcPt = np.array([maxLim*-1, rotM * (maxLim*-1) + rotC])
            
        #Calculate the directional cosines for the axes
        #This notes that the vector starts at the axes origin
        
        #Calculate the length of the vector between the origin and end point
        vecLen = np.sqrt(srCalcPt[0]**2 + srCalcPt[1]**2)
        
        #Calculate the directional cosines for the X and Y axes
        cosX = srCalcPt[0] / vecLen
        cosY = srCalcPt[1] / vecLen
        
        #Calculate stability ratio
        stabRatio = cosY / np.abs(cosX)
            
        #Store current stability ratio calculation in list
        srList.append(stabRatio)
        
        #Print stability ratio in bootom right corner of axes
        lnAx.text(1-0.025, 0.1,
                  'SR = '+str(np.round(stabRatio,3)),
                  ha = 'right', va = 'center', transform = lnAx.transAxes)
        
    #Tight layout
    plt.tight_layout()
    
    #Save figure
    plt.savefig('linesOfAction.png', format = 'png', dpi = 300)
    
    #Close figure
    plt.close()
    
    #Convert moment arm data to dataframe and export
    pd.DataFrame(list(zip(subscapNames, loaList)),
                 columns = ['subscapLine', 'lineOfAction']).to_csv(
                     'lineOfActionData.csv',
                     index = False)

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

#Create labels

#Position names and labels
posNames = ['abd0', 'abd30', 'abd60', 'abd90', 'ABER', 'APP']
posLabels = ['0'+u'\u00b0'+' Abd.', '30'+u'\u00b0'+' Abd.',
             '60'+u'\u00b0'+' Abd.', '90'+u'\u00b0'+' Abd.',
             'Abd. Ext. Rot.', 'Apprehension']

#Load names and labels
loadNames = ['_0N', '_10N', '_20N', '_30N', '_40N']
loadLabels = ['0N', '10N', '20N', '30N', '40N']

#Plane names and labels
planeNames = ['_SP', '_TP']
planeLabels = ['Scapular Plane', 'Transverse Plane']

#Subscapularis lines
subscapNames = ['ss1', 'ss2', 'ss3', 'ss4']

#Titles for muscle lines
muscleTitles = ['Superior', 'Middle-Superior', 'Middle-Inferior', 'Inferior']

#Set code directory to return to
homeDir = os.getcwd()

#Navigate to data directory
os.chdir(f'..\\Data\\{specimen}\\{latarjetCond}')

#Get list of directories to work through
dirList = glob(os.getcwd()+'\\*\\')

# %% Get reference data for estimating glenoid planes

#Get global directory list for each plane
dirList_SP = [ii for ii in dirList if '_SP_' in ii]
dirList_TP = [ii for ii in dirList if '_TP_' in ii]

#Loop through directory lists and identify which have defined glenoid planes
#If they do they are added as a key in the dictionary and set as the current
#definition directory. Any subsequent directories without a glenoid plane are
#mapped to that key to use it in defining it's glenoid plane

#There are cases though where the initial directories may be searched through
#and not have a glenoid plane definition directory selected already. Hence we
#first need to loop through the directories to find the first glenoid definition
#directory to set as the starting 'currDefDir'

#Scapular plane
glenoidDefCond_SP = {}

#Check for starting directory
for currDir in dirList_SP:
    #Navigate to directory
    os.chdir(currDir)
    #Check if glenoid plane file present
    if os.path.isfile('gp.csv'):
        #Add as key to dictionary
        glenoidDefCond_SP[currDir] = []
        #Set as the starting definition directory
        currDefDir = currDir
        #Exit loop with break
        break
#Return to home directory
os.chdir('..')

#Loop through directories
for currDir in dirList_SP:
    #Navigate to directory
    os.chdir(currDir)
    #Check if glenoid plane file present
    if os.path.isfile('gp.csv'):
        #Add as key to dictionary
        #Check first if it's already in there from initial definition
        if currDir not in glenoidDefCond_SP.keys():
            glenoidDefCond_SP[currDir] = []
        #Set as current definition directory
        currDefDir = currDir
    else:
        #Append the current directory to the current definition key
        glenoidDefCond_SP[currDefDir].append(currDir)
#Return to home directory
os.chdir('..')

#Transverse plane
glenoidDefCond_TP = {}

#Check for starting directory
for currDir in dirList_TP:
    #Navigate to directory
    os.chdir(currDir)
    #Check if glenoid plane file present
    if os.path.isfile('gp.csv'):
        #Add as key to dictionary
        glenoidDefCond_TP[currDir] = []
        #Set as the starting definition directory
        currDefDir = currDir
        #Exit loop with break
        break
#Return to home directory
os.chdir('..')

#Loop through directories
for currDir in dirList_TP:
    #Navigate to directory
    os.chdir(currDir)
    #Check if glenoid plane file present
    if os.path.isfile('gp.csv'):
        #Add as key to dictionary
        #Check first if it's already in there from initial definition
        if currDir not in glenoidDefCond_TP.keys():
            glenoidDefCond_TP[currDir] = []
        #Set as current definition directory
        currDefDir = currDir
    else:
        #Append the current directory to the current definition key
        glenoidDefCond_TP[currDefDir].append(currDir)
#Return to home directory
os.chdir('..')
        
#Navigate to the definition directories and read in glenoid data

#Scapular plane

#Create dictionary to store gp, pts and img data
gpData_SP = {val:{'gp': [], 'pts':[], 'refImg':[]} for val in list(glenoidDefCond_SP.keys())}

#Loop through the definition directories and grab data
for defDir in gpData_SP.keys():
    #Navigate to directory
    os.chdir(defDir)
    #Load glenoid plane data and append to dictionary
    gpData_SP[defDir]['gp'].append(pd.read_csv('gp.csv'))
    #Load points data and append to dictionary
    gpData_SP[defDir]['pts'].append(pd.read_csv('pts.csv'))
    #Load ref image data and append to dictionary
    gpData_SP[defDir]['refImg'].append(io.imread(glob('*_proc.tif')[0]))
    #Return to home directory
    os.chdir('..')

#Transverse plane
    
#Create dictionary to store gp, pts and img data
gpData_TP = {val:{'gp': [], 'pts':[], 'refImg':[]} for val in list(glenoidDefCond_TP.keys())}

#Loop through the definition directories and grab data
for defDir in gpData_TP.keys():
    #Navigate to directory
    os.chdir(defDir)
    #Load glenoid plane data and append to dictionary
    gpData_TP[defDir]['gp'].append(pd.read_csv('gp.csv'))
    #Load points data and append to dictionary
    gpData_TP[defDir]['pts'].append(pd.read_csv('pts.csv'))
    #Load ref image data and append to dictionary
    gpData_TP[defDir]['refImg'].append(io.imread(glob('*_proc.tif')[0]))
    #Return to home directory
    os.chdir('..')

# %% Create scaling factor based on phantom size
    
##### TODO: consider if a plane doesn't have any phantoms?????    

#Taking an average became necessary given the phantom was not viewable in every image

#Scapular plane
    
#Set list to store phantom sizes in
phantomScaleAll = []

#Loop through the directories
for currDir in dirList_SP:
    
    #Navigate to directory
    os.chdir(currDir)
    
    #Load the phantom if available
    if os.path.isfile('phantom.csv'):
        
        #Load the points
        phantomPoints = pd.read_csv('phantom.csv')
        
        #Fit circle to phantom bead
        _, _, phantomRadius = fitCircle(phantomPoints['X'].to_numpy(),
                                        phantomPoints['Y'].to_numpy())
        
        #Add to scaling size list            
        phantomScaleAll.append(phantomRadius / (phantomDiameter/2))
        
    #Return to home directory
    os.chdir('..')
    
#Take average scale size to apply
#Print mean and SD to check consistency
phantomScale_SP = np.mean(phantomScaleAll)
print(f'Scapular plane phantom scale: {np.round(phantomScale_SP,2)} \u00B1 {np.round(np.std(phantomScaleAll),2)}')

#Transverse plane

#Set list to store phantom sizes in
phantomScaleAll = []

#Loop through the directories
for currDir in dirList_TP:
    
    #Navigate to directory
    os.chdir(currDir)
    
    #Load the phantom if available
    if os.path.isfile('phantom.csv'):
        
        #Load the points
        phantomPoints = pd.read_csv('phantom.csv')
        
        #Fit circle to phantom bead
        _, _, phantomRadius = fitCircle(phantomPoints['X'].to_numpy(),
                                        phantomPoints['Y'].to_numpy())
        
        #Add to scaling size list            
        phantomScaleAll.append(phantomRadius / (phantomDiameter/2))
        
    #Return to home directory
    os.chdir('..')
    
#Take average scale size to apply
#Print mean and SD to check consistency
phantomScale_TP = np.mean(phantomScaleAll)
print(f'Scapular plane phantom scale: {np.round(phantomScale_TP,2)} \u00B1 {np.round(np.std(phantomScaleAll),2)}')    

# %% Perform calculations

#Loop through directories
for currDir in dirList:
    
    #Navigate to current directory
    os.chdir(currDir)
    
    #Search for the presence of digitised data
    #Will use the 'ss4.csv' file for this
    if not glob('ss4.csv'):
        #Set to not analyse this folder
        analyseDir = False
    else:
        analyseDir = True
		
	##### TODO: check for all necessary files - flag error or print notification???
        
    #Run calculations if data is available
    if analyseDir:
        
        # %% Import data
        
        #Identify the processed .tif file
        imageFile = glob('*_proc.tif')[0]
        
        #Load image
        img = io.imread(imageFile)
        
        #Identify analysis conditions
        
        #Position
        currPos = []
        for pos in posNames:
            if re.search(pos, imageFile):
                #Set the current position to this label
                currPos.append(pos)
        #Check for presence of only one position
        if len(currPos) != 1:
            #Raise an error
            raise ValueError('Unable to identify current position.')
        else:
            #Set position name
            currPosName = currPos[0]
            #Set current position label
            currPosInd = posNames.index(currPosName)
            currPosLabel = posLabels[currPosInd]
            
        #Load
        currLoad = []
        for load in loadNames:
            if re.search(load, imageFile):
                #Set the current position to this label
                currLoad.append(load)
        #Check for presence of only one position
        if len(currLoad) != 1:
            #Raise an error
            raise ValueError('Unable to identify current load.')
        else:
            #Set position name
            currLoadName = currLoad[0]
            #Set current position label
            currLoadInd = loadNames.index(currLoadName)
            currLoadLabel = loadLabels[currLoadInd]
            
        #Plane
        currPlane = []
        for plane in planeNames:
            if re.search(plane, imageFile):
                #Set the current position to this label
                currPlane.append(plane)
        #Check for presence of only one position
        if len(currPlane) != 1:
            #Raise an error
            raise ValueError('Unable to identify current plane.')
        else:
            #Set position name
            currPlaneName = currPlane[0]
            #Set current position label
            currPlaneInd = planeNames.index(currPlaneName)
            currPlaneLabel = planeLabels[currPlaneInd]
        
        #Import the digitised data
        
        #Subscapularis points
        #Create dictionary to store in
        subscapPoints = {}
        #Get data
        for ss in range(len(subscapNames)):
            #Import .csv as dataframe and store to dataframe
            subscapPoints[subscapNames[ss]] = pd.read_csv(subscapNames[ss]+'.csv')
            
        #Humeral head points
        hhPoints = pd.read_csv('hh.csv')
        
        #Glenoid plane
        #Check if this condition is the one where the plane was defined
        if currDir in list(glenoidDefCond_SP.keys()) or currDir in list(glenoidDefCond_TP.keys()):
            
            #Simply load in glenoid points
            gpPoints = pd.read_csv('gp.csv')
            
        else:
            
            #Estimate glenoid plane based on image transformation
            
            #Import the digitised reference points from the current image
            refPts = pd.read_csv('pts.csv')
            
            #Identify which reference image should be used for current directory
            if 'SP' in currPlaneName:
                #Search through scapular plane dictionary
                for checkDir in glenoidDefCond_SP.keys():
                    if currDir in glenoidDefCond_SP[checkDir]:
                        #Set the reference directory to current option
                        refDir = checkDir
                        #Break the loop to exit if done
                        break
            elif 'TP' in currPlaneName:
                #Search through scapular plane dictionary
                for checkDir in glenoidDefCond_SP.keys():
                    if currDir in glenoidDefCond_TP[checkDir]:
                        #Set the reference directory to current option
                        refDir = checkDir
                        #Break the loop to exit if done
                        break

            #Calculate procrustes transformation between points
                    
            #Extract the appropriate points to an array
            if 'SP' in currPlaneName:
                refDigPts = np.transpose(np.array([gpData_SP[refDir]['pts'][0]['X'].to_numpy(),
                                                   gpData_SP[refDir]['pts'][0]['Y'].to_numpy()]))
            elif 'TP' in currPlaneName:
                refDigPts = np.transpose(np.array([gpData_TP[refDir]['pts'][0]['X'].to_numpy(),
                                                   gpData_TP[refDir]['pts'][0]['Y'].to_numpy()]))

            #Set points to transform
            transDigPts =  np.transpose(np.array([refPts['X'].to_numpy(),
                                                  refPts['Y'].to_numpy()]))
                
            #Apply procrustes algorithm
            resErr, newPts, tform = procrustes(refDigPts, transDigPts)
            
            #Apply the transformation to the reference glenoid plane points
            #Get gp points in array based on current plane
            if 'SP' in currPlaneName:
                gpPts = np.transpose(np.array([gpData_SP[refDir]['gp'][0]['X'],
                                               gpData_SP[refDir]['gp'][0]['Y']]))
            elif 'TP' in currPlaneName:
                gpPts = np.transpose(np.array([gpData_TP[refDir]['gp'][0]['X'],
                                               gpData_TP[refDir]['gp'][0]['Y']]))
            
            #Transform reference glenoid points to new array
            gpPts_trans = np.zeros(gpPts.shape)
            for pp in range(gpPts.shape[0]):
                gpPts_trans[pp] = tform['rotation'].dot(gpPts[pp]) - tform['translation']
             
            #Extract the relevant reference image
            if 'SP' in currPlaneName:
                refImg = gpData_SP[refDir]['refImg'][0]
            elif 'TP' in currPlaneName:
                refImg = gpData_TP[refDir]['refImg'][0]            
                
            #Visualise new points on image alongside point alignment comparison
            visGlenoidFit()
            #Save figure
            plt.savefig('estimatedGlenoidPlaneFit.png', format = 'png', dpi = 300)
            #Close figure
            plt.close()
            
            #Convert fitted glenoid plane to same dataframe format as digitised
            gpPoints = pd.DataFrame(gpPts_trans, columns = ['X','Y'])
            
            #Print out residual error to folder for reference
            with open('glenoidPlaneFit_resErr.txt', 'w') as f:
                f.write('%f' % resErr)
            f.close()

        # %% Calculate moment arms
        
        #Fit circle to humeral head
        hhCentreX, hhCentreY, hhRadius = fitCircle(hhPoints['X'].to_numpy(),
                                                   hhPoints['Y'].to_numpy())
        
        #Visualise imported points
        visDigitised()
        #Save figure
        plt.savefig('digitisedPointsFitted.png', format = 'png', dpi = 300)
        #Close figure
        plt.close()
        
        #Moment arm calculations and visualisation
        calcMomentArms()
                             
        # %% Calculate lines of action & stability ratios
                             
        #Line of action calculations and visualisation
        calcLoA()
        
    # %% Finish up in current directory
    
    #Navigate back up to data directory
    os.chdir('..')

# %% -----

    