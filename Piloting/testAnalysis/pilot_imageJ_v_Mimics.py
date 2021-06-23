# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 09:43:17 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Code compares the outputs of the original Mimics approach vs. the ImageJ
    approach to determine feasibility and ease of use.
    
"""

# %% Import packages

from matplotlib import pyplot as plt
import csv
import numpy as np
import pandas as pd
import math
from skimage import io

# %% TODO

##### Clean up throughout --- number of repetitious code getting the same things

# %% Import adjusted ImageJ image

#Create image object from file
img = io.imread('imageJ_adjustedImage.tif')

#View image
plt.imshow(img, cmap = 'gray', origin = 'upper')

# %% Import digitised items from ImageJ

#Subscapularis points
#Create labels
subscapNames = ['ss1', 'ss2', 'ss3', 'ss4']
#Create dictionary to store in
subscapPoints = {}
#Get data
for ss in range(len(subscapNames)):
    #Import .csv as dataframe and store to dataframe
    subscapPoints[subscapNames[ss]] = pd.read_csv(subscapNames[ss]+'.csv')

#Humeral head points
hhPoints = pd.read_csv('hh.csv')

#Glenoid planes
gpPoints = pd.read_csv('gp.csv')

#Phantom points
phantomPoints = pd.read_csv('phantom.csv')

### TODO: need NSEW points to help with image orientation if it's rotated?

# %% Fit a circle to the humeral head points

##### TODO: clean up function...define as function...

##### TODO: consider number of points on humerus?

#Calculate means of XY points
x_m = np.mean(hhPoints['X'].to_numpy())
y_m = np.mean(hhPoints['Y'].to_numpy())

#Calculation of the reduced coordinates
u = hhPoints['X'].to_numpy() - x_m
v = hhPoints['Y'].to_numpy() - y_m

# linear system defining the center (uc, vc) in reduced coordinates:
#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
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
hhCentreX = x_m + uc
hhCentreY = y_m + vc

#Calculate distance to centre and radius
Ri = np.sqrt((hhPoints['X'].to_numpy()-hhCentreX)**2 + (hhPoints['Y'].to_numpy()-hhCentreY)**2)
hhRadius = np.mean(Ri)
# residuals = sum((Ri-R)**2)

# %% Fit circle to phantom bead

### TODO: same as above...map to function...

#Calculate means of XY points
x_m = np.mean(phantomPoints['X'].to_numpy())
y_m = np.mean(phantomPoints['Y'].to_numpy())

#Calculation of the reduced coordinates
u = phantomPoints['X'].to_numpy() - x_m
v = phantomPoints['Y'].to_numpy() - y_m

# linear system defining the center (uc, vc) in reduced coordinates:
#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
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
phantomCentreX = x_m + uc
phantomCentreY = y_m + vc

#Calculate distance to centre and radius
Ri = np.sqrt((phantomPoints['X'].to_numpy()-phantomCentreX)**2 + (phantomPoints['Y'].to_numpy()-phantomCentreY)**2)
phantomRadius = np.mean(Ri)
# residuals = sum((Ri-R)**2)

# %% Visualise imported ImageJ points

#Plot image
fig, ax = plt.subplots()
ax.imshow(img, cmap = 'gray', origin = 'upper')

#Plot humeral head points
ax.scatter(hhPoints['X'], hhPoints['Y'], s = 10, c = 'green')

#Plot humeral head fitted circle
ax.add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                         edgecolor = 'green', facecolor = 'none'))
ax.scatter(hhCentreX, hhCentreY, s = 500, c = 'green', marker = '+')

#Plot subscap points and fitted lines
#Set point to extend subscapularis lines to
#Currently half-way between image edge and humeral head centre
extendX = hhCentreX / 2
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
    ax.plot(np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))),
            m * np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))) + c,
            c = 'red', lw = 1, ls = '--')
    
#Plot phantom points
ax.scatter(phantomPoints['X'], phantomPoints['Y'], s = 10, c = 'blue')

#Plot humeral head fitted circle
ax.add_artist(plt.Circle((phantomCentreX, phantomCentreY), phantomRadius,
                         edgecolor = 'blue', facecolor = 'none'))
ax.scatter(phantomCentreX, phantomCentreY, s = 125, c = 'blue', marker = '+')

#Plot glenoid plane points and fit line
ax.scatter(gpPoints['X'], gpPoints['Y'], s = 10, c = 'yellow')
m,c = np.polyfit(gpPoints['X'], gpPoints['Y'], 1)
ax.plot(gpPoints['X'], m * gpPoints['X'] + c, c = 'yellow', lw = 1)

# %% Moment arm calculations - scapular plane

#Calculate the scaling factor to make the bead radius = 3
phantomScale = phantomRadius / 3

##### TODO: clean up function...

#Visualise moment arm line of the subscapularis lines
fig, ax = plt.subplots(figsize = (10,10), nrows = 2, ncols = 2)

#Set variable for axes
whichAx = [[0,0], [0,1],
           [1,0], [1,1]]

#Set titles for muscle lines
muscleTitles = ['Superior', 'Middle-Superior', 'Middle-Inferior', 'Inferior']

#Loop through four for subscap lines
for ss in range(len(subscapNames)):
    
    #Display image on current axes
    ax[whichAx[ss][0],whichAx[ss][1]].imshow(img, cmap = 'gray', origin = 'upper')
    
    #Display humeral head points
    ax[whichAx[ss][0],whichAx[ss][1]].scatter(hhPoints['X'], hhPoints['Y'], s = 10, c = 'green')
    ax[whichAx[ss][0],whichAx[ss][1]].add_artist(plt.Circle((hhCentreX, hhCentreY), hhRadius,
                             edgecolor = 'green', facecolor = 'none'))
    ax[whichAx[ss][0],whichAx[ss][1]].scatter(hhCentreX, hhCentreY, s = 500, c = 'green', marker = '+')
    
    #Display subscapularis line
    #Plot points
    ax[whichAx[ss][0],whichAx[ss][1]].scatter(subscapPoints[subscapNames[ss]]['X'],
                                              subscapPoints[subscapNames[ss]]['Y'],
                                              s = 10, c = 'red')
    #Fit line
    m,c = np.polyfit(subscapPoints[subscapNames[ss]]['X'],
                     subscapPoints[subscapNames[ss]]['Y'],
                     1)
    #Plot fitted line
    ax[whichAx[ss][0],whichAx[ss][1]].plot(subscapPoints[subscapNames[ss]]['X'],
                                           m * subscapPoints[subscapNames[ss]]['X'] + c,
                                           c = 'red', lw = 1)
    #Extend line of action
    ax[whichAx[ss][0],whichAx[ss][1]].plot(np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))),
                                           m * np.array((extendX,np.min(subscapPoints[subscapNames[ss]]['X']))) + c,
                                           c = 'red', lw = 1, ls = '--')
    
    #Add title
    ax[whichAx[ss][0],whichAx[ss][1]].set_title(muscleTitles[ss])
    
    #Calculate the point on the extended subscapularis line that is intercepted
    #by a perpendicular line from the humerus
    
    #Set the points at the start and the end of the subscapularis line
    x1 = extendX
    y1 = m * extendX + c
    x2 = max(subscapPoints[subscapNames[ss]]['X'])
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
    #This is effectively the moment arm
    ax[whichAx[ss][0],whichAx[ss][1]].plot(
        np.array((x3,x4)),np.array((y3,y4)),
        c = 'yellow', lw = 1, ls = '--')
    
    #Calculate moment arm distance
    #Includes scale by bead size
    ma = (math.sqrt(((x3 - x4)**2) + ((y3 - y4)**2))) / phantomScale
    #Check if moment arm is positive (above humeral head) or negative (below humeral head)
    #This perhaps considers that humerus is aligned in image appropriately.
    #### TODO: may need to include some extra points to calculate a rotation
    #### that places the humerus exactly vertically
    if y3 < y4:
        #Humeral head centre is above moment arm intersection
        ma = ma*-1
    
    #Print moment arm output for frontal plane
    print(f'{muscleTitles[ss]} subscapularis moment arm: {round(ma,2)}')
    
### Only a very slight difference between ImageJ and Mimics. Likely due to slight
### difference in humeral head fitting procedures, perhaps slightly different digitisation
### of points --- but not a big enough difference to cause issues...
    
# %% Calculate lines of action - scapular plane

#Line intersection code
def lineIntersect(Ax1, Ay1, Ax2, Ay2, Bx1, By1, Bx2, By2):
    """ returns a (x, y) tuple or None if there is no intersection """
    d = (By2 - By1) * (Ax2 - Ax1) - (Bx2 - Bx1) * (Ay2 - Ay1)
    uA = ((Bx2 - Bx1) * (Ay1 - By1) - (By2 - By1) * (Ax1 - Bx1)) / d
    uB = ((Ax2 - Ax1) * (Ay1 - By1) - (Ay2 - Ay1) * (Ax1 - Bx1)) / d
    x = Ax1 + uA * (Ax2 - Ax1)
    y = Ay1 + uA * (Ay2 - Ay1)
 
    return x, y

#### TODO: visualisation below will only really work for scapular plane situation
#### with how the lines are extended

#Create figure to visualise on
fig = plt.figure(figsize = (6,12))
# fig, ax = plt.subplots(figsize = (10,10), nrows = 2, ncols = 2)

#Loop through subscapularis lines
for ss in range(len(subscapNames)):
    
    #Add axes to show image on
    imAx = plt.subplot2grid((8,8), (ss*2,0), rowspan = 2, colspan = 4)
    
    #Show image
    imAx.imshow(img, cmap = 'gray', origin = 'upper')
    
    #Turn off axes labels
    imAx.axis('off')
    
    #Display subscapularis line
    #Plot points
    imAx.scatter(subscapPoints[subscapNames[ss]]['X'], subscapPoints[subscapNames[ss]]['Y'],
                 s = 5, c = 'red')
    #Fit line
    subscapM,subscapC = np.polyfit(subscapPoints[subscapNames[ss]]['X'],
                                   subscapPoints[subscapNames[ss]]['Y'],
                                   1)
    #Plot fitted line
    imAx.plot(subscapPoints[subscapNames[ss]]['X'],
              subscapM * subscapPoints[subscapNames[ss]]['X'] + subscapC,
              c = 'red', lw = 1)
    #Extend line of action to min and max x-values
    #Get current axes limits
    retXLim = imAx.get_xlim()
    retYLim = imAx.get_ylim()
    #Plot lines
    imAx.plot(np.array((retXLim[0],np.min(subscapPoints[subscapNames[ss]]['X']))),
              subscapM * np.array((retXLim[0],np.min(subscapPoints[subscapNames[ss]]['X']))) + subscapC,
              c = 'red', lw = 1, ls = '--')
    imAx.plot(np.array((retXLim[1],np.max(subscapPoints[subscapNames[ss]]['X']))),
              subscapM * np.array((retXLim[1],np.max(subscapPoints[subscapNames[ss]]['X']))) + subscapC,
              c = 'red', lw = 1, ls = '--')
    #Reset axes (note that axes don't seem to be changed, but here just in case)
    imAx.set_xlim(retXLim)
    imAx.set_ylim(retYLim)
    
    #Add title
    imAx.set_title(muscleTitles[ss])
    
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
    # Ax1, Ay1, Ax2, Ay2, Bx1, By1, Bx2, By2 = gx1, gy1, gx2, gy2, sx1, sy1, sx2, sy2
    intX, intY = lineIntersect(gx1, gy1, gx2, gy2, sx1, sy1, sx2, sy2)
    imAx.scatter(intX, intY, c = 'blue', s = 5, zorder = 4)
    
    #Determine slope (i.e. negative reciprocal) of perpendicular line to glenoid plane
    mPerp = 1 / (-glenoidM)
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
    
    #Determine angle to rotate clockwise about intersection point based on glenoid
    #plane gradient
    #Negative gradient needs to be rotated clockwise by 90 - degrees
    if glenoidFlipM < 0:
        #Convert gradient to radians
        rotAng = np.radians(90) - np.arctan(glenoidM)
    # elif glenoidFlipM > 0:
    #     rotDeg = 
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
    lnAx.plot([0,maxLim/2*-1], [0,0], lw = 2, c = 'black') #x-axis
    lnAx.plot([0,0], [0,maxLim/2], lw = 2, c = 'black') #y-axis
    
    #Fit muscle line
    rotM, rotC = np.polyfit(rx, ry, 1)
    lnAx.plot(np.array(rx), rotM * np.array(rx) + rotC, c = 'red', lw = 1, zorder = 0)
    #Extend to origin
    lnAx.plot(np.array((0,np.min(rx))),
              rotM * np.array((0,np.min(rx))) + rotC,
              c = 'red', lw = 1, ls = '--', zorder = 0)
    
    #Turn off axes labels
    lnAx.axis('off')

    #Convert gradient of line to angle in degrees
    if rotM > 0:
        lineOfAction = 180 + np.degrees(np.arctan(rotM))
    elif m < 0:
        #Inferiorly directed line of action (i.e. < 180 degrees)
        lineOfAction = 180 - (np.degrees(np.arctan(rotM))*-1)
    elif m == 0:
        #180 degree (i.e. straight compression) line of action
        lineOfAction = 180
        
    #Print line of action in top left corner of axes
    lnAx.text(0.025, 0.9,
              'LoA = '+str(np.round(lineOfAction,2))+u'\u00b0',
              ha = 'left', va = 'center', transform = lnAx.transAxes)
        
    #### TODO: store line of action calculations

#Tight layout
plt.tight_layout()

