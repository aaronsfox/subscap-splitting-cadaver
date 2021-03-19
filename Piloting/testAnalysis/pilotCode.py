# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:03:22 2020

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script does a test run of processing a digitised image of the subscapularis 
    lines of action. Notes are provided in comments throughout.
    
"""

# %% Import packages

from PIL import Image
from matplotlib import pyplot as plt
import csv
import numpy as np
import pandas as pd
import math

from skimage import io

# %% Import image

#Create image object from files
# I = plt.imread('UDF_S19-009R_specimen1_D0910_T1435_condNP_0abd_40N_SP.tif')
img = io.imread('UDF_S19-009R_specimen1_D0910_T1435_condNP_0abd_40N_SP.tif')
# img = Image.open('UDF_S19-009R_specimen1_D0910_T1435_condNP_0abd_40N_SP.tif')
# imgArray = np.array(img)
# imgDig = Image.open('sampleDigitisation.bmp')

#The image in this case needs to be rotated 180 degrees to match Mimics
img = img.rotate(180)

##### NOTE: might need to check whether the image took any cropping in Mimics.
##### This might cause mismatch with the overlay. 
    ##### Importing the digitised Mimics bitmap seems to be identical after the
    ##### rotation, so this doesn't seem to be a problem.

#View image
plt.imshow(np.flipud(img), cmap = 'gray', origin = 'lower')
# plt.imshow(np.flipud(imgDig), cmap = 'gray', origin = 'lower')

# %% Import digitised items

#Read in raw text items
itemList = []
with open('digitisationExport.txt') as txtFile:                                                                                          
	itemReader = csv.reader(txtFile, delimiter = '\t')
	for item in itemReader:
		itemList.append(item)

#Set-up dictionary to store point data in
itemDict = {'pointName': [], 'pointType': [], 'X': [], 'Y': [], 'radius': []}

##### NOTE: in this case XZ represents the XY points and hence the item list
##### indices reflect this
        
#Search through the item list for the various items
for ii in range(len(itemList)):
    
    #Extract current item
    currItem = itemList[ii]
    
    #First check if not empty
    if len(currItem) > 0:
    
        #Check the first item (i.e. item name) against criteria to extract
        #Top subscapularis
        if currItem[0].startswith('ss1'):
            #Set point type
            itemDict['pointType'].append('ss1')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Middle-top subscapularis
        if currItem[0].startswith('ss2'):
            #Set point type
            itemDict['pointType'].append('ss2')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Middle-bottom subscapularis
        if currItem[0].startswith('ss3'):
            #Set point type
            itemDict['pointType'].append('ss3')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Bottom subscapularis
        if currItem[0].startswith('ss4'):
            #Set point type
            itemDict['pointType'].append('ss4')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Humeral head
        if currItem[0].startswith('hh') and len(currItem) == 5:
            #Set point type
            itemDict['pointType'].append('hh')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Glenoid plane
        if currItem[0].startswith('gp'):
            #Set point type
            itemDict['pointType'].append('gp')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius as NaN
            itemDict['radius'].append(np.nan)
        #Phantom bead
        if currItem[0].startswith('phantom'):
            #Set point type
            itemDict['pointType'].append('phantom')
            #Extract point name
            itemDict['pointName'].append(currItem[0].strip())
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius
            itemDict['radius'].append(float(currItem[4].strip()))
        #Humeral head circle (3pt method)
        if currItem[0].startswith('hh3pt'):
            #Set point type
            itemDict['pointType'].append('hhCircle')
            #Extract point name
            itemDict['pointName'].append('hhCircle')
            #Extract X coordinate
            itemDict['X'].append(float(currItem[1].strip()))
            #Extract Y coordinate
            itemDict['Y'].append(float(currItem[3].strip()))
            #Set radius
            itemDict['radius'].append(float(currItem[4].strip()))

#Convert dictionary to dataframe
df_items = pd.DataFrame.from_dict(itemDict)

# %% Processing

#Fit a circle to the humeral head points

##### TODO: clean up function...

##### TODO: probably need more points on humerus to round it out...

#Get humeral head points
hhPts = df_items.loc[df_items['pointType'] == 'hh',['X','Y']].to_numpy()
x = hhPts[:,0]
y = hhPts[:,1]

#Calculate means
x_m = np.mean(x)
y_m = np.mean(y)

#Calculation of the reduced coordinates
u = x - x_m
v = y - y_m

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
xc = x_m + uc
yc = y_m + vc

#Calculate distance to centre and radius
Ri = np.sqrt((x-xc)**2 + (y-yc)**2)
R = np.mean(Ri)
residuals = sum((Ri-R)**2)

#### NOTE: I prefer the circular fitting in Mimics manually, rather than by points
#### provides an opportunity to fit it yourself with experimentation

# %% Visualise

#Visualise points on image
fig, ax = plt.subplots()
#Plot image
# plt.imshow(img, cmap = 'gray')
plt.imshow(np.flipud(img), cmap = 'gray', origin = 'lower')
# plt.imshow(imgDig, cmap = 'gray')
# plt.imshow(np.flipud(imgDig), cmap = 'gray', origin = 'lower')

# #Set axes if not loading image
# ax.set_xlim([0,img.size[0]])
# ax.set_ylim([0,img.size[1]])

##### seems necessary to invert the points but not the image
##### it seems the points come across as you would think from an XY perspective
##### from Mimics - in that going up = increasing Y, but the image doesn't 
##### import that way --- it's probably just a visualisation thing then...
    ##### the points overlaid on the image look slightly off...
    ##### I think it's OK and it's just a weird visual hitch...when you zoom in
    ##### they look accurate

#Plot subscapularis points and lines
#Set point labels
ssLabels = ['ss1','ss2','ss3','ss4']
#Identify point to extend lines of action to
#Currently half-way between image edge and humeral head centre
extendX = df_items.loc[df_items['pointType'] == 'hhCircle',['X']].values[0][0] / 2
for pp in range(len(ssLabels)):
    #Get points
    pts = df_items.loc[df_items['pointType'] == ssLabels[pp],['X','Y']].to_numpy()
    #Plot points
    plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'red')
    #Fit and plot line
    m,c = np.polyfit(pts[:,0], pts[:,1], 1)
    plt.plot(pts[:,0], m * pts[:,0] + c, c = 'red', lw = 1)
    #Extend line of action
    plt.plot(np.array((extendX,np.min(pts[:,0]))),
             m * np.array((extendX,np.min(pts[:,0]))) + c,
             c = 'red', lw = 1, ls = '--')

# #Plot humeral head points
# pts = df_items.loc[df_items['pointType'] == 'hh',['X','Y']].to_numpy()
# plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'green')
    
#Plot humeral head circle
pts = df_items.loc[df_items['pointType'] == 'hhCircle',['X','Y','radius']].to_numpy()
ax.add_artist(plt.Circle((pts[0][0], pts[0][1]), pts[0][2],
                         edgecolor = 'green', facecolor = 'none'))
#Include the humeral head centre
plt.scatter(pts[0][0], pts[0][1], s = 500, c = 'green', marker = '+')

#Plot glenoid plane points and fit line
pts = df_items.loc[df_items['pointType'] == 'gp',['X','Y']].to_numpy()
plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'yellow')
m,c = np.polyfit(pts[:,0], pts[:,1], 1)
plt.plot(pts[:,0], m * pts[:,0] + c, c = 'yellow', lw = 1)

# #Plot bead circle
# pts = df_items.loc[df_items['pointType'] == 'phantom',['X','Y','radius']].to_numpy()
# ax.add_artist(plt.Circle((pts[0][0], pts[0][1]), pts[0][2],
#                          edgecolor = 'blue', facecolor = 'none'))

# #Plot multi-point calculated humeral head circle and centre
# plt.scatter(xc, yc, s = 500, c = 'green', marker = '+')
# ax.add_artist(plt.Circle((xc, yc), R, edgecolor = 'green', facecolor = 'none'))

##### NOTE: 3pt fits pretty similarly to this example - probably easier option

# #Equal axes
# plt.gca().set_aspect('equal', adjustable = 'box')

# %% Moment arm calculations

#Calculate the scaling factor to make the bead radius = 3
beadImgRadius = df_items.loc[df_items['pointType'] == 'phantom',['radius']].to_numpy()[0][0]
beadScale = beadImgRadius / 3

##### TODO: clean up function...

#Visualise moment arm line of one of the subscapularis lines
fig, ax = plt.subplots()
plt.imshow(np.flipud(img), cmap = 'gray', origin = 'lower')
#Humeral head
pts = df_items.loc[df_items['pointType'] == 'hhCircle',['X','Y','radius']].to_numpy()
plt.scatter(pts[0][0], pts[0][1], s = 500, c = 'green', marker = '+')
ax.add_artist(plt.Circle((pts[0][0], pts[0][1]), pts[0][2],
                         edgecolor = 'green', facecolor = 'none'))
#Subscapularis line
pp = 3
pts = df_items.loc[df_items['pointType'] == ssLabels[pp],['X','Y']].to_numpy()
plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'red')
m,c = np.polyfit(pts[:,0], pts[:,1], 1)
plt.plot(pts[:,0], m * pts[:,0] + c, c = 'red', lw = 1)
plt.plot(np.array((extendX,np.min(pts[:,0]))),
         m * np.array((extendX,np.min(pts[:,0]))) + c,
         c = 'red', lw = 1, ls = '--')

#Calculate the point on the extended subscapularis line that is intercepted
#by a perpendicular line from the humerus

#Set the points at the start and the end of the subscapularis line
x1 = extendX
y1 = m * extendX + c
x2 = max(pts[:,0])
y2 = m * x2 + c

#Set the humeral head point
x3 = df_items.loc[df_items['pointType'] == 'hhCircle',['X']].values[0][0]
y3 = df_items.loc[df_items['pointType'] == 'hhCircle',['Y']].values[0][0]

#Calculate intersecting point with perpendicular line
#SEE: https://stackoverflow.com/questions/1811549/perpendicular-on-a-line-from-a-given-point
k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)**2 + (x2-x1)**2)
x4 = x3 - k * (y2-y1)
y4 = y3 + k * (x2-x1)

#Plot the perpendicular line between the points
#This is effectively the moment arm
plt.plot(np.array((x3,x4)),np.array((y3,y4)),
         c = 'yellow', lw = 1, ls = '--')

#Calculate moment arm distance
##### NOTE: this isn't currently scaled by the bead size, as I believe
##### the units might already be accurate --- but could be incorporated.
##### In addition, the direction of the moment arm (+ve ve. -ve) should
##### also be considered.
ma = (math.sqrt(((x3 - x4)**2) + ((y3 - y4)**2))) / beadScale

##### The units do look incorrect, as the radius of the bead is ~30, rather
##### than the ~3 it should be. The moment arm calculation above divided by
##### 10 also lines up with the Ackland et al. study

# %% Calculate lines of action

##### TODO: for lines of action we should simply rotate the 
##### coordinates/lines so that the normal of the glenoid plane
##### is horizontal, which will make the angle calculation easier...
