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

# %% Import image

#Create image object from files
img = Image.open('UDF_S19-009R_specimen1_D0910_T1435_condNP_0abd_40N_SP.tif')
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

#The Y-axis coordinates are inverted due to the way image coordinates are imported
#To fix this so that the XY axes originates in the lower left corner we invert
#them relative to the length of the axis (i.e. image height). We do this in getting
#the y-values below
    ##### NOT DOING THIS ANYMORE...
        
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
        if currItem[0].startswith('hh'):
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

# %% Visualise

#Visualise points on image
fig, ax = plt.subplots()
#Plot image
# plt.imshow(img, cmap = 'gray')
plt.imshow(np.flipud(img), cmap = 'gray', origin = 'lower')
# plt.imshow(imgDig, cmap = 'gray')
# plt.imshow(np.flipud(imgDig), cmap = 'gray', origin = 'lower')

##### seems necessary to invert the points but not the image
##### it seems the points come across as you would think from an XY perspective
##### from Mimics - in that going up = increasing Y, but the image doesn't 
##### import that way --- it's probably just a visualisation thing then...
    ##### the points overlaid on the image look slightly off...
    ##### I think it's OK and it's just a weird visual hitch...when you zoom in
    ##### they look accurate

#Plot subscapularis points
ssLabels = ['ss1','ss2','ss3','ss4']
for pp in range(len(ssLabels)):
    pts = df_items.loc[df_items['pointType'] == ssLabels[pp],['X','Y']].to_numpy()
    plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'red')

#Plot humeral head points
pts = df_items.loc[df_items['pointType'] == 'hh',['X','Y']].to_numpy()
plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'green')

#Plot glenoid plane points
pts = df_items.loc[df_items['pointType'] == 'gp',['X','Y']].to_numpy()
plt.scatter(pts[:,0], pts[:,1], s = 10, c = 'yellow')

#Plot bead circle
pts = df_items.loc[df_items['pointType'] == 'phantom',['X','Y','radius']].to_numpy()
ax.add_artist(plt.Circle((pts[0][0], pts[0][1]), pts[0][2],
                         edgecolor = 'blue', facecolor = 'none'))

#Plot humeral head circle and centre
plt.scatter(xc, yc, s = 500, c = 'green', marker = '+')
ax.add_artist(plt.Circle((xc, yc), R, edgecolor = 'green', facecolor = 'none'))
