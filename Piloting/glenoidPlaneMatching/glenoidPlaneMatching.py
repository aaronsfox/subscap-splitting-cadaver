# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 16:02:24 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Pilot script to test out whether we can use a single glenoid plane digitisation
    and map glenoid planes on other images based on digitised landmarks. The general
    premise is to identify the translation/rotation matrix that will align the points,
    and then apply that to the glenoid plane from the original image to map it on
    the new image.
    
"""

# %% Import packages

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from skimage import io

# %% Define functions

def procrustes(X, Y, scaling = True, reflection = 'best'):
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

def applyTransform(Y, tform):
    """
    A follow on to apply the transform extracted from the procrustes function
    to a new set of points

    Inputs:
    ------------
    Y    
        input of XY coordinates to be transformed
        rows = points; columns = X,Y

    tform 
        the transform dictionary extracted from the procrustes function

    Outputs
    ------------
    Z
        the matrix of transformed values

    """

    # n,m = X.shape
    ny,my = Y.shape

    # muX = X.mean(0)
    muY = Y.mean(0)

    # X0 = X - muX
    Y0 = Y - muY

    # ssX = (X0**2.).sum()
    ssY = (Y0**2.).sum()

    # centred Frobenius norm
    # normX = np.sqrt(ssX)
    normY = np.sqrt(ssY)

    # scale to equal (unit) norm
    # X0 /= normX
    Y0 /= normY

    # # optimum rotation matrix of Y
    # A = np.dot(X0.T, Y0)
    # U,s,Vt = np.linalg.svd(A,full_matrices=False)
    # V = Vt.T
    # T = np.dot(V, U.T)

    

    traceTA = s.sum()

    Z = normY*np.dot(Y0, T) + muX

    # transformation matrix
    if my < m:
        T = T[:my,:]
    c = muX - b*np.dot(muY, T)
    
    #transformation values 
    tform = {'rotation':T, 'scale':b, 'translation':c}
   
    return d, Z, tform

# %% Scapular plane test

### TODO: Considerations for upper coordinate system of images?

#Load images
img1_SP = io.imread('img1_SP_abd0_0N.tif')
img2_SP = io.imread('img2_SP_abd60_0N.tif')
img3_SP = io.imread('img3_SP_abd90_40N.tif')

#Visualise images
fig, ax = plt.subplots(ncols = 3, figsize = (12,4))
ax[0].imshow(img1_SP, cmap = 'gray', origin = 'upper')
ax[1].imshow(img2_SP, cmap = 'gray', origin = 'upper')
ax[2].imshow(img3_SP, cmap = 'gray', origin = 'upper')

#Import and visualise glenoid plane
gp = pd.read_csv('img1_SP_gp.csv')
ax[0].scatter(gp['X'], gp['Y'], s = 10, c = 'yellow')
m,c = np.polyfit(gp['X'], gp['Y'], 1)
ax[0].plot(gp['X'], m * gp['X'] + c, c = 'yellow', lw = 1)

#Import and visualise points
#Image 1
img1_SP_pts = pd.read_csv('img1_SP_pts.csv')
ax[0].scatter(img1_SP_pts['X'].to_numpy(), img1_SP_pts['Y'].to_numpy(),
              c = 'red', s = 10)
#Image 2
img2_SP_pts = pd.read_csv('img2_SP_pts.csv')
ax[1].scatter(img2_SP_pts['X'].to_numpy(), img2_SP_pts['Y'].to_numpy(),
              c = 'green', s = 10)
#Image 3
img3_SP_pts = pd.read_csv('img3_SP_pts.csv')
ax[2].scatter(img3_SP_pts['X'].to_numpy(), img3_SP_pts['Y'].to_numpy(),
              c = 'magenta', s = 10)

#Apply the procrustes algorithm to get the transforms and points for img2 and img3
#Image 2
X = np.transpose(np.array([img1_SP_pts['X'].to_numpy(),
                           img1_SP_pts['Y'].to_numpy()]))
Y = np.transpose(np.array([img2_SP_pts['X'].to_numpy(),
                           img2_SP_pts['Y'].to_numpy()]))
img2_resErr, img2_newPts, img2_tform = procrustes(X, Y, scaling = False)
#Plot new image 2 points on image 1 to look at error
ax[0].scatter(img2_newPts[:,0], img2_newPts[:,1], c = 'green', s = 5)
#Image 3
Y = np.transpose(np.array([img3_SP_pts['X'].to_numpy(),
                           img3_SP_pts['Y'].to_numpy()]))
img3_resErr, img3_newPts, img3_tform = procrustes(X, Y, scaling = False)
#Plot new image 2 points on image 1 to look at error
ax[0].scatter(img3_newPts[:,0], img3_newPts[:,1], c = 'magenta', s = 5)

#Apply the transform to the glenoid plane for each of the two images to get the
#new glenoid planes for image2 and image3
#Get gp points in array
gpPts = np.transpose(np.array([gp['X'], gp['Y']]))
#Image 2 transform
img2_gpPts = np.zeros(gpPts.shape)
for pp in range(gpPts.shape[0]):
    img2_gpPts[pp] = img2_tform['rotation'].dot(gpPts[pp]) - img2_tform['translation']
#Image 3 transform
img3_gpPts = np.zeros(gpPts.shape)
for pp in range(gpPts.shape[0]):
    img3_gpPts[pp] = img3_tform['rotation'].dot(gpPts[pp]) - img3_tform['translation']

#Visualise new gp points on images
#Image 2
ax[1].scatter(img2_gpPts[:,0], img2_gpPts[:,1], s = 10, c = 'yellow')
m,c = np.polyfit(img2_gpPts[:,0], img2_gpPts[:,1], 1)
ax[1].plot(img2_gpPts[:,0], m * img2_gpPts[:,0] + c, c = 'yellow', lw = 1)
#Image 3
ax[2].scatter(img3_gpPts[:,0], img3_gpPts[:,1], s = 10, c = 'yellow')
m,c = np.polyfit(img3_gpPts[:,0], img3_gpPts[:,1], 1)
ax[2].plot(img3_gpPts[:,0], m * img3_gpPts[:,0] + c, c = 'yellow', lw = 1)

#Add axes titles
ax[0].set_title('Digitised Glenoid Plane & Aligned Points',
                fontsize = 12, fontweight = 'bold', fontfamily = 'Arial')
ax[1].set_title('Digitised Points & Transformed Glenoid',
                fontsize = 12, fontweight = 'bold', fontfamily = 'Arial')
ax[2].set_title('Digitised Points & Transformed Glenoid',
                fontsize = 12, fontweight = 'bold', fontfamily = 'Arial')

#Turn off axes
ax[0].axis('off')
ax[1].axis('off')
ax[2].axis('off')

#Tight layout
plt.tight_layout()

#Save image
plt.savefig('estimatedTransformedGlenoidPlanes_Scapular.png',
            format = 'png', dpi = 300)

#Close figure
plt.close()

# %% Transverse plane test


