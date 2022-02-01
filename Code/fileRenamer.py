# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:48:06 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script reorganises the folder and file names across the different specimen,
    conditions, and testing configurations. The variability in the naming practices
    means this needs to be a little more manual than it could be. The process goes
    as follows:
        
        - Navigate to the top level folder for a specimen
        - Set the current condition to rename based on the folder name
        - Set the labels across the various condition categories
        - Run the script to rename folders and files for the current iteration
        - Re-run again from the top level of the specimens folder structure
        
    If using the 'run' function this needs to run from it's own directory, and 
    hence needs to navigate to the participants top level data folder.
    ***THIS IS THE BETTER OPTION***
    
    Note that it is best to get all the details for the folder, but then close
    Windows Explorer to avoid any shared open issues.
    
"""

# %% Import packages

import os
import shutil
from glob import glob
import re

# %% ----- CHANGE THESE -----

#Set the speciment ID name
specimenName = 'S9_l'

#Set the condition folder name and what to rename it to
conditionName = 'split50'
conditionRename = 'split50' #split50, split25_upper, split25_lower

#Set the position labels and what to rename to
positionNames = ['_0deg_', '_30deg_', '_60deg_',
                 '_90deg_', '_ABER_', '_Apreh_'] #need underscore as 0abd is in X0abd
#Sometimes position names are modified within the same section, and hence need to use a 2nd option sometimes
positionNames2 = ['_0abd_', '_30abd_', '_60abd_',
                 '_90abd_', '_ABER_', '_apreh_'] #need underscore as 0abd is in X0abd

#Set the loading names
loadNames = ['_0N', '_10N', '_20N', '_30N', '_40N'] #need underscore as 0N is in 10N

#Set the plane names
planeNames = [' SP', ' TP']
planeNames2 = ['_SP', '_TP']

#Set preface for scanning ID (there are sometimes multiples of certain scans to consider)
prefaceLabel = 'D0507_T'

# %% Set-up

# ----- Ensure these align to above names -----

#Set position, load and plane renaming (shouldn't change!)
positionRenames = ['abd0', 'abd30', 'abd60', 'abd90', 'ABER', 'APP']
loadRenames = ['0N', '10N', '20N', '30N', '40N']
planeRenames = ['SP', 'TP']

#Get starting directory
homeDir = os.getcwd()

#Navigate to the condition folder
os.chdir('..\\Data\\'+specimenName+'\\'+conditionName)

#Get all the folder names to work through
dirList = glob(os.getcwd()+'\\*\\')

# %% Rename files and folders

#Loop through directories
for dd in range(len(dirList)):
    
    #Navigate to current folder
    os.chdir(dirList[dd])
    
    #Set current label to scan
    scanLabel = dirList[dd].split('\\')[-2] #2nd last based on way glob splits
    
    #Identify current parameters
    
    #Position
    #Set position list checker to start in first list
    positionListChecker = 1
    currPosition = []
    for position in positionNames:
        if re.search(position, scanLabel):
            #Set the current position to this label
            currPosition.append(position)
    #Check if empty
    if len(currPosition) < 1:
        #Check in second set of position labels
        for position in positionNames2:
            if re.search(position, scanLabel):
                #Set the current position to this label
                currPosition.append(position)
        #Reset position checker variable to look in second list
        positionListChecker = 2
    #Check if still empty
    if len(currPosition) < 1:
        raise ValueError('No position identified in string')
    #Check for list longer than 1 and raise error, otherwise flatten to single string
    if len(currPosition) > 1:
        raise ValueError('More than 1 position identified in string')
    else:
        #Flatten list
        currPosition = currPosition[0]
        #Identify the index of the position in the names list
        if positionListChecker == 1:
            currPositionInd = positionNames.index(currPosition)
        elif positionListChecker == 2:
            currPositionInd = positionNames2.index(currPosition)
        #Set the rename position based on the index
        currPositionRename = positionRenames[currPositionInd]
    
    #Loading
    currLoad = []
    for load in loadNames:
        if re.search(load, scanLabel):
            #Set the current load to this label
            currLoad.append(load)
    #Check if empty
    if len(currLoad) < 1:
        raise ValueError('No load identified in string')
    #Check for list longer than 1 and raise error, otherwise flatten to single string
    if len(currLoad) > 1:
        raise ValueError('More than 1 load identified in string')
    else:
        #Flatten list
        currLoad = currLoad[0]
        #Identify the index of the position in the names list
        currLoadInd = loadNames.index(currLoad)
        #Set the rename position based on the index
        currLoadRename = loadRenames[currLoadInd]
    
    #Plane
    #Set plane list checker to start in first list
    planeListChecker = 1
    currPlane = []
    for plane in planeNames:
        if re.search(plane, scanLabel):
            #Set the current load to this label
            currPlane.append(plane)
    #Check if empty
    if len(currPlane) < 1:
        #Check in second set of position labels
        for plane in planeNames2:
            if re.search(plane, scanLabel):
                #Set the current position to this label
                currPlane.append(plane)
        #Reset position checker variable to look in second list
        planeListChecker = 2
    #Check if empty
    if len(currPlane) < 1:
        raise ValueError('No plane identified in string')
    #Check for list longer than 1 and raise error, otherwise flatten to single string
    if len(currPlane) > 1:
        raise ValueError('More than 1 plane identified in string')
    else:
        #Flatten list
        currPlane = currPlane[0]
        #Identify the index of the plane in the names list
        if planeListChecker == 1:
            currPlaneInd = planeNames.index(currPlane)
        elif planeListChecker == 2:
            currPlaneInd = planeNames2.index(currPlane)
        #Set the rename position based on the index
        currPlaneRename = planeRenames[currPlaneInd]
        
    #Identify the preface label
    #Find the position of the preface label in the string
    prefaceInd = scanLabel.find(prefaceLabel)
    #Check if not found
    if prefaceInd == -1:
        raise ValueError('Preface label not found')
    #Extract the actual preface label by taking from start index to next underscore
    prefaceRename = scanLabel[prefaceInd:-1].split('_')[0]+'_'+scanLabel[prefaceInd:-1].split('_')[1]
        
    #Create the entire renaming label
    renameLabel = currPositionRename+'_'+currLoadRename+'_'+currPlaneRename+'_'+prefaceRename
    
    #Find the tif files in the current folder
    imageList = glob('*.tif')
    
    #Loop through image list if necessary and rename accordingly
    for imageNo in range(len(imageList)):
        #Set new name
        newFileName = renameLabel+'_'+str(imageNo+1)+'.tif'
        #Rename file
        os.rename(imageList[imageNo], newFileName)
        
    #Navigate back up to condition directory
    os.chdir('..')
    
    #Lastly rename the folder to match the rename label
    #Double check if there is a folder already with this name
    if renameLabel in os.listdir():
        #Raise an error
        raise ValueError('Proposed renamed label already present in directory.')
    else:
        #Rename the folder
        os.rename(dirList[dd], renameLabel)
    
#Jump back up to the main directory
os.chdir('..')

#Rename the condition directory
#Double check if there is a folder already with this name
if conditionRename in os.listdir():
    #Raise an error
    raise ValueError('Proposed renamed condition already present in directory.')
else:
    #Rename the folder
    #Need to use shutil here rather than os to get around permission errors
    shutil.move(conditionName, conditionRename)
    
# %% Finish up
    
#Return to home directory
os.chdir(homeDir)
        
# %% ----- end of fileRenamer.py -----