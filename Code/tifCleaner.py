# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:42:07 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Script that cleans out the 'CAM0' tif files that are black
    
    MUST CHECK WHAT THE PREFACE OF THE BLANK FILE IS AND CHANGE!
    
    Also ensure that you're in the condition directory
    
"""

# %% Import packages

import os
from glob import glob

# %% Delete unwanted files

#Get the directories in the current folder
dirList = glob(os.getcwd()+'\\*\\')

#Set blank preface
blankPreface = 'CAM0'

#Loop through directories
for dd in range(len(dirList)):
    
    #Navigate to current directory
    os.chdir(dirList[dd])
    
    #Get file list
    imageList = glob('*.tif')
    
    #Loop through and delete if blank preface files present
    for ii in range(len(imageList)):
        if imageList[ii].startswith(blankPreface):
            #Delete it
            os.remove(imageList[ii])
            
    #Navigate back up a directory
    os.chdir('..')
        
# %% ----- end of tifCleaner.py -----
    
