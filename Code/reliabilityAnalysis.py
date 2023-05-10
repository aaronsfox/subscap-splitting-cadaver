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

# %% Set-up

# #Set matplotlib parameters
# from matplotlib import rcParams
# # rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = 'Arial'
# rcParams['font.weight'] = 'bold'
# rcParams['axes.labelsize'] = 12
# rcParams['axes.titlesize'] = 16
# rcParams['axes.linewidth'] = 1.5
# rcParams['axes.labelweight'] = 'bold'
# rcParams['legend.fontsize'] = 10
# rcParams['xtick.major.width'] = 1.5
# rcParams['ytick.major.width'] = 1.5
# rcParams['legend.framealpha'] = 0.0
# rcParams['savefig.dpi'] = 300
# rcParams['savefig.format'] = 'pdf'

#Set specimen
specimenNames = ['S1_r', 'S2_r', 'S3_r', 'S4_l', 'S5_l', 'S6_r', 'S7_r', 'S8_l', 'S9_l']

#Set conditions
conditionNames = ['split50', 'split25_upper', 'split25_lower']

#Set loads
loadNames = ['0N', '20N', '40N']

#Set muscle lines
muscleNames = ['ss1', 'ss2', 'ss3', 'ss4']

#Set arm positions
posNames = ['0abd', '90abd', 'ABER', 'APP']

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

# %% ----- End of reliabilityAnalysis.py ----- %% #