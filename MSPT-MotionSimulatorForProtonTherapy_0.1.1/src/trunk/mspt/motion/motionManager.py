########################################################################
#
# motionSuperClass.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# October 2013
# 
#
# Copyright 2011-2014 Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
#
# This file is part of MSPT- Motion Simulator for Proton Therapy.
#
#    MSPT- Motion Simulator for Proton Therapy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MSPT- Motion Simulator for Proton Therapy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MSPT- Motion Simulator for Proton Therapy.  If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

import numpy as np
from motionCos import MotionCos
from motionNoisyCos import MotionNoisyCos
from motionBreathHiccup import MotionBreathHiccup
from motionBreathCough import MotionBreathCough
from motionRecorded import MotionRecorded
import os


motionTypes = ['motionNoisyCos','motionCos','motionBreathHiccup','motionBreathCough','motionRecorded']

class Motion(object):
    '''Class Motion: used to define motions. They can be of different types: 
    
        * *'motionCos'*: Motion defined as a cosine wave along X,Y,Z axis.
        * *'motionNoisyCos'*: Similar to 'motionCos' with random noise.
        * *'motionBreathHiccup'*: 'NoisyCos' plus hiccup. Hiccup is defined by a starting and ending date and a periodical motion of random amplitude is added to the initial noisy cosine.
        * *'motionBreathCough'*: 'NoisyCos' plus cough. Cough is defined by a single random motion occuring at a random time.
        * *'motionRecorded'*: A list of dates and X,Y,Z displacement vectors
        
    .. note::
    
        The user could create new motions. To do so, in motionManager.py add the name of the created class (the first letter must be in lower case) in the global variable 'motionTypes'.\
        The new class must implement 'getDisplacementVectorAtTime(timer)' to define the displacement vector at time = timer. One could \
        decide to return a 3D matrix containing local displacement vectors (In such case patient.patient.py and scanning.deliverySimulation_Static_Dynamic.py\
        should be adapted for deformable transformations. The 'history' attribute of MotionManager should also be adapted). It must also implement '__str__' to provide information about the motion.
        
    :param argsDict: dictionary of arguments used to initialize the types of motion. The key 'typeMotion' must be present. It it equal to \
    one of theses types: 'motionNoisyCos','motionCos','motionBreathHiccup','motionBreathCough'
    :param typeFloat: type of numpy float to use: 'float32' or 'float64'.
    
        .. note: 
    
        The units for the stdandard deviations and the averages must be in cm and the time in sec.
   
    '''
    def __init__(self,argsDict,typeFloat):
        '''

        '''
        if 'typeMotion' not in argsDict.keys():
            raise ValueError ("Motion initialization impossible: no 'typeMotion' value defined")
        
        if argsDict['typeMotion'] in motionTypes:
            code  = 'self._motion = %s(argsDict,typeFloat)'%(argsDict['typeMotion'][0].capitalize()+argsDict['typeMotion'][1:])
            exec code
        else:
            strErr = "Unknown motion type: %s"%argsDict['typeMotion']
        self._history = {'time':[],'x':[],'y':[],'z':[]}
        self._typeFloat = typeFloat
     
    def __str__(self):
        '''
        :returns: The string information about the current motion.
        
        
        '''
        return str(self._motion)
   
    
    def getDisplacementVectorAtTime(self,timer):
        '''
        :param timer: Time at when the displacement vector must be computed.    
        
        :returns: Displacement vector defined in the current motion.
        
        '''
        vec = self._motion.getDisplacementVectorAtTime(timer)
        vecTest = np.zeros((3),dtype=self._typeFloat)
        if timer not in self._history['time'] and vec.shape == vecTest.shape:
            self._history['time'].append(timer)
            self._history['x'].append(vec[0])
            self._history['y'].append(vec[1])
            self._history['z'].append(vec[2])
        elif timer not in self._history['time']:
            print "Warning: History of motion won't be saved. Please adapt the code for the history in motion.motionManager."
        return vec
    
    def getHistory(self):
        '''Get the history of the motions. 
        
        :returns: A dictionary with keys: 'time', 'x', 'y' and 'z'. Each value of these keys is a list.
        
        '''
        return self._history

    
    
    
    