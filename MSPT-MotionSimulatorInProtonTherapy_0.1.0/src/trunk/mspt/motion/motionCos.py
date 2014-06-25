########################################################################
#
# motionCos.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# October 2013
# 
#
#
# Copyright 2011-2014 Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
#
# This file is part of MSPT- Motion Simulator in Proton Therapy.
#
#    MSPT- Motion Simulator in Proton Therapy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MSPT- Motion Simulator in Proton Therapy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MSPT- Motion Simulator in Proton Therapy.  If not, see <http://www.gnu.org/licenses/>.
#
#
########################################################################

import numpy as np


# typeFloat = 'float32'
# #typeFloat = 'float64'

class MotionCos(object):
    '''Class modeling breathing motions as a cosine function. The model is based on models used by:
        
        * Lujan et al. 1999: "A method for incorporating organ motion due to breathing into 3D dose calculations"
        * George et al. 2005: "The application of the sinusoidal model to lung cancer patient respiratory motion"
        * Seco et al. 2009: "Breathing interplay effects during proton beam scanning: simulation and statistical analysis."
    
    The model of the breathing motion based on a cosine function was introduced by Lujan et al. and slightly modified later by George et al.\
    We use this last version that can be expressed as:
    
    .. math::
        
        f(t) = z0 + b cos( (2*Pi*t / tau)  - phi) \
        where z0 is the average position (it is set to 0),t the time, b the amplitude, and tau the motion period. \
        Lujan and Seco also added a phase phi that can be added in the cosine. 
        
    :param argsDict: Dictionary configuring the motion. Keys that must be in this dictionary are:
        
        * *'breathingPeriod'*: numpy array of 3 elements: [ x_period, y_period, z_period]
        * *'magnitude'*: numpy array of 3 elements: [ x_magnitude, y_magnitude, z_magnitude]
        * *'initialPhase'*: numpy array of 3 elements: [ x_phase, y_phase, z_phase]
    
    If no motion is one direction, set the 3 parameters for this direction (x,y or z) to 0.
    
    :param typeFloat: type of numpy float to use: 'float32' or 'float64'.
     
    '''
    def __init__(self,argsDict,typeFloat):
    
        '''
        
        '''
        
        requiredKeys = ['breathingPeriod' , 'magnitude', 'initialPhase']
        for key in requiredKeys:
            if key not in argsDict.keys():
                strErr = 'motionCos should be initialized with 3 arrays of parameters: %s , %s is missing.'%(requiredKeys,key)
                raise ValueError(strErr)
        
        self._typeFloat = typeFloat
        
        test=np.zeros((3),dtype=self._typeFloat)
        typetest= type(test)
        for idx,name in enumerate(['breathingPeriod' , 'initialPhase', 'magnitude']):
            item = argsDict[name]
            if type(item) != typetest:
                print item
                stringError = 'In init motion simulator, '+str(name) +' is not *NumPy* array'
                raise Exception(stringError)
            if len(np.shape(item)) != 1:
                print item
                stringError = 'In init motion simulator, '+str(name) +' has not only one row'
                raise Exception(stringError)
            if item.size != test.size:
                print item
                stringError = 'In init motion simulator, '+str(name) +' has not a size of 3.'
                raise Exception(stringError)
            if item.dtype != test.dtype:
                print item
                stringError = 'In init motion simulator, '+str(name) +' is not of type *Float*.'
                raise Exception(stringError)
        
        tabDirection = ['X','Y','Z']
        for idx,value in np.ndenumerate(argsDict['breathingPeriod']):
            if value < 0:
                strError = "In init motion Cos, %s direction has a breathing period < 0!"
                raise ValueError(strError)
        self._breathingPeriod = argsDict['breathingPeriod']
        self._indDirection = np.where(self._breathingPeriod > 0) #Indices for motion directions
        self._magnitude = argsDict['magnitude']
        self._phase = argsDict['initialPhase']
        
    
     
    def __str__(self):
        strValue = "Motion: p0 + ampl. x cos( 2*PI*t / tau  - phi ):\n\tdirections:[X, Y, Z]\n\t-Breathing period: %s sec\n\t-Amplitude: %s cm\n\t-Initial Phase: %s rad\n"%(str(self._breathingPeriod),str(self._magnitude),str(self._phase))
        return strValue
   
    
    def getDisplacementVectorAtTime(self,timer):
        '''Computes the displacement vector according to the equation:
        
        .. math::
        
            f(t) = z0 + b cos( 2*Pi*t / tau  + phi) 
        
        :param timer: time in sec. 
        
        :returns: Displacement vector : numpy array with 3 elements.
        '''
        if timer < 0 :
            raise ValueError("Time < 0 in get displacement vector")
        vec = np.zeros((3),dtype=self._typeFloat)
        vec[self._indDirection] = self._magnitude[self._indDirection] * np.cos( ((timer * 2 * np.pi ) / self._breathingPeriod[self._indDirection]) + self._phase[self._indDirection] )
        
        return vec
    
    
        