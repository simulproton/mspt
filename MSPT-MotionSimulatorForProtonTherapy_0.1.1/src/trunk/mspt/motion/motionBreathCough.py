########################################################################
#
# motionBreathCough.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# 19 February 2014
# 
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
#
########################################################################

import numpy as np
import motionNoisyCos
from motionNoisyCos import MotionNoisyCos


class MotionBreathCough(MotionNoisyCos):
    '''
    The class MotionBreathCough, inherit from MotionNoisyCos. In addition to noise added to the motion, this class add a cough to the motion.\
    The cough is modeled as a single event happening at a random date. The amplitude of a cough is the initial amplitude (wihtout cough)\
    minus a random value defined by a normal distribution (mean = 'magnitude'/2 and stdev = 'variationsMag'/2). The reason why we are using a "minus" instead of\
    a "plus" is that usually when someone coughs his rib cage is compressed.
    
    :param argsDict: Dictionary configuring the motion. Keys that must be in this dictionary are:
        
        * *'breathingPeriod'*: numpy array of 3 elements: [ x_period, y_period, z_period]
        * *'magnitude'*: numpy array of 3 elements: [ x_magnitude, y_magnitude, z_magnitude]
        * *'initialPhase'*: numpy array of 3 elements: [ x_phase, y_phase, z_phase]
        * *'variationsMag'*: numpy array of 3 elements representing the standard deviation of the normal distribution applied to\
        the motion amplitude of each coordinate. We consider that the mean of the distribution is the value defined in 'magnitude'.\
        This represents the noise applied to the motion magnitudes.
        * *'variationsPeriod'*: numpy array of 3 elements representing the standard deviation of the normal distribution applied to\
        the length of the period for each coordinate. We consider that the mean of the distribution is the value defined in 'breathingPeriod'.\
        This represents the noise applied to the length of the periods.
    
    :param typeFloat: type of numpy float to use: 'float32' or 'float64'. 
    
    
        
    '''
    def __init__(self,argsDict,typeFloat):
    
        '''
        
        '''
        MotionNoisyCos.__init__(self,argsDict,typeFloat)
        
        self._coughTime = np.random.rand()*30 # Time when happens the coughing
        self._coughOccured = False
        self._CoughFunc = [motionNoisyCos.buildNormDistrFunction(self._magnitude[i]/2.0 , self._stdevMag[i]/2.0) for i in range (3)]
        
     
    def __str__(self):
        strValue = "Motion: p0 + ampl. x cos( 2*PI*t / tau  - phi ):\n\tdirections:[X, Y, Z]\n\t-Breathing period: %s sec\n\t-Amplitude: %s cm\n\t-Initial Phase: %s rad\n"%(str(self._breathingPeriod),str(self._magnitude),str(self._phase))
        strVariations = "\tMagnitude PDF type:%s , stdev = %s\n\tPeriod PDF type:%s , stdev = %s\n"%(self._distrTypMag,str(self._stdevMag),self._distrTypPeriod,str(self._stdevPeriod))
        strCough = "Coughing should happen after around %f sec\n"%self._coughTime
        return strValue+strVariations + strCough
   
    
    def getDisplacementVectorAtTime(self,timer):
        '''Computes the displacement vector. It uses the getDisplacementVectorAtTime() function of MotionNoisyCos.
        
        We add a "Cough" when the random date has just occured. The random amplitude is then subtracted from the initial motion amplitude. \
        The amplitude of the cough is random using a Gaussian pdf whose mean amplitude is half the amplitude used for the noise\
        and the stdev is half the stdev used for noise in the motion for each coordinate.

        :param  timer: Time in sec. 
        
        :returns: Displacement vector : numpy array with 3 elements. 

        '''
        vec = super(  MotionBreathCough,self).getDisplacementVectorAtTime(timer)
        if timer >= self._coughTime and not self._coughOccured:
            self._coughOccured = True
            print "Cough at %f sec."%timer
            return updateVecDispl( self._CoughFunc,vec)
        else: 
             return vec
        
    

def updateVecDispl(CoughFunc,vec):
    '''
    Update a displacement vector, by subtraction the cough amplitude to each coordinate.\
    We consider the cough as reducing the amplitude.
    
    :param CoughFunc: list of normal distribution functions
    :param vec: Initial displacement vector
    
    :returns: Displacement vector : numpy array with 3 elements. 
    
    '''
    for i in range(3):
        vec[i] = vec[i] - CoughFunc[i]()
    return vec
    
        