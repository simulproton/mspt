########################################################################
#
# motionBreathHiccup.py
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

nbMaxHiccup = 10

class MotionBreathHiccup(MotionNoisyCos):
    '''
    The class MotionBreathHiccup, inherit from MotionNoisyCos. In addition to noise added to the motion, \
    this class adds some hiccups to the motion.
    
    The hiccup is modeled as a periodical event starting at a random date, having a random period, \
    a random number of period (for ease we set this number between 1 and 10: the value \
    can be changed in motion.motionBreathHiccup with the global variable nbMaxHiccup) which determines the \
    ending date of the hiccup. The amplitude of a hiccup is the sum of the initial amplitude (wihtout hiccup)\
    plus a random value defined by a normal distribution (mean = 'magnitude'/2 and stdev = 'variationsMag'/2).
    
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
        
        self._startHiccupTime = np.random.rand()*10 # Hiccup starts between 0 and 10 seconds after the beginning of the delivery
        self._hiccupPeriod = np.random.rand()*3 # Hiccup period between 0 and 3 seconds
        self._endHiccupTime = self._startHiccupTime + np.random.randint(1,nbMaxHiccup)*self._hiccupPeriod # Hiccup end: start time + nb of periods * period duration 
        self._HiccupFunc = [motionNoisyCos.buildNormDistrFunction(self._magnitude[i]/2.0 , self._stdevMag[i]/2.0) for i in range (3)]
        self._lastHiccup = None

    
     
    def __str__(self):
        strValue = "Motion: p0 + ampl. x cos( 2*PI*t / tau  - phi ):\n\tdirections:[X, Y, Z]\n\t-Breathing period: %s sec\n\t-Amplitude: %s cm\n\t-Initial Phase: %s rad\n"%(str(self._breathingPeriod),str(self._magnitude),str(self._phase))
        strVariations = "\tMagnitude PDF type:%s , stdev = %s\n\tPeriod PDF type:%s , stdev = %s\n"%(self._distrTypMag,str(self._stdevMag),self._distrTypPeriod,str(self._stdevPeriod))
        strHiccup = "Hiccup should happen after around %f sec, every %f sec, until %f sec\n"%(self._startHiccupTime,self._hiccupPeriod,self._endHiccupTime)
        return strValue+strVariations + strHiccup
   
    
    def getDisplacementVectorAtTime(self,timer):
        '''Computes the displacement vector. It uses the getDisplacementVectorAtTime() function of MotionNoisyCos.\
        We add a "Hiccup" when using a random starting date, a random hiccup period duration and a random number of periods.\
        The amplitude of the hiccup is random using a Gaussian pdf whose mean amplitude is half the amplitude used for the noise\
        and the stdev is half the stdev used for noise in the motion for each coordinate.

        :param  timer: Time in sec. 
        
        :returns: Displacement vector : numpy array with 3 elements. 

        '''
        vec = super( MotionBreathHiccup,self).getDisplacementVectorAtTime(timer)
        if self._endHiccupTime >= timer >= self._startHiccupTime:
             if self._lastHiccup is not None:
                deltaT = timer - self._lastHiccup
                if deltaT < self._hiccupPeriod:
                    return vec
                else: 
                    self._lastHiccup = timer
                    print "Hiccup at %f sec."%timer
                    return updateVecDispl(self._HiccupFunc,vec)
             else:
                self._lastHiccup = timer
                print "Hiccup at %f sec."%timer
                return updateVecDispl(self._HiccupFunc,vec)
        else: 
             return vec
        
    

def updateVecDispl(HiccupFunc,vec):
    '''Update a displacement vector, by adding the hiccup amplitude to each coordinate.
    
    :param HiccupFunc: list of normal distribution functions
    :param vec: Initial displacement vector
    
    :returns: Displacement vector : numpy array with 3 elements. 
    
    '''
    for i in range(3):
        vec[i] = vec[i] + HiccupFunc[i]()
    return vec
    
        