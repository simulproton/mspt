########################################################################
#
# motionMonitor.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# November 2013
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


class MotionMonitor(object):
    ''' The class MotionMonitor aims to simulate a motion monitoring system. A motion (motionManager) is assigned to this monitor. This motion \
    provides the true measurement to the monitor. Some noise is added to this measurement in order to make it more realistic. The noise \
    is modeled by a normal distribution. 
    
    :param motion: The motion being used. It should be a motion manager so that it can be used for any type of motion. 
    :param measureStdErr: Numpy array with 3 values corresponding to the standard deviation of the normal distribution for each coordinate.
    :param measureAvgErr: Numpy array with 3 values corresponding to the mean of the normal distribution for each coordinate. 
    
    .. note: 
    
        The units for the stdandard deviation and the average are in cm and the time in sec.
    
    '''


    def __init__(self, motion, measureStdErr, measureAvgErr):
        '''
        Initialistion of the motion monitoring. 
        '''

        self._motion = motion
        if len(measureStdErr) != 3:
            raise ValueError ("In MotionMonitor, measureStdErr should have a length of 3")
        self._funcErrors = []
        self._measureStdErr = measureStdErr
        self._measureAvgErr = measureAvgErr
        for (stdev,avg) in zip(measureStdErr,measureAvgErr):   
            self._funcErrors.append( buildNormDistrFunction(avg , stdev))
        
        self._history = {'time':[],'x': {'PatientMotion':[] , 'MeasuredMotion':[]}\
                ,'y': {'PatientMotion':[] , 'MeasuredMotion':[]}\
                ,'z': {'PatientMotion':[] , 'MeasuredMotion':[]}}
        
    def __str__(self):
        '''
        String information about the monitoring system
        '''
        strMonitor = "Monitoring system:\nMeasure noise simulated by normal distribution centered on zero."\
        +"\t-StDev along X IEC Fixed: %f cm\n"%(self._measureStdErr[0])\
        +"\t-StDev along Y IEC Fixed: %f cm\n"%(self._measureStdErr[1])\
        +"\t-StDev along Z IEC Fixed: %f cm\n"%(self._measureStdErr[2])\
        +"\t-Avg along X IEC Fixed: %f cm\n"%(self._measureAvgErr[0])\
        +"\t-Avg along Y IEC Fixed: %f cm\n"%(self._measureAvgErr[1])\
        +"\t-Avg along Z IEC Fixed: %f cm\n"%(self._measureAvgErr[2])
        return strMonitor
    
    def monitorPatientMotionAtTime(self, timer):
        '''Function that simulates the acquisition of the measurement. 
        
        :param timer: time in sec. 
        
        :returns: Dictionary: keys:'PatientMotion' ( patient True motion) ,  'MeasuredMotion' (motion measurement)
        
        '''
        
        patientMotion = self._motion.getDisplacementVectorAtTime(timer)
        measuredMotion = self.getMotionMeasurement(patientMotion)
        if timer not in self._history['time']:
            self._history['time'].append(timer)
            self._history['x']['PatientMotion'].append(patientMotion[0])
            self._history['y']['PatientMotion'].append(patientMotion[1])
            self._history['z']['PatientMotion'].append(patientMotion[2])
            
            self._history['x']['MeasuredMotion'].append(measuredMotion[0])
            self._history['y']['MeasuredMotion'].append(measuredMotion[1])
            self._history['z']['MeasuredMotion'].append(measuredMotion[2])

            
        return {'PatientMotion':patientMotion , 'MeasuredMotion':measuredMotion}
    
    def getHistory(self):
        '''
        :returns: Dictionary: keys:
            
            * 'time' : list of time values.
            
            * 'x' (dictionary): keys:
                
                * 'PatientMotion': list of true patient position along the X axis
                * 'MeasuredMotion': list of patient position measurement along the X axis

            * 'y' (dictionary): keys:
                
                * 'PatientMotion': list of true patient position along the Y axis
                * 'MeasuredMotion': list of patient position measurement along the Y axis

            * 'z' (dictionary): keys:
                
                * 'PatientMotion': list of true patient position along the Z axis
                * 'MeasuredMotion': list of patient position measurement along the Z axis

        '''
        return self._history
    

    def getMotionMeasurement(self,patientMotion):
        '''Computes the measurement vector. It adds Gaussian noise to true patient motion.
        
        :param patientMotion: Numpy array with 3 values corresponding to the patient true motion. 
        
        :returns: Numpy array with 3 values corresponding to the motion measurement.
            
        '''
        errVec = np.zeros(3)
        for i in range(3):
            errVec[i] = self._funcErrors[i]()
        return patientMotion+errVec
        

def buildNormDistrFunction(mean , stdev):
    '''Build a function that computes the values of a normal distribution.
    
    :param mean: Mean value of the distribution
    :param stdev: Standard deviation of the distribution.
    
    :returns: A callable function.
    
    '''
    def func():
        if stdev == 0:
            return mean
        return np.random.normal(loc = mean, scale = stdev,size=1)[0]
    return func