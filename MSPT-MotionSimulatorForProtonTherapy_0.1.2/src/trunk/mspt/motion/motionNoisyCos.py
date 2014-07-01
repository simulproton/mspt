########################################################################
#
# motionNoisyCos.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# October 2013
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

requiredKeys = ['breathingPeriod' , 'magnitude', 'initialPhase','variationsMag','variationsPeriod']
# requiredKeys = ['breathingPeriod' , 'magnitude', 'initialPhase','variationsMag','variationsPeriod', 'distributionMag','distributionPeriod']
# typeDistr = ['normal' , 'lognormal']


class MotionNoisyCos(object):
    '''Class modeling breathing motions as a cosine function. The model is based on models used by:
        
        * Lujan et al. 1999: "A method for incorporating organ motion due to breathing into 3D dose calculations"
        * George et al. 2005: "The application of the sinusoidal model to lung cancer patient respiratory motion"
        * Seco et al. 2009: "Breathing interplay effects during proton beam scanning: simulation and statistical analysis."
    
    The model of the breathing motion based on a cosine function was introduced by Lujan et al. and slightly modified later by George et al.\
    We use this last version that can be expressed as:
    
        
        f(t) = z0 + b cos( (2*Pi*t / tau)  - phi)
        
        
    where z0 is the average position (it is set to 0),t the time, b the amplitude, and tau the motion period.
        
    Lujan and Seco also added a phase phi that can be added in the cosine. We add a Gaussian noise on each coordinate to introduce uncertainties.
        
    We represent the noisy cosine function as: f(t) = z0 + b(t) cos( (2*Pi*t / tau(Cylcle))  + phi) where b(t) is a value of a normal distribution\
    whose mean is 'magnitude' and whose standard deviation is 'distributionMag'. tau(Cylcle) is the period which depends on the number of cylce, i.e.\
    every new cycle (Cycle = round(timer/(2*Pi)) ) the period tau is changed following the normal distribution. The mean is 'breathingPeriod' and the \
    standard deviation is 'distributionPeriod'.
        
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
        
        
        for key in requiredKeys:
            if key not in argsDict.keys():
                strErr = 'motionCos should be initialized with 3 arrays of parameters: %s , %s is missing.'%(requiredKeys,key)
                raise ValueError(strErr)
        self._typeFloat = typeFloat
        
        test=np.zeros((3),dtype=self._typeFloat)
        typetest= type(test)
        for idx,name in enumerate(requiredKeys):
            # if name in ['distributionMag','distributionPeriod']:
#                 continue
#             # item = eval(name)
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
        
        self._tabDirection = ['X','Y','Z']
        for idx,value in np.ndenumerate(argsDict['breathingPeriod']):
            if value < 0:
                strError = "In init motion Cos, %s direction has a breathing period < 0!"
                raise ValueError(strError)
        self._breathingPeriod = argsDict['breathingPeriod']
        self._magnitude = argsDict['magnitude']
        self._phase = argsDict['initialPhase']
        self._stdevMag = argsDict['variationsMag']
        self._stdevPeriod = argsDict['variationsPeriod']
#         for key, var, typ in zip( ['distributionMag','distributionPeriod'] ,['self._distrMag','self._distrPeriod'],['self._distrTypMag','self._distrTypPeriod'] ):
#             if  argsDict[key] in typeDistr:
#                 if argsDict[key] == 'normal':
#                     code = '%s = buildNormDistrFunction'%(var)
# #                     var = buildNormDistrFunction
#                     exec code
#                 elif argsDict[key] == 'lognormal':
#                     code = '%s = buildLogNormDistrFunction'%(var)
# #                     var = buildLogNormDistrFunction
#                     exec code
#                 else:
#                     strErr = "Distribution (for %s) %s has not been implemented yet... sorry... "%(key,argsDict[key])
#                     raise ValueError(strErr)
#                 code = '%s = argsDict[key]'%(typ)
# #                 typ = argsDict[key]
#                 exec code
#             else:
#                 strErr = "Unkown distribution type for %s. Should be in %s. But is %s"%(key,typeDistr,argsDict[key])
#                 raise ValueError(strErr)
        self._distrMag = buildNormDistrFunction
        self._distrPeriod = buildNormDistrFunction
        
        self._magFunc = []
        self._periodFunc = []
        for idx in range( len(self._tabDirection)):
            
            meanMag = self._magnitude[idx]
            stdMag = self._stdevMag[idx]
            funcMag = self._distrMag(meanMag,stdMag)
            self._magFunc.append(funcMag)

            
            meanPeriod = self._breathingPeriod[idx]
            stdPeriod = self._stdevPeriod[idx]
            funcPeriod = self._distrPeriod(meanPeriod,stdPeriod)
            self._periodFunc.append(funcPeriod)
        
        self._currCycle = None 
        self._currPeriodNoise = [None,None,None]
        self._currMagnitudeNoise = [None,None,None]

     
    def __str__(self):
        strValue = "Motion: p0 + ampl. x cos( 2*PI*t / tau  - phi ):\n\tdirections:[X, Y, Z]\n\t-Breathing period: %s sec\n\t-Amplitude: %s cm\n\t-Initial Phase: %s rad\n"%(str(self._breathingPeriod),str(self._magnitude),str(self._phase))
        strVariations = "\tMagnitude PDF type:%s , stdev = %s\n\tPeriod PDF type:%s , stdev = %s\n"%(self._distrTypMag,str(self._stdevMag),self._distrTypPeriod,str(self._stdevPeriod))
        return strValue+strVariations 
   
    
    def getDisplacementVectorAtTime(self,timer):
        '''Computes the displacement vector according to the equation f(t) = z0 + b(t) cos( (2*Pi*t / tau(Cylcle))  + phi) 
        
        :param  timer: Time in sec. 
        
        :returns: Displacement vector : numpy array with 3 elements. 
        '''
        if timer < 0 :
            raise ValueError("Time < 0 in get displacement vector")
        vec = np.zeros((3),dtype=self._typeFloat,order='C')
        currCycle = np.round(timer/(2*np.pi))
        if currCycle!=  self._currCycle:
            self._currCycle = currCycle
            for idx in range( len(self._tabDirection)):
                period = self._periodFunc[idx]()
                magnitude = self._magFunc[idx]()
                self._currPeriodNoise[idx] = period
                self._currMagnitudeNoise[idx] = magnitude
        for idx in range( len(self._tabDirection)):
            period = self._currPeriodNoise[idx]
            magnitude = self._currMagnitudeNoise[idx]
            if period == 0:
                vec[idx] = magnitude
            else:
                vec[idx] = magnitude* np.cos( ((timer * 2 * np.pi ) / period) + self._phase[idx] )
        return vec
    

# def buildLogNormDistrFunction(mean , stdev):
#     def func():
#         if stdev == 0 :
#             return mean
# #         print "m:%f"%mean
#         if mean == 0:
# #             print "m***:%f"%mean
#             newMean = 1e-6
#         else:
#             newMean = mean
#         variance = stdev * stdev
#         mu = calcMuForLogNorm(newMean , variance )
#         sigma = calcSigmaForLogNorm( newMean , variance )
#         return np.random.lognormal(mean = mu, sigma = sigma,size=1)[0]
#     return func
    
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

    
    
# def calcMuForLogNorm( mean , variance ):
#     return np.log( (mean*mean)/np.sqrt(variance+mean*mean))
#     
# def calcSigmaForLogNorm( mean , variance ):
#     return np.sqrt(np.log(1 + (variance/(mean*mean))))
#         
# def calcMeanForLogNorm( mu, sigma):
#     return np.exp(mu + (sigma * sigma)/2)
#     
# def calcVarianceForLogNorm( mu, sigm):
#     return np.exp(2*mu+sigma*sigma)*(np.exp(sigma*sigma)-1)
    