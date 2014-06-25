########################################################################
#
# File: tools.py
#  
# Proton Therapy Simulator Project
# Written by by Paul Morel, 
#               LIGM, Universite Paris-Est Marne La Vallee (France)
#               paul.morel@univ-mlv.fr
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


########################################################################
#
# Description:
#  
# 
# 
#               
#               
# 
#
########################################################################
import numpy as np


# Import C extensions (wrapped in C++ code) using float values (float32)
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSet import _fillDoseArray_PythonExt as _fillDoseArray
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSet import _fillDensityArray_PythonExt as _fillDensityArray
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSet import _fillCTArray_PythonExt as _fillCTArray
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSet import _fillStopPwrArray_PythonExt as _fillStopPwrArray
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeff import _calcRadiolDepth_PythonExt as _calcRadDeff 
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeff import _singleRay_PythonExt as _singleRay
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeff import _fewRays_RadDepth_PythonExt as _fewRaysRadDepth


# Import C extensions (wrapped in C++ code) using double values (float64)
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeffDouble import _calcRadiolDepth_PythonExtDouble as _calcRadDeffDouble 
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeffDouble import _singleRay_PythonExtDouble as _singleRayDouble
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionCalcDeffDouble import _fewRays_RadDepth_PythonExtDouble as _fewRaysRadDepthDouble
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSetDouble import _fillDoseArray_PythonExtDouble as _fillDoseArrayDouble
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSetDouble import _fillDensityArray_PythonExtDouble as _fillDensityArrayDouble
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSetDouble import _fillCTArray_PythonExtDouble as _fillCTArrayDouble
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBuildSpotSetDouble import _fillStopPwrArray_PythonExtDouble as _fillStopPwrArrayDouble

            
# def wrapFunction(func, *args):
#     '''Wrapper function:
#     
#     '''
#     def wrapper( val ):
#         ret = func(val,*args)
#         return ret
#     return wrapper
    
    
def fillPatientDose(dimNewGrid,originalDoseGrid,imagePosition,imageEndPosition,imageOrientation,spacingInfoCT,spacingInfoRD,imagePositionRD,typeFloat):
    '''Calls the C extension: fillDoseArray_PythonExt and returns its result.
        
    :param dimNewGrid: dimensions of the 3D array to fill with the dose contained in originalDoseGrid. Its size is m x n x q.
    :param originalDoseGrid: 3D array containing the planned dose extracted from RD dicom file. Its size can be different than m x n x q, \
    it can be smaller or greater.  
    :param imagePosition: coordinates (in mm) of the top left corner on the first slice of the volume representing the field of view. Equal to \
    imagePositionRD if the field of view of dose grid is smaller than the field of view of the CT image.
    :param imageEndPosition: coordinates (in mm) of the bottom right corner on the last slice of the field of view
    :param imageOrientation: axis information in orginialDoseGrid: 1D array with 9 elements
    :param spacingInfoCT: spacing information (mm) in newDoseGrid: x:spacingInfoCT[0] ,y:spacingInfoCT[1] ,z:spacingInfoCT[2] ,
    :param spacingInfoRD: spacing information (mm) in originalDoseGrid (from RD File) : x:spacingInfoRD[0] ,y:spacingInfoRD[1] ,z:spacingInfoRD[2] 
    :param imagePositionRD: coordinates (in mm) of the top left corner on the first slice of the original dose grid
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The dose grid filled with the values of the original dose grid and resampled to match the spacing of the CT image.

    '''

    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    
    if type(originalDoseGrid) != typetest:
        raise Exception('In fillSposetDose, originalDoseGrid is not *NumPy* array')
    if len(np.shape(originalDoseGrid)) != 3:
        raise Exception('In fillSposetDose, originalDoseGrid is not NumPy *3D array*')
    if originalDoseGrid.dtype != datatest:
        raise Exception('In fillSposetDose, originalDoseGrid is not *Float* NumPy array')
    
    
    for idx,item in enumerate([imagePosition,imageEndPosition,spacingInfoCT,spacingInfoRD,dimNewGrid,imagePositionRD]):
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In fillSposetDose, itemVec('+str(idx) +') is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In fillSposetDose, itemVec('+str(idx) +') is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In fillSposetDose, itemVec('+str(idx) +') is not *Float* NumPy array'
            raise Exception(stringError)
        
        test=np.zeros((9),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        if type(imageOrientation) != typetest:
            raise Exception('In fillSposetDose, imageOrientation is not *NumPy* array')
        if len(np.shape(imageOrientation)) != 1:
            raise Exception('In fillSposetDose, imageOrientation is not NumPy *1D array*')
        if imageOrientation.dtype != datatest:
            raise Exception('In fillSposetDose, imageOrientation is not *Float* NumPy array')

    for idx,name in enumerate(['dimNewGrid','originalDoseGrid','imagePosition','imageEndPosition','imageOrientation','spacingInfoCT','spacingInfoRD','imagePositionRD']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In fillSposetDose, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)

    if typeFloat == 'float32':
        return _fillDoseArray.fillDoseGrid(dimNewGrid,originalDoseGrid,imagePosition,imageEndPosition,imageOrientation,spacingInfoCT,spacingInfoRD,imagePositionRD)
    else:
        return _fillDoseArrayDouble.fillDoseGrid(dimNewGrid,originalDoseGrid,imagePosition,imageEndPosition,imageOrientation,spacingInfoCT,spacingInfoRD,imagePositionRD)

            
def fillPatientDensity(ctGrid,conversionTable,typeFloat):
    '''Calls the C extension: fillDensityArray_PythonExt and returns its result.
    
    :param ctGrid: 3D numpy array of CT values
    :param conversionTable: Table used for the conversion lits of tuples : [ ( ,) ,( ,) , ( ,) ,( ,)], \
    each tuple is composed of 2 elements, the CT value and the corresopnding mass density.
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
    
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The density grid filled with the density values from the original CT image.
       
    '''
        
        
    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    
    if type(ctGrid) != typetest:
        raise Exception('In fillSpotsetDensity, ctGrid is not *NumPy* array')
    if len(np.shape(ctGrid)) != 3:
        raise Exception('In fillSpotsetDensity, ctGrid is not NumPy *3D array*')
    if ctGrid.dtype != datatest:
        raise Exception('In fillSpotsetDensity, ctGrid is not *Float* NumPy array')
    
    
    if len(conversionTable) == 0:
        raise Exception('In fillSpotsetDensity, conversion table is not of length >0')
    for tup in conversionTable:
        if len(tup)!= 2:
            raise Exception('In fillSpotsetDensity, conversion table has tuple of size != 2')
        
    convTable = np.zeros((2,len(conversionTable)),dtype=typeFloat)
    for row in range(2):
        for col in range(len(conversionTable)):
            convTable[row,col] = conversionTable[col][row]
            
    for idx,name in enumerate(['ctGrid','convTable']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In fillSpotsetDensity, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)
        
    if typeFloat == 'float32':
        return _fillDensityArray.fillDensityGrid(ctGrid,convTable)
    else:
        return _fillDensityArrayDouble.fillDensityGrid(ctGrid,convTable)
        
def fillPatientCT(ctGrid,newShape,indStart,typeFloat):
    '''Calls the C extension: fillCTArray_PythonExt and returns its result.


    :param ctGrid: 3D numpy array of CT values
    :param newShape: shape of the new CT grid
    :param indStart: Indices where to start taking the CT values from the original CT grid in order to fill the new grid. 
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The reshaped CT grid.
        
    '''
        
        
    # .... Check arguments
    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype

    
    if type(ctGrid) != typetest:
        raise Exception('In fillSpotsetCT, ctGrid is not *NumPy* array')
    if len(np.shape(ctGrid)) != 3:
        raise Exception('In fillSpotsetCT, ctGrid is not NumPy *3D array*')
    if ctGrid.dtype != datatest:
        print "type ctGrid: %s"%ctGrid.dtype
        raise Exception('In fillSpotsetCT, ctGrid is not *Float* NumPy array')
    
    
    for idx,item in enumerate([newShape,indStart ]):
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In fillSposetCT, itemVec('+str(idx) +') is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In fillSposetCT, itemVec('+str(idx) +') is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In fillSposetCT, itemVec('+str(idx) +') is not *Float* NumPy array'
            raise Exception(stringError)
    
        for idx,name in enumerate(['ctGrid','newShape','indStart']):
            item = eval(name)
            if not item.flags['C_CONTIGUOUS']:
                strErr = 'In fillSposetCT, '+str(name) +' is not *C-Contiguous* NumPy array'
                raise ValueError(strErr)
        
    if typeFloat == 'float32':
        return _fillCTArray.fillCTGrid(ctGrid,newShape,indStart)
    else:
        return _fillCTArrayDouble.fillCTGrid(ctGrid,newShape,indStart)
        
            
def fillPatientStopPwr(densGrid,conversionTable,typeFloat):
    '''Calls the C extension: fillStopPwrArray_PythonExt and returns its result.

    :param densGrid: 3D grid of densiy values
    :param conversionTable: table used for the conversion from density to relative stopping power. Array of 3 rows and 4 columns. Row 0 : lowest bound of the density values for each medium, \
    Row 1 : upper bound of the density values for each medium , Row 3: corresponding Stopping Power value. Each column correspond to 1 type of medium.
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The relative stopping power grid.
         
     '''
    # .... Check arguments
    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype

    
    if type(densGrid) != typetest:
        raise Exception('In fillSpotsetStopPwr, densGrid is not *NumPy* array')
    if len(np.shape(densGrid)) != 3:
        raise Exception('In fillSpotsetStopPwr, densGrid is not NumPy *3D array*')
    if densGrid.dtype != datatest:
        raise Exception('In fillSpotsetStopPwr, densGrid is not *Float* NumPy array')
    
    test=np.zeros((3,conversionTable.shape[1]),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype    
    if type(conversionTable) != typetest:
        raise Exception('In fillSpotsetStopPwr, conversionTable is not *NumPy* array')
    if conversionTable.shape  != test.shape :
        strErr = 'In fillSpotsetStopPwr, conversionTable is not a *3 x %i* Numpy array'%conversionTable.shape[1]
        raise Exception(strErr)
    if conversionTable.dtype != datatest:
        raise Exception('In fillSpotsetStopPwr, conversionTable is not *Float* NumPy array')
    
    for idx,name in enumerate(['densGrid','conversionTable']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In fillSpotsetStopPwr, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)
        
    if typeFloat == 'float32':
        return _fillStopPwrArray.fillStopPwrGrid(densGrid,conversionTable)
    else:
        return _fillStopPwrArrayDouble.fillStopPwrGrid(densGrid,conversionTable)
    
def calcRadiolDepth(stPwrGrid,startVec,incVec,sourceVec,typeFloat):
    '''Calls the C extension: calcRadiolDepth_PythonExt and returns its result.
    
    :param stPwrGrid: 3D grid of relative stopping power values
    :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
    :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
    :param sourceVec: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The radiological depth grid.

    
    '''

    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    if type(stPwrGrid) != typetest:
        raise Exception('In calcRadiolDepth, densityGrid is not *NumPy* array')
    if len(np.shape(stPwrGrid)) != 3:
        raise Exception('In calcRadiolDepth, densityGrid is not NumPy *3D array*')
    if stPwrGrid.dtype != datatest:
        raise Exception('In calcRadiolDepth, densityGrid is not *Float* NumPy array')

    
    for idx,name in enumerate(['startVec','incVec','sourceVec']):
        item = eval(name)
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In calcRadiolDepth, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In calcRadiolDepth, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In calcRadiolDepth, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
        
    if len(np.where(stPwrGrid < 0)[0]) > 0:
        raise  ValueError("Stopping Power grid: values < 0 when going to compute Radio Depth")
    
    
    for idx,name in enumerate(['stPwrGrid','startVec','incVec','sourceVec']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In calcRadiolDepth, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)
                
        
    if typeFloat == 'float32':
        return _calcRadDeff.calcRadDeff(stPwrGrid,startVec,incVec,sourceVec)
    else:
        return _calcRadDeffDouble.calcRadDeff(stPwrGrid,startVec,incVec,sourceVec)


def singleRay(stPwrGrid,startVec,incVec,sourceVec, targetVec,typeFloat):
    '''Calls the C extension: _singleRay_PythonExt and returns its result.
    
    :param stPwrGrid: 3D grid of relative stopping power values
    :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
    :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
    :param sourceVec: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
    :param targetVec: (x,y,z) coordinates of the beam's target in cm in the in the IEC fixed coordinate system.
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The radiological depth grid.  
    
    '''
    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest= test.dtype
    if type(stPwrGrid) != typetest:
        raise Exception('In singleRay, densityGrid is not *NumPy* array')
    if len(np.shape(stPwrGrid)) != 3:
        raise Exception('In singleRay, densityGrid is not NumPy *3D array*')
    if stPwrGrid.dtype != datatest:
        raise Exception('In singleRay, densityGrid is not *Float* NumPy array')

    
    for idx,name in enumerate(['startVec','incVec','sourceVec','targetVec']):
        item = eval(name)
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In singleRay, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In singleRay, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In singleRay, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
        
    if len(np.where(stPwrGrid < 0)[0]) > 0:
        raise  ValueError("Stopping Power grid: values < 0 when going to compute Radio Depth")
    
    
    for idx,name in enumerate(['stPwrGrid','startVec','incVec','sourceVec','targetVec']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In singleRay, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)
    
    if typeFloat == 'float32':
        return _singleRay.singleRayTrace(stPwrGrid,startVec,incVec,sourceVec,targetVec)
    else:
        return _singleRayDouble.singleRayTrace(stPwrGrid,startVec,incVec,sourceVec,targetVec)


def calcRadiolDepthFewRays(stPwrGrid,startVec,incVec,sourceVec, targetVec,auxTargets,typeFloat):
    '''Calls the C extension:  _fewRays_RadDepth_PythonExt and returns its result.
    
    :param stPwrGrid: 3D grid of relative stopping power values
    :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
    :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
    :param sourceVec: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
    :param targetVec: (x,y,z) coordinates of the beam's target in cm in the in the IEC fixed coordinate system.
    :param auxVec: array of (x,y,z) coordinates of the beam's auxiliary targets (i.e., targets around the main target targetVec)\
    in cm in the in the IEC fixed coordinate system.
    :param typeFloat: 'float32' or 'float64'. It it the numpy float type to be used.
    
    .. note::
        
        All the arguments of the function (except typeFloat) must be numpy arrays of float values (typeFloat) and stored internally \
        according to the 'C' order.
        
    :returns: The radiological depth grid. 
    
    '''

    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    if type(stPwrGrid) != typetest:
        raise Exception('In calcRadiolDepthFewRays, densityGrid is not *NumPy* array')
    if len(np.shape(stPwrGrid)) != 3:
        raise Exception('In calcRadiolDepthFewRays, densityGrid is not NumPy *3D array*')
    if stPwrGrid.dtype != datatest:
        raise Exception('In calcRadiolDepthFewRays, densityGrid is not *Float* NumPy array')

    
    for idx,name in enumerate(['startVec','incVec','sourceVec','targetVec']):
        item = eval(name)
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
    
    for idx,name in enumerate(['auxTargets']):
        item = eval(name)
        test=np.zeros((1,3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if np.shape(item)[1] != np.shape(test)[1]:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In calcRadiolDepthFewRays, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
        
    if len(np.where(stPwrGrid < 0)[0]) > 0:
        raise  ValueError("Stopping Power grid: values < 0 when going to compute Radio Depth")
    
    
    for idx,name in enumerate(['stPwrGrid','startVec','incVec','sourceVec','targetVec','auxTargets']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In calcRadiolDepthFewRays, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)
                
        
    if typeFloat == 'float32':
        return _fewRaysRadDepth.fewRaysRadDepth(stPwrGrid,startVec,incVec,sourceVec, targetVec,auxTargets)
    else:
        return _fewRaysRadDepthDouble.fewRaysRadDepthDouble(stPwrGrid,startVec,incVec,sourceVec, targetVec,auxTargets)

