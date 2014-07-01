########################################################################
#
# dicomRDReader.py
# Dicom Reader Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On June, 5 2012
# 
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


import dicom 
import numpy as np
import sys, os

import dicom.config 
dicom.config.enforce_valid_values = False

import dicomReader
from dicomReader import DicomReader

class DicomRDReader(DicomReader):
    '''Dicom reader for dicom dose grids (RD dicom file). It represents a 3D dose distribution. 

    :param path: Path to the RD dicom file.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    
    '''
    
    def __init__(self , path =None, typeFloat = 'float32'):

        self._typeFloat = typeFloat
        self._dosePlan = None
        if path is not None:
            self.loadRTDosePlanFile(path)
            self._dcmType = 'RT Dose Storage'
            print "%s loaded"%self._dcmType

        else:
            raise ValueError ("Path shouldn't be None to initialize a RD dicom reader")        
        
    def __str__(self):
        '''

        :returns: A string containing information about the dicom dose grid loaded
    
        '''

        strToPrint = "[RT Dose:\n\t-Frame Of Reference UID:%s\n\t-Dose units:%s\n\t\
        -Num. Frames:%i\n\t-Num. Rows:%i\n\t-Num. Cols:%i ]"\
        %(self.dataForAttribute('FrameOfReferenceUID'),str(self.DoseUnits),self.PixelArray.shape[0],self.PixelArray.shape[1],self.PixelArray.shape[2])
        return strToPrint
        
    def loadRTDosePlanFile(self, path):
        '''Load a RT Plan Dose from a Dicom file located at "path". It initialize all the attributes of the class.
        
        :param path: Path where the dicom RD file is located.
        
        '''
                
        fileName, fileExtension = os.path.splitext(path)
        if fileExtension == ".dcm":
            self._dosePlan = dicom.read_file(path)
            if str(self._dosePlan.SOPClassUID) != 'RT Dose Storage':
                strErr = "RD file is not 'RT Dose Storage' , but %s"%self._dosePlan.SOPClassUID
                raise ValueError(strErr)
            self._pixel_array = None
            self._nCols = int(self.dataForAttribute('Columns'))
            self._nRows = int(self.dataForAttribute('Rows'))
            self._nFrames = len( self.dataForAttribute('GridFrameOffsetVector'))
            exec "self._imagePosition = [np.%s(x) for x in self.dataForAttribute('ImagePositionPatient')]"%self._typeFloat
            exec "self._frameOffsetVector = [np.%s(x) for x in self.dataForAttribute('GridFrameOffsetVector')]"%self._typeFloat
            exec "self._imageOrientation = [np.%s(x) for x in self.dataForAttribute('ImageOrientationPatient')]"%self._typeFloat
            self._imageOrientation.append(0.0)
            self._imageOrientation.append(0.0)
            self._imageOrientation.append(1.0)
            assert(len(self._imageOrientation) == 9)
            exec "self._spacing = [np.%s(x) for x in self.dataForAttribute('PixelSpacing')]"%self._typeFloat
            self._spacingInfo = {'rows':self._spacing[0] , 'cols':self._spacing[1] , 'frames':self._frameOffsetVector[1] }
    	    self._path = path
        else:
            print "File extension : %s"%fileName
            strErr = "Error : %s files are not supported by DicomRDReader"%fileExtension
            raise ValueError(strErr)

    def getPath(self):
        '''Get the path to the RD dicom file.
        
        :returns: Path where the dicom RD file is located.
        
        '''
        return self._path      
  
  
    
    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader: 'RT Dose Storage'.
    
        '''
        return self._dcmType
    
    @property
    def PixelArray(self):
        '''Get the 3D numpy array of the dose distribution.
        
        :returns: 3D numpy array.
        
        '''    
        if self._pixel_array is None:
            exec "scale = np.%s(self.dataForAttribute('DoseGridScaling'))"%self._typeFloat
            self._pixel_array = scale * self.dataForAttribute('pixel_array')
        return self._pixel_array
    
    @property
    def DoseUnits(self):
        '''Get the dose units of the dose distribution.
        
        :returns: The value for the attribute 'DoseUnits'
        
        '''  
        return self.dataForAttribute('DoseUnits')
        
        
    @property
    def SpacingInfo(self):
        '''Get the spacing information of the dose grid.
        
        :returns: A dictionary with keys: 'rows', 'cols','frames' and values the spacings defined in the units defined in the dicom file.
        
        '''  
        return self._spacingInfo
    
    
    @property
    def dicomData(self):
        '''Get direct access to the dicom data.
        
        :returns: The dicom data loaded using the package pydicom.
        
        '''  
        return self._dosePlan
    
    
    def dataForAttribute(self ,  attribute):
        '''Get the value of a given attribute ( attribute should be a string,e.g. 'DoseUnits') for the dose planned and
        returns the value. 
        
        :param  attribute: The attribute is a string defined in the dicom documentation (data loaded with pydicom) or specific strings:
            
            * **PixelArray**: returns the 3D numpy array of the dose distribution.
            * **DoseUnits**: returns the dose units used in the dose distributions.
            * **SpacingInfo**: returns the voxel spacing as a dictionary (spacing between 'rows', 'cols' and 'frames').
            * **EndImagePositionPatient**: Coordinates of the center of the last voxel of the volume. For example if the volume shape is (m,n,p), it returns the coordinates if the voxel (m-1,n-1,p-1) in the units defined in the spacing information.
            * **ImageOrientation**: Orientation of the X,Y,Z axis in the dicom image expressed according to the columns,rows,frames. \
            For example: [1,0,0,0,1,0,0,0,1] for X from left to right, Y from top to bottom and Z from front to back. This represents the unit vectors of the 3D volume in the numpy array.
            
        :returns: Attribute value.
        
        '''        
        if  attribute == 'PixelArray':
            return self.PixelArray
        elif  attribute == 'DoseUnits':
            return self.DoseUnits
        elif  attribute ==  'SpacingInfo':
            return self.SpacingInfo
        
        elif  attribute == 'EndImagePositionPatient':
            exec  "imagePosition = [np.%s(x) for x in self.dataForAttribute('ImagePositionPatient')]"%self._typeFloat
            exec "frameOffsetVector = [np.%s(x) for x in self.dataForAttribute('GridFrameOffsetVector')]"%self._typeFloat
            if frameOffsetVector[0] != 0:
                raise ValueError('Error - GridFrameOffsetVector type not implemented')
            exec "spacing = [np.%s(x) for x in self.dataForAttribute('PixelSpacing')]"%self._typeFloat
            rows = int(self.dataForAttribute('Rows'))
            cols = int(self.dataForAttribute('Columns'))
            xEnd = (imagePosition[0] + (cols-1) * spacing[1])
            yEnd = (imagePosition[1] + (rows-1) * spacing[0])
            zEnd = (imagePosition[2] + frameOffsetVector[-1])
            return (xEnd,yEnd,zEnd)
        elif  attribute == 'ImageOrientation':
            assert(len(self._imageOrientation) == 9)
            return self._imageOrientation
        else:
            return getattr(self._dosePlan, attribute, None) 
            
