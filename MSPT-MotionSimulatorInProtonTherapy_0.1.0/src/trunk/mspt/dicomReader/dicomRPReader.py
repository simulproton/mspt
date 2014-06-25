########################################################################
#
# dicomRPReader.py
# Dicom Reader Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On Oct, 17 2012
# 
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

import dicom
import numpy
import sys, os
from scanningPath import ScanningPathMultipleBeams

import dicomReader
from dicomReader import DicomReader

class DicomRPReader(DicomReader):
    '''Dicom reader for dicom RT plans (RP dicom file). 

    :param path: Path to the RP dicom file.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    
    '''

    def __init__(self, path = None, typeFloat = 'float32'):
        self._RTPlan = None
        if path is not None:
            self.loadRTPlanFile(path)
            self._RTPlanPath = path
        self._typeFloat = typeFloat
    
    
    def loadRTPlanFile(self, path):
        '''Load a RT Plan from a Dicom file located at "path".
        
        :param path: Path where is located the dicom RP file.
        
        '''
        fileName, fileExtension = os.path.splitext(path)
        if fileExtension == ".dcm":
            self._RTPlan = dicom.read_file(path)
            if str( self._RTPlan.SOPClassUID) != 'RT Ion Plan Storage':
                strErr = "RP file is not 'RT Ion Plan Storage' , but %s"%self._dosePlan.SOPClassUID
                raise ValueError(strErr)
            self._FractionGroupSeq = self._RTPlan.FractionGroupSequence
            self._ScanningPath = None
            try:
                self._NbFractions = self._RTPlan.FractionGroupSequence[0].NumberOfFractionsPlanned
            except:
                print "Warning: FractionGroupSequence[0].NumberOfFractionsPlanned not found in RP plan"
                self._NbFractions = 1
            self._dcmType = 'RT Ion Plan Storage'
            print "%s loaded"%self._dcmType

        else:
            print "File extension : %s"%fileName
            strErr = "Error : %s files are not supported by DicomRPReader"%fileExtension
            raise ValueError(strErr)
        

    @property
    def fractionGroupSequence(self):
        '''Get access to the Fraction Group Sequence of the RP dicom file.
        
        :returns: Fraction Group Sequence dicom data from pydicom.
        
        '''    
        return self._FractionGroupSeq
        
    def dataForAttribute(self , attr):
        '''Get the value of a given attribute  for the RT plan and
        returns the value. 
        
        :param  attribute: The attribute is a string defined in the dicom documentation (data loaded with pydicom) or specific strings:
            
            * **FractionGroupSequence**: Fraction Group Sequence dicom data from pydicom.
            * **NumberOfFractionsPlanned**: Number of fractions used in the plan.
             
        :returns: Attribute value.
        
        '''        
        if attr == 'FractionGroupSequence':
            return self.fractionGroupSequence
        if attr == 'NumberOfFractionsPlanned':
            return self._NbFractions
        else:
            return getattr(self._RTPlan, attr, None)
    
        
    def getScanningPath(self , dictSettings):
        '''Get the ScanningPathMultipleBeams instance created with the RP file.
        
        :param dictSettings: Dictionary of settings for the scanning path
        
        '''
    
        plan = ScanningPathMultipleBeams(self._RTPlan, dictSettings,self._typeFloat)
        return plan    
    
    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader: 'RT Ion Plan Storage'.
    
        '''
        return self._dcmType  
        
        
    @property
    def rtPlanPath(self):
        '''
         :returns: The RP dicom file loaded with pydicom
    
        '''
        return self._RTPlanPath
 