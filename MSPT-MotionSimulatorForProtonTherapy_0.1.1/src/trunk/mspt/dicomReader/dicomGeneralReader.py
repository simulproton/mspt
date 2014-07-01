########################################################################
#
# dicomGeneralReader.py
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
import sys, os, time
import tools

import dicomReader
from dicomReader import DicomReader


class DicomGeneralReader(DicomReader):
    '''A general dicom reader.
    
    :param path: Path to a directory where the CT slices are stored.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.

    '''
    
    def __init__(self , path = None, typeFloat = 'float32'):
        self._typeFloat = typeFloat
        self._dicomData = None
        if path is not None:
            fileName, fileExtension = os.path.splitext(path)
            if fileExtension == ".dcm":
                self._dicomData = dicom.read_file(path)
                try:
                    self._dcmType = str(self._dicomData.SOPClassUID)
                except: 
                    self._dcmType = "General Dicom Storage"
                
            print "%s loaded"%self._dcmType
    

    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader.
    
        '''
        return self._dcmType    
    

        
    
    def getFrameOfReferenceUID(self):
        '''
        Access the field providing the Frame of reference UID. 
        '''
        if  self._dicomData:
            try:
                frameOfRef =  self._dicomData.FrameOfReferenceUID
            except:
                frameOfRef = 'No_FrameOfReferenceUID'
        else:
            raise AttributeError("No Frame Of Reference UID because dicom data not loaded")
        return frameOfRef
        


        
    def dataForAttribute(self , attribute):
    
        if attribute == 'FrameOfReferenceUID':
            return self.getFrameOfReferenceUID()
        else:
            return getattr(self._dicomData, attribute, None)
            
