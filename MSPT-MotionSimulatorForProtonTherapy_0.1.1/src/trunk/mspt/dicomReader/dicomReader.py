########################################################################
#
# dicomReader.py
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
import sys, os, glob


class DicomReader(object):
    '''Dicom reader: general class to read dicom files. 

    :param path: Path to the dicom file or dicom directory (if CT series).
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    
    '''
    def __init__(self, path,typeFloat = 'float32' ):
        '''
        
        '''
        self._dcmType = None # Types available: 'CT Image Storage' , 'RT Dose Storage',
                            # 'RT Ion Plan Storage',  'RT Structure Set Storage'

     
        
    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader: 'RT Dose Storage', 'RT Ion Plan Storage', 'RT Structure Set Storage' and 'CT Image Storage' or 'General Dicom Storage'.
    
        '''
        raise NotImplementedError('dcmType has not been implemented')
        
        
    def dataForAttribute(self , attribute):
        '''Get the data for an attribute. 
    
        :param attribute: String corresponding to the attibute. Attributes are either those defined in the dicom documentation of specific attributes defined by each specific class
        :returns: Value for the attribute if this is a dicom series, None otherwise.
    
        '''  
        raise NotImplementedError('dataForAttribute has not been implemented')

