########################################################################
#
# dicomManager.py
# Dicom Reader Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On Nov, 8 2012
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

from dicomCTReader import DicomCTReader
from dicomRDReader import DicomRDReader
from dicomRPReader import DicomRPReader
from dicomRSReader import DicomRSReader
from dicomGeneralReader import DicomGeneralReader
import sys,os, glob
import dicom

class DicomManager(object):
    '''Dicom manager manage many DicomReaders which link to different types of dicom files. The dicom manager initialization creates a dicom reader for each path provided.
    
        More information on dicom standards can be found at Dicom `NEMA <http://medical.nema.org/standard.html>`_
    
    :param paths: A list of paths to the dicom files to open.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    
    '''
    def __init__(self, paths,typeFloat):
        ''' Initialization
        '''
        self._dicomReaders = []
        self._typeFloat = typeFloat
        for path in paths:
            if os.path.isdir(path):
                for infile in glob.glob( os.path.join(path, '*.dcm') ):
                    dicomFile = dicom.read_file(infile)
                    dcmFileType = str(dicomFile.SOPClassUID)
                    if dcmFileType == 'CT Image Storage':
                        reader = DicomCTReader(path,typeFloat)
                        break
                    else:
                        raise ValueError("Directory: Case not defined yet in dicom reader")
            else:
                dicomFile = dicom.read_file(path)
                dcmFileType = str(dicomFile.SOPClassUID)
                if dcmFileType == 'RT Dose Storage':
                    reader = DicomRDReader(path,typeFloat)
                elif dcmFileType == 'RT Ion Plan Storage':
                    reader = DicomRPReader(path,typeFloat)
                elif dcmFileType == 'RT Structure Set Storage':
                    reader = DicomRSReader(path,typeFloat)
                else:
                    reader = DicomGeneralReader(path,typeFloat)
                    

            self._dicomReaders.append(reader)
        self._fullPlan = self.checkFullRTPlan()
        print "Full RT Plan (CT, RP, RD, RS): %s"%str(self._fullPlan)
        self._sameFrameOfRefUID = self.checkFrameOfReferenceUID()
        print "Same Frame Of Reference UID: %s"%str(self._sameFrameOfRefUID)
        self._sameStudyInstanceUID = self.checkSameStudyInstanceUID()
        print "Same Study Instance UID: %s"%str(self._sameStudyInstanceUID)
        
        
    def checkFullRTPlan(self):
        '''Check if the dicom files of the dicom manager correspond to an entire RT plan: 1 CT image, 1 RS file, 1 RP file, 1 RD file.
        
        :returns: True if the dicom files correspond to an entire RT Plan, False otherwise
        
        '''
        self._idxCT = -1
        self._idxRS = -1
        self._idxRP = -1
        self._idxRD0 = -1

        for idx,reader in enumerate(self._dicomReaders):
            if reader.dcmType == 'CT Image Storage':
                self._idxCT = idx
            elif reader.dcmType == 'RT Structure Set Storage':
                self._idxRS = idx
            elif reader.dcmType == 'RT Ion Plan Storage':
                self._idxRP = idx
            elif reader.dcmType == 'RT Dose Storage':
                self._idxRD0 = idx
            else:
                print "%s : not treated yet"%reader.dcmType
            
            
        if self._idxCT != -1 and self._idxRS != -1 and self._idxRD0 !=-1 and self._idxRP != -1:
            return True
        return False
    
    def checkFrameOfReferenceUID(self):
        '''Check if all dicom files of the dicom manager have the same Frame Of Reference UID. 
        
        :returns: True if the dicom files all have the same Frame Of Reference UID, False otherwise
        
        '''
        flag = 0
        listFrameOfRef = []
        
        for idx,reader in enumerate(self._dicomReaders):
            listFrameOfRef.append(reader.dataForAttribute('FrameOfReferenceUID'))
        frameOfRef = listFrameOfRef[0]
        
        for idx,item in enumerate(listFrameOfRef[1:]):
            frameRef = listFrameOfRef[1+ idx].replace("'",'')
            if frameRef != frameOfRef:
                flag = 1
            
        if flag == 1:
            print "Warning: all the files don't have the same Frame Of Reference UID:"
            print listFrameOfRef
            print"--------"
            return False
        return True
        
    def checkSameStudyInstanceUID(self):
        '''Check if all dicom files of the dicom manager have the same Study Instance UID. 
        
        :returns: True if the dicom files all have the same Study Instance UID, False otherwise
        
        '''  
        flag = 0
        listStudyInstanceUID = []
        
        for idx,reader in enumerate(self._dicomReaders):
            listStudyInstanceUID.append(reader.dataForAttribute( 'StudyInstanceUID'))
        studyInstanceUID =  listStudyInstanceUID[0]

        for idx,item in enumerate(listStudyInstanceUID[1:]):
            studyInst = listStudyInstanceUID[1+ idx].replace("'",'')
            if studyInst != studyInstanceUID:
                flag = 1
     
        if flag == 1:
            print "Warning: all the files don't have the same Study Instance UID:"
            print listStudyInstanceUID
            print"--------"
            return False
        return True

    
    def fullPlan(self):
        '''Get the full RT plan dicom readers.
        
        :returns: None if the dicom manager do not have a full RT plan (CT, RS, RP, RD dicom files) or if only RD is missing.\
        Otherwise returns a dictionary with keys : 'CT Image Storage', 'RT Structure Set Storage', 'RT Ion Plan Storage', 'RT Dose Storage'. \
        The values are the dicom readers.\
        If RD is missing the value for 'RT Dose Storage' is None.\
        '''
        if self._fullPlan:
            return {'CT Image Storage':self._dicomReaders[self._idxCT],\
            'RT Structure Set Storage':self._dicomReaders[self._idxRS],\
            'RT Ion Plan Storage':self._dicomReaders[self._idxRP],\
            'RT Dose Storage':self._dicomReaders[self._idxRD0]}

        elif self._idxRD0 == -1 and  self._idxCT != -1 and self._idxRS != -1 and self._idxRP != -1:
            return {'CT Image Storage':self._dicomReaders[self._idxCT],\
            'RT Structure Set Storage':self._dicomReaders[self._idxRS],\
            'RT Ion Plan Storage':self._dicomReaders[self._idxRP],\
            'RT Dose Storage':None}        
        
        return None
    
    def getDicomReadersForPatient(self):
        '''Get the dicom readers needed to represent the patient, i.e. CT,RS,RD dicom files.
        
        :returns: None if the dicom manager do not have a full RT plan (CT, RS, RP, RD dicom files) or if only RD is missing.\
        Otherwise returns a dictionary with keys : 'CT Image Storage', 'RT Structure Set Storage', 'RT Dose Storage'. The values are the dicom readers.\
        If RD is missing the value for 'RT Dose Storage' is None
 
        '''
        if self._fullPlan:
            return {'CT Image Storage':self._dicomReaders[self._idxCT],\
            'RT Structure Set Storage':self._dicomReaders[self._idxRS],\
            'RT Dose Storage':self._dicomReaders[self._idxRD0]}
        
        elif self._idxRD0 == -1 and  self._idxCT != -1 and self._idxRS != -1 and self._idxRP != -1:
            return {'CT Image Storage':self._dicomReaders[self._idxCT],\
            'RT Structure Set Storage':self._dicomReaders[self._idxRS],\
            'RT Dose Storage':None}

        return None
    
    def getDicomReaderForScanningPath(self):
        '''Get the dicom reader needed to build the scanning path, i.e. RP dicom file.
        
        :returns: None if the dicom manager do not have a full RT plan (CT, RS, RP, RD dicom files)or if only RD is missing.\
        Otherwise returns dicom reader corresponding to the RP dicom file ('RT Ion Plan Storage').
                
        '''
        if self._fullPlan or (self._idxRD0 == -1 and  self._idxCT != -1 and self._idxRS != -1 and self._idxRP != -1):
            return self._dicomReaders[self._idxRP]
        return None
                    
    
    def getDicomReaders(self):
        '''Get the list of the dicom readers.
        
        :returns: A list of dicom readers.
        
        '''
        return self._dicomReaders
        
    