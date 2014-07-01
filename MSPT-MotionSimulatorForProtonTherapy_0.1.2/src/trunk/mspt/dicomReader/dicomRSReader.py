########################################################################
#
# dicomRSReader.py
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


class DicomRSReader(DicomReader):
    '''Dicom reader for dicom Structure Sets (RS dicom file). 

    :param path: Path to the RS dicom file.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    

    '''
    
    def __init__(self , path = None, typeFloat = 'float32'):
        self._typeFloat = typeFloat
        self._RTStructureSet = None
        if path is not None:
            self.loadRTStrucSetFile(path)
        self._dcmType = 'RT Structure Set Storage'
        print "%s loaded"%self._dcmType
    
        
    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader: 'RT Structure Set Storage'.
    
        '''
        return self._dcmType    
    
    def loadRTStrucSetFile(self, path):
        '''Load a Structure Set from a Dicom file located at "path".
        
        :param path: Path where the dicom RS file is located.
        
        '''

        fileName, fileExtension = os.path.splitext(path)
        if fileExtension == ".dcm":
            self._RTStructureSet = dicom.read_file(path)
            if str( self._RTStructureSet.SOPClassUID) != 'RT Structure Set Storage':
                strErr = "RS file is not 'RT Structure Set Storage' , but %s"%self._dosePlan.SOPClassUID
                raise ValueError(strErr)
        else:
            print "File extension : %s"%fileName
            strErr = "Error : %s files are not supported by DicomRSReader"%fileExtension
            raise ValueError(strErr)
        
    
    def getFrameOfReferenceUID(self):
        '''Access the field providing the Frame of reference UID. 
            
        :returns: The Frame of reference UID defined for the first region of interest (ROI). 
        '''
        if self._RTStructureSet:
            return self._RTStructureSet.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
        else:
            raise AttributeError("No Frame Of Reference UID because RT Struct. Set file loaded")
    
    def getNumberOfROIs(self):
        ''' Access and return 
        
        :returns: The number of regions of interest (ROIs) defined in the dicom file. 
        '''
        if self._RTStructureSet:
            return len(self._RTStructureSet.ROIContourSequence)
        else:
            raise AttributeError("No ROIs because RT Struct. Set file loaded")
        
    
        
    def dataForAttribute(self , attribute):
        '''Get the value of a given attribute for the structure set and
        returns the value. 
        
        :param  attribute: The attribute is a string defined in the dicom documentation (data loaded with pydicom) or specific strings:
            
            * **FrameOfReferenceUID**: Frame Of Reference UID
             
        :returns: Attribute value.
        
        '''        
        if attribute == 'FrameOfReferenceUID':
            return self.getFrameOfReferenceUID()
        else:
            return getattr(self._RTStructureSet, attribute, None)
            

    def getMasksForAssociatedDicomImage(self,imagePosition, imageOrientation, pixelSpacing, imageShape):
        '''The goal of this function is to create a 3D mask for each ROI contained in the RS file.
        
        It is assumed the current RS file is assiciated to a dicom image whose information is received by this function.

        For each ROI :
        
            #. First create an empty mask
            #. Then go through each contour defined in ContourSequence.ContourData
            
            For each contour:
                
                #. Convert each of these contour coordinates into indices in the final volume.
                #. Send thes indices to the function tools.getMaskForContour()
                #. Receive a mask from tools.getMaskForContour() and set the slice (according to slice index) of the 3D mask to the value of the one received by the function. 

        
        :param imagePosition: Coordinates of the first voxel (top left corner of the first slice) of the dicom image
        :param imageOrientation: axis orientation of the dicom image
        :param pixelSpacing: spacing between frames, cols and rows in the image
        :param imageShape: 3D shape of the dicom image(frames, rows, columns)
        
        
        :returns: A list of tuple in which are stored as: ( mask 3D Numpy array , ROI name, ROI observation, ROI index, ROI type)  - Type = 'CLOSED_PLANAR' or 'POINT' 
        
        '''
    
        roisList = list()
        triplets = tools.tripletwise
        nbRois = len(self._RTStructureSet.ROIContourSequence)
        count = 0
        for roi in self._RTStructureSet.ROIContourSequence:
            count = count + 1
            num = roi.ReferencedROINumber
            name = self.getROINameForRefROINumber(num)
            print "\tROI %i / %i : %s"%(count,nbRois,name)
            obs = self.getROIObservationForRefROINumber(num)
            typeROI = roi.ContourSequence[0].ContourGeometricType
            if typeROI == 'CLOSED_PLANAR':
                mask = np.zeros(imageShape, dtype= "int8", order = 'C')
#                 t0 = time.clock() 
                for contour in roi.ContourSequence:
                    listPoints = triplets(contour.ContourData)
                    if listPoints == []:
                        raise ValueError("RSReader: listPoints empty")
                    listIndices = []
                    frame = -1
                    for point in listPoints:
                        (f,r,c) = tools.coordinatesToIndexFromImagePosition(point, imagePosition, imageOrientation, pixelSpacing)
                        listIndices.extend((c,r))
                        frame = f
                    if listIndices == []:
                        raise ValueError("RSReader:listIndices empty")
                    if len(listIndices)%2 != 0:
                        raise ValueError("RSReader: listIndices has not an even number of values")
                    maskSlice = tools.getMaskForContour(listIndices , imageShape[1:] , outline = 1 , fill = 1 , background = 0 )
                    if 0 < frame < imageShape[0]:
                        mask[frame,:,:] = mask[frame,:,:] | maskSlice[:,:]
                    else:
                        print "Mask: frame %i skipped"%frame

#                 t1 = time.clock()
                roisList.append( [mask,name,obs,num,typeROI] )
#                 print "Mask %s (%s : %i) done in %e sec" %(name,obs,num,(t1-t0))
            elif typeROI == 'POINT':
                roisList.append( [roi.ContourSequence[0].ContourData,name,obs,num,typeROI] )
            else:
                print "ROI type: %s"%str(typeROI)
                raise ValueError("Unknown roi type. Should be POINT or CLOSED_PLANAR")
        return roisList
    
    def getROINameForRefROINumber(self,index):
        '''
        
        :param index: Index of an ROI
        
        :returns: The ROI name from the given ROI index 
        '''
        for roi in self._RTStructureSet.StructureSetROISequence:
            if index == roi.ROINumber:
                return roi.ROIName
        print "Error - no ROI number %i"%index
        raise ValueError("No ROI name found for ROI number")


    def getROIObservationForRefROINumber(self,index):
        '''
        
        :param index: Index of an ROI
        
        :returns: The ROI observation from the given ROI index
        '''
        for roi in self._RTStructureSet.RTROIObservationsSequence:
            if index == roi.ReferencedROINumber:
                return roi.RTROIInterpretedType
        print "Error - no ROI number %i"%index
        raise ValueError("No ROI observation type found for ROI number")

    