########################################################################
#
# dicomCTReader.py
# Dicom Reader Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On June, 5 2012
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
import numpy as np
import sys, os, glob
import dicomReader
from dicomReader import DicomReader
import tools as dcmTools


class DicomCTReader(DicomReader):
    '''Dicom reader for dicom series :CT images contained in the directory located at "path". 

    :param path: Path to a directory where the CT slices are stored.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'.
    
    '''
    def __init__(self, path, typeFloat = 'float32'):
        
        self._typeFloat = typeFloat
        self._dicomSeries = None
        if path is not None:
            self.loadDicomSeriesFromDirectory( path, sort = True, checkSeries = True)
            self._dcmType = 'CT Image Storage'
            print "%s loaded"%self._dcmType
        else:
            raise ValueError ("Path shouldn't be None to initialize a CT dicom reader")
        
    def __str__(self):
        '''
        :returns: A string containing information about the dicom image loaded
    
        '''

        if self._dicomSeries is None:
            return "No dicom series loaded"
        else:
            strDcmSeries="[dcmSeries:%s\n\t-Frame Of Reference UID:%s\n\t-Number Of Frames:%i\n\t-Number \
            Of Rows:%i\n\t-Number Of Columns:%i ]"\
            %(self._isSeries,self.getFrameOfReferenceUID(),self.numberOfFrames,self.numberOfRows,self.numberOfCols)
            return strDcmSeries
    
    def loadDicomSeriesFromDirectory(self, path, sort = True, checkSeries = True):
        '''Loads a Dicom series (A DICOM Series is a set of images with the same Series UID) from a given directory. It initializes all the attributes of the class.
                
        :param path: Path the directory where the CT slices are stored.
        :param sort: (True /False) sort the slices according to the slice location in the volume when loading the data.
        :param checkSeries: (True /False) check if all the slices belong to the same series

        '''
    
        DicomSeries = []
        self._frames = 0
        for infile in glob.iglob( os.path.join(path, '*.dcm') ):
            ds = dicom.read_file(infile)
            if str(ds.SOPClassUID) != 'CT Image Storage':
                continue
            DicomSeries.append(ds)
            self._frames = self._frames + 1
        
        self._isSorted = False
        self._isSeries = False
        self._3DPixelArrayForDisplay = None
        self._spacing = None
        self._rescaling = None
        self._3DPixelArrayRescaled = None
        self._dicomSeries = DicomSeries
        
        if sort:
            self.sortDcmSlicesInSeries()
        if checkSeries:
            listSeriesUID = self.getAttributeInSeries('SeriesInstanceUID')
            if len(set(listSeriesUID)) == 1:
                print "One dicom series has been loaded"
                self._isSeries = True
            else:
                print set(listSeriesUID)
                print "Warning: Provided dicom files are from different series"
        
      
    def getAttributeInSeries(self, attribute):
        '''Get a given attribute ( attribute should be a string,e.g. 'SliceLocation') for each slice in the dicom series.
        The attributes are the "tags" defined in the dicom documentation.
        
        :param attribute: String corresponding to the desired attribute.
        
        :returns: A list of values corresponding to the attribute for each slice in the series
        
        '''
        listAttr = []
        for idx,dcmSlice in enumerate(self._dicomSeries):
            val = getattr(dcmSlice, attribute, None)
            listAttr.append(val)
        return listAttr
        
    def getAttributeForSlice(self, attribute, sliceIdx = 0):
        '''Get a given attribute ( attribute should be a string,e.g. 'SliceLocation') for a slice in the dicom series.
        The attributes are the "tags" defined in the dicom documentation.
        
        :param attribute: String corresponding to the desired attribute.
        :param sliceIdx: slice index
        
        :returns: The value of the attribute.
        
        '''
        val = getattr(self._dicomSeries[sliceIdx], attribute, None)
        return val
        
    def printAttributeInSeries(self,attribute):
        '''Get a given attribute ( attribute should be a string,e.g. 'SliceLocation') for each slice in the dicom series and print it.
        The attributes are the "tags" defined in the dicom documentation.
        
        :param attribute: String corresponding to the desired attribute.
        
        '''
        listAttr = self.getAttributeInSeries(attribute)
        for idx,value in enumerate(listAttr):
            print "%s - Slice %i : %s" % (attribute ,idx , str(value))
        
    def sortDcmSlicesInSeries(self):
        '''Sort slices in a series by slice location.
        
        '''
        if self._dicomSeries:
            newList = sorted(self._dicomSeries , key = lambda ds: ds.SliceLocation)
            self._dicomSeries = newList
            self._isSorted = True


    def get3DPixelArray(self, rescale = True):
        '''Creates a 3D pixel array from the series if slices are sorted and from one same series
        
        :param rescale: True if one wants to use the 'rescale' 'intercept' attributes of the dicom files to obtain the intensity on Hounsfield Units. False otherwise.
        
        :returns: The 3D numpy array created.
        
        '''
        if self._3DPixelArrayForDisplay is None:
            if self._isSorted and self._isSeries:
                self._3DPixelArrayForDisplay = []
                nbFrames = len(self._dicomSeries)
                for idx,dcmSlice in enumerate(self._dicomSeries):
                    self._3DPixelArrayForDisplay.append(dcmSlice.pixel_array)
                print "CT loaded: %i slices"%len(self._dicomSeries)
        self._3DPixelArrayForDisplay = np.array(self._3DPixelArrayForDisplay,dtype=self._typeFloat)
        if rescale:
            if self._3DPixelArrayRescaled is None:
                if self._rescaling is None:
                    self.getRescalingPixelInfoSeries()  
                slope = self._rescaling['slope']
                intercept = self._rescaling['intercept']
                self._3DPixelArrayRescaled = slope * self._3DPixelArrayForDisplay + intercept
            return self._3DPixelArrayRescaled
        else:
            return self._3DPixelArrayForDisplay
            
        
        
    def getSpacingSeries(self):
        '''Get spacing information. 
        
        :returns: A dictionary with 'units' , spacing between 'rows', 'cols' and 'frames'
        
        '''
        print "Get spacing info"
        if not self._spacing:
            if self._isSorted and self._isSeries:
                self._spacing = dict()
                self._spacing['units'] = 'mm'
                spacing = self.getAttributeForSlice('PixelSpacing',0)
                exec "self._spacing['rows'] = np.%s(spacing[0])"%self._typeFloat #Vertical Spacing
                exec "self._spacing['cols'] = np.%s(spacing[1])"%self._typeFloat #Horizontal Spacing
                exec "self._spacing['frames'] = np.%s(self.getAttributeForSlice('SliceThickness',0))"%self._typeFloat #Depth Spacing

        if not self._spacing:            
            print "Pixel Spacing in Dicom Series couldn't be created because:"
            if not self._isSorted:
                print "\t - Slices are not sorted"
            if not self._isSeries:
                print "\t - Slices are from different series"
                    
        return self._spacing
    
    def getRescalingPixelInfoSeries(self):
        '''Get rescaling information. 
        
        :returns: a dictionary with 'units' (HU: Hounsfield Units, OD: optical Density and SU: unspecified) 'slope' for rescaling slope and 'intercept' for rescaling intercept. 
        
        '''
        print "Get rescaling info"
        if not self._rescaling:
            if self._isSorted and self._isSeries:
                self._rescaling = dict()
                rescaleType = self.getAttributeForSlice('RescaleType',0)
                if rescaleType:
                    self._rescaling['units'] = rescaleType
                else:
                    self._rescaling['units'] = 'HU'
                
                exec "self._rescaling['slope'] = np.%s(self.getAttributeForSlice('RescaleSlope',0))"%self._typeFloat
                exec "self._rescaling['intercept'] =np.%s(self.getAttributeForSlice('RescaleIntercept',0))"%self._typeFloat

                
                if not self._rescaling:            
                    print "Pixel Value Rescaling in Dicom Series couldn't be created because:"
                    if not self._isSorted:
                        print "\t - Slices are not sorted"
                    if not self._isSeries:
                        print "\t - Slices are from different series"
                        
            return self._rescaling       
        
    def getFrameOfReferenceUID(self):
        '''Get and return the Frame Of Reference UID for each slice or for the series
        
        :returns: The Frame of Reference UID if all the dicom slices loaded belong to the same dicom series. Otherwise returns a list  of all the frame of reference UIDs.
        
        '''
        if self._isSeries:
            return self.getAttributeForSlice('FrameOfReferenceUID',0)
        else:
            return self.getAttributeInSeries('FrameOfReferenceUID')
            
    
    @property
    def numberOfFrames(self):
        '''Get the number of frames (i.e. number of slices)
    
        :returns: The number of frames if the slices correspond to a series otherwise -1
    
        '''
        if self._isSeries:
            return int(self._frames)
        else:
            return -1

        
    @property
    def numberOfRows(self):
        '''Get the number of rows 
    
        :returns: The number of rows if the slices correspond to a series otherwise -1
    
        '''

        if self._isSeries:
            return int(self.getAttributeForSlice('Rows',0))
        else:
            return -1
    @property
    def numberOfCols(self):
        '''Get the number of columns 
    
        :returns: The number of columns if the slices correspond to a series otherwise -1
        
        '''

        if self._isSeries:
            return int(self.getAttributeForSlice('Columns',0))
        else:
            return -1    


        
    @property
    def dcmType(self):
        '''
         :returns: The type of dicom reader: 'CT Image Storage'.
    
        '''
        return self._dcmType
            
    def dataForAttribute(self , attr):
        '''Get the data for an attribute. 
    
        :param attr: String corresponding to the attibute. Attributes are either those defined in the dicom documentation of specific attributes:
        
            * **Shape**: shape of the volume represented in the CT series (-> returns: Dictionary with keys 'Frames', 'Rows', 'Cols')
            * **PixelArray**: get the 3D numpy array. (-> returns: 3D pixel array where the voxel values are scaled to be in the units defined in the dicom files.)
            * **FrameOfReferenceUID**: returns the frame of reference UID of the series.
            * **SpacingInfo**: return the voxel spacing as a dictionary ( 'units' , spacing between 'rows', 'cols' and 'frames' ).
            * **EndImagePositionPatient**: Coordinates of the center of the last voxel of the volume. For example if the volume shape is (m,n,p), it returns the coordinates if the voxel (m-1,n-1,p-1) in the units defined in the spacing information.
            * **ImageOrientationPatient**: Orientation of the X,Y,Z axis in the dicom image expressed according to the columns,rows,frames. \
            For example: [1,0,0,0,1,0,0,0,1] for X from left to right, Y from top to bottom and Z from front to back. This represents the unit vectors of the 3D volume in the numpy array.

        :returns: Value for the attribute if this is a dicom series, None otherwise.
    
        '''

        if self._isSeries:
            if attr == 'Shape':
                return { 'Frames':self.numberOfFrames , 'Rows':self.numberOfRows , 'Cols':self.numberOfCols}
            if attr == 'PixelArray':
                return self.get3DPixelArray( rescale = True)
            if attr == 'RescaleInfo':
                return self.getRescalingPixelInfoSeries()
            if attr == 'FrameOfReferenceUID':
                return self.getFrameOfReferenceUID()
            if attr == 'SpacingInfo':
                return self.getSpacingSeries()

            if attr == 'EndImagePositionPatient':
                if self._typeFloat == 'float64':
                    imagePosition = [np.float64(x) for x in self.dataForAttribute('ImagePositionPatient')]
                else:
                    imagePosition = [np.float32(x) for x in self.dataForAttribute('ImagePositionPatient')]
            
                spacing = self.dataForAttribute('SpacingInfo')
                rows = self.numberOfRows
                cols = self.numberOfCols
                frames = self.numberOfFrames
                xEnd = (imagePosition[0] + (cols-1) * spacing['cols'])
                yEnd = (imagePosition[1] + (rows-1) * spacing['rows'])
                zEnd = (imagePosition[2] + (frames-1) * spacing['frames'])
                return (xEnd,yEnd,zEnd)
            elif attr == 'ImageOrientationPatient':
                imageOrientation = self.getAttributeForSlice( attr, 0)
                frameDirection = dcmTools.crossProduct( imageOrientation[:3] , imageOrientation[3:])
                imageOrientation.append(frameDirection[0])
                imageOrientation.append(frameDirection[1])
                imageOrientation.append(frameDirection[2])
                return np.array(imageOrientation,dtype = self._typeFloat, order = 'C')
                
                
            else:
                return self.getAttributeForSlice( attr, 0)
        else:
            return None
        
        