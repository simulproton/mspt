########################################################################
#
# patient.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On 04/06/2012
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
##
########################################################################
import tools as patientTools
import math
import numpy as np
import dicom, dicom.UID
from dicom.dataset import Dataset, FileDataset
from ..dicomReader import ctToDensityConv  as ctToDensConv
from ..dicomReader import tools as dcmTools
from ..mathTools import tools as mathTools
import os, sys, time, gc
import cPickle
from ..dataViewer import array2DViewer as viewer2D
from operator import itemgetter, attrgetter
import itertools
from ..simulator.tools import Tools as simTools

templateRDFile = '/RefData/RD1_TemplateFile.dcm'


class Patient(object):

    '''
    The class "Patient" aims to represent the patient. The patient is represented by a set of 3D matrices. It is built from the CT dicom image and the RS (Structure set)\
     dicom file. The 3D matrices created are:
     
        * **CT Matrix**: 3D matrix representing the CT image in Hounsfield Units.
        * **Density Matrix**: 3D matrix representing the densities at each voxel of the CT Matrix. The conversion from CT number to density is \
        made by the dicomReader.ctToDensityConv data. 
        * **Relative Stopping Power matrix** : 3D matrix representing the relative (to water) stopping power at each voxel of the CT Matrix. It is built using the \
        data provided by physics.stoppingPower and the density matrix.
        * **Expected dose distribution** (optional): It corresponds to the 3D matrix of the dose distribution (RD dicom file) computed beforehand when planning the treatment.\
        If the RD dicom file is provided to the simulator the CT will be resized to cover the same volume as the RD dicom file. Moreover,\
        the matrix containing the RD data is resampled to have the same spacing as the CT matrix.
        * **Volumes of Interest (VOIs) 3D masks**: For each volume structure defined in the RS dicom file, a 3D matrix is created to store a binary mask \
        (1 in the volume and 0 outside) of each structure.
        * **Three 3D arrays** : X,Y,Z corresponding to the x,y,z coordinates in the dicom image of each voxel.
 
    :param stopPwrDataRef: Stopping power data ( physics.stoppingPower.StoppingPower) to be used for the conversion from density to relative stopping power.
    :param dcmData: Dictionary containing the dicom readers. The keys are:'CT Image Storage','RT Structure Set Storage','RT Dose Storage'  
    :param configData:  Dictionary configuring the motion. Mandatory keys are:
    
    
        * *'mainPath'* : corresponds to the path where any simulation will be saved.
        * *'nameSimulation'* : corresponds to the path where the current simulation will be saved.
        * *'typeFloat'* : The type of numpy float to be used: 'float32' or 'float64'
        * *'staticDelivery'* : True if the patient is not moving, or False otherwise
        * *'ctPadding'* (if 'staticDelivery' is False) : True if the user wants to add margins of a specified density (see 'paddedCTValue' ) around\
        the current CT image. This is useful, when the moving patient can be outside of the bounding box formed by the borders of the CT image. 
        * *'padValueFRC'* (if 'ctPadding' is True) : Numpy array of 3 integer values : [ number of added frames : nbF, number of added rows : nbR, \
        number of added columns : nbC]. For example : [ 0,0,10] will add 10 columns to the right of the CT image and 10 columns to the left of the CT image.
        * *'paddedCTValue'* (if 'ctPadding' is True): The value in CT number to use to pad the CT image. For example: air = -1000, water = 0.
        * *'storeMasksAsCPKL'* : True if the user wants to store the masks of the VOIs as numpy array in a cPickle structure (this allows to have access to the \
        3D numpy array outdise the simulation). False otherwise. Note: whether the value is set to True or False the masks are saved as png images.
        * *'skipSaveDensityImage'* : True to skip the storage of the 3D density matrix as png images. False otherwise. Moreover, if False and the patient is moving\
        it will stored the 3D density array every time the array is updated. 
        * *'listSaveDensImAlong'* (if 'skipSaveDensityImage' is False) : list of the axis along which the user wants to save the density images. Possible values are\
        'x' (i.e. the images are stored along the column axis of the dicom image), 'y' (i.e. the images are stored along the rows axis of the dicom image), \
        'z' (i.e. the images are stored along the frame axis of the dicom image). This can be used when the patient us moving and the user wants to see the motions\
        from the side, the top or the front of the dicom image. Examlpes: listSaveDensImAlong = [] , or listSaveDensImAlong = ['x'] or listSaveDensImAlong = ['x','z'] ... 
        * *'saveDensCPKL'* (if 'skipSaveDensityImage' is False) : True if the user wants to store the density 3D array in a cPickle structure (this allows to have access to the \
        3D numpy array outdise the simulation). False otherwise. Moreover, if False and the patient is moving\
        it will stored the 3D density array every time the array is updated. 
        * *'skipSaveStpPwrImage'* : True to skip the storage of the 3D relative stopping power matrix as png images. False otherwise, however, if the patient is moving it only stores it once.
        * *'listSaveStPwrImAlong'* (if 'skipSaveStpPwrImage' is False) :  Similar to 'listSaveDensImAlong', to save the images of the 3D relative stopping power matrix.
        * *'saveStpPwrCPKL'* (if 'skipSaveStpPwrImage' is False): Similar to 'saveDensCPKL' to save the 3D array of the relative stopping power. However, if the patient is moving it only stores it once.
        * *'skipSaveDeff'* : Similar to 'skipSaveStpPwrImage' to save the 3D array of the radiologcial depth. If the patient is moving it stored the array at every update.  
        * *'listSaveDeffAlong'* (if 'skipSaveDeff' is False) : Similar to 'listSaveDensImAlong', to save the images of the 3D radiological depth matrix.  
        * *'saveDeffCPKL'* (if 'skipSaveDeff' is False) :  Similar to 'saveDensCPKL' to save the 3D radiological depth matrix. If the patient is moving it stores the array at every update. 
        * *'MeV'* : Constant to convert 1MeV to Joules: 1MeV = 1.602e-13J
        * *'protPerMU'* : Number of protons per MU to use.

    
    
  
        
    
    
    '''

    
    def __init__(self,stopPwrDataRef, dcmData, configData ):
        '''
        Creates a patient from dicom data
        '''
        self._globalVar = configData
        self._pathToOutput = self._globalVar.mainPath + self._globalVar.nameSimulation
        if dcmData is not None and stopPwrDataRef is not None:
            self._spacingInfo = dict()
            self.initPatientWithDicomData(stopPwrDataRef, dcmData)
            
        else:
            raise ValueError("No input data for Patient")

            
        
    def initPatientWithDicomData(self,stopPwrDataRef, dcmData):
        ''' Create a patient from Dicom Data. Dicom Data should be provided by a dictionnary with 3 keys:
        'CT Image Storage','RT Structure Set Storage','RT Dose Storage'. Their values are respectively dicomCTReader, dicomRSReader,
        dicomRDReader.
        
        :param stopPwrDataRef: Stopping power data ( physics.stoppingPower.StoppingPower) to be used for the conversion from density to relative stopping power.
        :param dcmData: Dictionary containing the dicom readers. The keys are:'CT Image Storage','RT Structure Set Storage','RT Dose Storage'  
        
        '''
        if dcmData is None:
            raise ValueError("DcmData is None in initPatientWithDicomData.")
        self._dcmCTReader = dcmData['CT Image Storage']
        self._dcmRSReader = dcmData['RT Structure Set Storage']
        dcmRDReader = dcmData['RT Dose Storage']
        self._stopPwrDataRef = stopPwrDataRef
        if dcmRDReader is not None:
            self._pathRDFile = dcmRDReader.getPath()
            self._RDFileLoaded = True
        else:
            self._pathRDFile = os.getcwd() + templateRDFile
            self._RDFileLoaded = False
        self.buildPatient(self._dcmCTReader,dcmRDReader,self._dcmRSReader,self._stopPwrDataRef )
        


        
        
    def saveDoseDistrAsRDDcmFile(self, newDose , filename):
        '''Store the "newDose" in a RD Dicom file. 
        
        :param newDose: 3D numpy array storing a dose distribution in Gy.
        :param filename: Name of the file to be saved. It will be saved in at 'mainPath'/'nameSimulation'/  (see class initialization).
        
        .. note::
        
            If the user provides a RD dicom file to the simulator, the new RD files created will used the values of the tags used in the input RD file,\
            otherwise it will rely on the tags and values used in the CT images and the RS dicom file and use a template of RD dicom file stored in the package \
            of the simulator. 
        
        '''
        rdDcmFile = dicom.read_file(self._pathRDFile)
        
        file_meta = Dataset()
        for tag in rdDcmFile.file_meta.dir():
            code = 'file_meta.%s'%tag+' = rdDcmFile.file_meta.%s'%tag
            try:
                exec code
            except:
                print "\tInfo: %s was not in the RD file tag"%tag

        
        if not filename.endswith('.dcm'):
            filename = filename +  '.dcm'
            
        pathToSave = self._pathToOutput #+ "NewRefDose"
            
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)
        
        
        newRD = FileDataset(pathToSave+filename, {},file_meta = file_meta,preamble=rdDcmFile.preamble,\
                    is_implicit_VR = rdDcmFile.is_implicit_VR , is_little_endian = rdDcmFile.is_little_endian)
        
        skippedTags = ['DoseSummationType','Columns','Rows','NumberOfFrames','PixelSpacing','SliceThickness','GridFrameOffsetVector','DoseGridScaling','PixelData']

        for tag in rdDcmFile.dir():
            if tag not in skippedTags:
                code = 'newRD.%s'%tag+' = rdDcmFile.%s'%tag
                try:
                    exec code
                except:
                    print "\tInfo:%s was not in the RD file tag"%tag
 
#         print "rd array type:%s"%(str(rdDcmFile.pixel_array.dtype))
#         print "rd dose shape: %s"%(str(newDose.shape))

        rdDcmFile.DoseSummationType = 'PLAN'

        typeArray = rdDcmFile.pixel_array.dtype
        if typeArray == 'uint16':
            maxVal = 65535
        elif typeArray == 'uint8':
            maxVal = 255
        else:
            strErr = "newType used in dicom not defined yet: %s"%typeArray
            raise ValueError(strErr)

        newRD.Columns = np.int(newDose.shape[2])
        newRD.Rows = np.int(newDose.shape[1])
        newRD.NumberOfFrames = np.int(newDose.shape[0])
        newRD.PixelSpacing = [self.getSpacingInfo('rows'),self.getSpacingInfo('cols')]
        newRD.SliceThickness = self.getSpacingInfo('frames')
        newRD.GridFrameOffsetVector =  [int(0+newRD.SliceThickness * i) for i in range(newRD.NumberOfFrames)]

        exec "gridScaling = np.%s(np.max(newDose))  / np.%s(maxVal)"%(self._globalVar.typeFloat,self._globalVar.typeFloat) 
#         print "grid scaling:%f"%gridScaling
        newRD.DoseGridScaling = gridScaling
            
        scaledData = np.zeros(newDose.shape,dtype=typeArray,order='C')
        scaledData[:,:,:] = ((newDose / gridScaling).astype(typeArray))[:,:,:]
        newRD.PixelData = scaledData.tostring()
        
        if not self._RDFileLoaded:
            newRD.FrameOfReferenceUID = self._dcmRSReader.dataForAttribute('FrameOfReferenceUID')
            newRD.StudyInstanceUID = self._dcmRSReader.dataForAttribute('StudyInstanceUID')
            newRD.ManufacturerModelName = self._dcmRSReader.dataForAttribute('ManufacturerModelName')
            newRD.Manufacturer = self._dcmRSReader.dataForAttribute('Manufacturer')
            newRD.SoftwareVersions = self._dcmRSReader.dataForAttribute('SoftwareVersions')
            newRD.PatientName = self._dcmRSReader.dataForAttribute('PatientName')
            newRD.PatientID = self._dcmRSReader.dataForAttribute('PatientID')
            newRD.InstanceCreationTime = self._dcmRSReader.dataForAttribute('InstanceCreationTime')
            newRD.InstanceCreationDate = self._dcmRSReader.dataForAttribute('InstanceCreationDate')
            newRD.StudyDate = self._dcmRSReader.dataForAttribute('StudyDate')
            newRD.StudyTime = self._dcmRSReader.dataForAttribute('StudyTime')
            newRD.PatientBirthDate = self._dcmRSReader.dataForAttribute('PatientBirthDate')
            newRD.ImagePositionPatient = self._dcmCTReader.dataForAttribute('ImagePositionPatient')
        
        newRD.save_as(pathToSave + filename)
        print "RD File: %s saved"%(pathToSave + filename)
        

    def buildPatient(self,dcmCTReader,dcmRDReader,dcmRSReader,stopPwrDataRef ):
        '''Function that "creates" the patient. It builds and initialize all the 3D arrays.
        
        :param dcmCTReader: dicomCTReader giving access to the CT dicom series.
        :param dcmRDReader: dicomRDReader giving access to the RD dicom data. Note: dcmRDReader can be None is the user did not provide a RD dicom file since\
        this dicom file is optional.
        :param dcmRSReader: dicomRSReader giving access to the RS dicom data.
        :param stopPwrDataRef: Stopping power data ( physics.stoppingPower.StoppingPower) to be used for the conversion from density to relative stopping power.
        
        '''
    
    
        print "############### Start building patient data ################"

        #Extract needed data from dcm CT image
        pixelArray = dcmCTReader.dataForAttribute('PixelArray')
        self._spacingInfo = dcmCTReader.dataForAttribute('SpacingInfo') #Note: spacing info in mm
        spacingInfoCT = self._spacingInfo
        imageOrientationCT = dcmCTReader.dataForAttribute('ImageOrientationPatient')
        imagePositionCT = dcmCTReader.dataForAttribute('ImagePositionPatient')
        endImagePositionCT = dcmCTReader.dataForAttribute('EndImagePositionPatient')
        frameDirection = dcmTools.crossProduct( imageOrientationCT[:3] , imageOrientationCT[3:])
        print "Image Orientation CT: %s"%str(imageOrientationCT)
        
        
        #Extract data from RD file
        if self._RDFileLoaded: # If a RD dicom has been provided, load the data
            imagePositionRD = dcmRDReader.dataForAttribute('ImagePositionPatient')
            endImagePositionRD = dcmRDReader.dataForAttribute('EndImagePositionPatient')
            spacingInfoRD = dcmRDReader.dataForAttribute('SpacingInfo')
            pixelArrayRD = dcmRDReader.PixelArray.astype(self._globalVar.typeFloat)
        
        else: # Otherwise use the CT data
            imagePositionRD = imagePositionCT
            endImagePositionRD = endImagePositionCT
            spacingInfoRD = self._spacingInfo
            pixelArrayRD = np.zeros(np.array(pixelArray).shape)

        
        #Find between the CT image and the RD dicom data the one that have the tightest field of view and set the values of "imagePosition" 
        #(i.e. coordinates of the first voxel (0,0,0) of the 3D arrays in the dicom coordinates system) and the values of "endImagePosition"
        #(i.e. coordinates of the last voxel (nb Frames - 1,nb Rows - 1, nbCols -1) of the 3D arrays in the dicom coordinates system)
        self._imagePosition = [0,0,0]
        for k in range(3):
            if imagePositionRD[k] < imagePositionCT[k]:
                self._imagePosition[k] = imagePositionCT[k]
            else:
                self._imagePosition[k] = imagePositionRD[k]
        
        endImagePosition = [0,0,0]
        for k in range(3):
            if endImagePositionRD[k] > endImagePositionCT[k]:
                endImagePosition[k] = endImagePositionCT[k]
            else:
                endImagePosition[k] = endImagePositionRD[k]
        
 
        #Build list of x,y,z values to build the matrices self._XArray, self._YArray ,self._ZArray
        xList = list(np.arange(self._imagePosition[0],endImagePosition[0]+spacingInfoCT['cols'],spacingInfoCT['cols']))
        yList = list(np.arange(self._imagePosition[1],endImagePosition[1]+spacingInfoCT['rows'],spacingInfoCT['rows']))
        zList = list(np.arange(self._imagePosition[2],endImagePosition[2]+spacingInfoCT['frames'],spacingInfoCT['frames']))
  
  
  
        #Set the shape of the 3D arrays. This considers the case where the user want to add padding values around the CT image.
        nbFramesRD = len(zList)
        nbRowsRD = len(yList)
        nbColsRD = len(xList)
        
        padFrames = 0
        padRows = 0
        padCols = 0
        if not self._globalVar.staticDelivery:
            if self._globalVar.ctPadding:
                padFrames = self._globalVar.padValueFRC[0]
                padRows = self._globalVar.padValueFRC[1]
                padCols = self._globalVar.padValueFRC[2]
                
        
        self._nFrames = nbFramesRD + 2 * padFrames
        self._nRows = nbRowsRD + 2 * padRows
        self._nCols = nbColsRD + 2 * padCols

        if not self._globalVar.staticDelivery:
            if self._globalVar.ctPadding:
                print "Padding with %s"%str(self._globalVar.padValueFRC)
                maskPadding = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
                startFrame = padFrames
                startRow = padRows
                startCol = padCols
                maskPadding[startFrame:startFrame + nbFramesRD, \
                            startRow:startRow + nbRowsRD,\
                            startCol:startCol + nbColsRD] = 1
                            
                indPadding = np.where(maskPadding > 0)
                newZ1 = [ zList[0] - i*self._spacingInfo['frames'] for i in range(1,startFrame+1)]
                newZ2 = [ zList[-1] + i*self._spacingInfo['frames']  for i in range(1,padFrames +1)]
                newZ1.reverse()
                zList = newZ1 + zList + newZ2
                
                newY1 = [ yList[0] - i*self._spacingInfo['rows'] for i in range(1,startRow+1)]
                newY2 = [ yList[-1] + i*self._spacingInfo['rows']  for i in range(1,padRows +1)]
                newY1.reverse()
                yList = newY1 + yList + newY2
                
                newX1 = [ xList[0] - i*self._spacingInfo['cols'] for i in range(1,startCol+1)]
                newX2 = [ xList[-1] + i*self._spacingInfo['cols']  for i in range(1,padCols +1)]
                newX1.reverse()
                xList = newX1+ xList + newX2
                
            else:
                indPadding = None
        else:
            indPadding = None
            
        print "Patient data: Num Frames: %i , num Rows: %i, num Cols: %i"%( self._nFrames,self._nRows,self._nCols)
        
        #Initialize all the 3D numpy arrays:
        self._pixelArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
        self._densityArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C') 
        self._stpPwrArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
        self._prescribDoseArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
        self._receivedDoseArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
        self._deff = np.nan*np.zeros((self._nFrames,self._nRows,self._nCols),dtype = self._globalVar.typeFloat,order='C')
        if not self._globalVar.staticDelivery:
            self._tmpPixelArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
            self._tmpDensityArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C') 
            self._tmpStpPwrArray = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
            self._tmpMaskVOIs = []
            
        #Matrices storing the coordinates of the arrays' voxels.
        self._XArray = np.ones((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat)
        self._YArray = np.ones((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat)
        self._ZArray = np.ones((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat)

        [x_mat,y_mat] = np.meshgrid( xList,yList)

        for i in xrange(self._nFrames):
            self._XArray[i,:,:] = x_mat
            self._YArray[i,:,:] = y_mat
        for i in xrange(self._nFrames):
            self._ZArray[i,:,:] = zList[i] * self._ZArray[i,:,:]
        
        
        #Attribute used to store the treatment plan isocenter.
        self._Isocenter = None


        #Indices of the first voxel and the last voxel used for the field of view in the CT image.
        indCTStart = dcmTools.coordinatesToIndexFromImagePosition(self._imagePosition, imagePositionCT, imageOrientationCT, spacingInfoCT)
        indCTEnd =  dcmTools.coordinatesToIndexFromImagePosition(endImagePosition, imagePositionCT, imageOrientationCT, spacingInfoCT)
        
        self._imageOrientation = imageOrientationCT


        ###Fill pixel array CT ############################################################################################
        if not self._globalVar.staticDelivery:
            if self._globalVar.ctPadding: #Padding
                self._pixelArray[:,:,:] =  self._globalVar.paddedCTValue
            
        pixelArray3D = np.array(pixelArray,dtype=self._globalVar.typeFloat)
        indStart = np.array([indCTStart[0] , indCTStart[1] , indCTStart[2]],dtype=self._globalVar.typeFloat)
        shape = np.array([nbFramesRD,nbRowsRD,nbColsRD],dtype=self._globalVar.typeFloat)
        print "CT Shape: %s"%shape
        newPixelArray = np.ascontiguousarray( patientTools.fillPatientCT(pixelArray3D,shape,indStart,self._globalVar.typeFloat),dtype=self._globalVar.typeFloat)
        if indPadding is None:
            self._pixelArray = newPixelArray
        else:
            
            self._pixelArray[startFrame:startFrame + nbFramesRD, \
                            startRow:startRow + nbRowsRD,\
                            startCol:startCol + nbColsRD]= newPixelArray[:,:,:]
         
     
        ### End -fill pixel array CT ############################################################################################
        
        
        ### Fill Dose Grid ############################################################################################
        if self._RDFileLoaded:
            spacingRD = np.array([spacingInfoRD['frames'],spacingInfoRD['rows'],spacingInfoRD['cols']],dtype=self._globalVar.typeFloat) 
            spacingCT = np.array([spacingInfoCT['frames'],spacingInfoCT['rows'],spacingInfoCT['cols']],dtype=self._globalVar.typeFloat)
            dimNewGrid = np.array([nbFramesRD,nbRowsRD,nbColsRD],dtype=self._globalVar.typeFloat)
            doseArray = np.ascontiguousarray(patientTools.fillPatientDose(dimNewGrid,pixelArrayRD,np.array(self._imagePosition,dtype=self._globalVar.typeFloat),np.array(endImagePosition,dtype=self._globalVar.typeFloat),np.array(imageOrientationCT,dtype=self._globalVar.typeFloat),spacingCT,spacingRD,np.array(imagePositionRD,dtype=self._globalVar.typeFloat),self._globalVar.typeFloat),dtype=self._globalVar.typeFloat)
        else:
            doseArray  = pixelArrayRD
        
        if indPadding is None:
            self._prescribDoseArray = doseArray
        else:
            self._prescribDoseArray[startFrame:startFrame + nbFramesRD, \
                            startRow:startRow + nbRowsRD,\
                            startCol:startCol + nbColsRD] = doseArray[:,:,:]
        
        ### End fill dose grid ############################################################################################
        
        
        ### Fill density grid ############################################################################################
        tableConversion = ctToDensConv.getConversionTableCTToMassDens()
        self._densityArray = np.ascontiguousarray(patientTools.fillPatientDensity(self._pixelArray,tableConversion,self._globalVar.typeFloat),dtype=self._globalVar.typeFloat)
        if not self._globalVar.skipSaveDensityImage:
            for axis in self._globalVar.listSaveDensImAlong:
                name = "density_along%s"%(axis.upper())
                viewer2D.storeMatrixAsImage(self._pathToOutput,self._densityArray,name,alongAxis = axis ,cpkl=self._globalVar.saveDensCPKL)
        ### End fill density grid ############################################################################################
              
        
        ### Create masks VOIs ############################################################################################
        self._maskVOIs = self._dcmRSReader.getMasksForAssociatedDicomImage(imagePositionCT, imageOrientationCT, spacingInfoCT, pixelArray3D.shape)

        for voi in self._maskVOIs:
            if voi[4] == "CLOSED_PLANAR":
                indStart = np.array([indCTStart[0] , indCTStart[1] , indCTStart[2]],dtype=self._globalVar.typeFloat)
                
                shape = np.array([nbFramesRD,nbRowsRD,nbColsRD],dtype=self._globalVar.typeFloat)
                reshapedMask= np.ascontiguousarray( patientTools.fillPatientCT(np.asarray(voi[0],dtype=self._globalVar.typeFloat,order="C"),shape,indStart,self._globalVar.typeFloat),dtype='int8')
                if indPadding is None:
                    voi[0] = reshapedMask
                else:
                    tmpMask = np.zeros((self._nFrames,self._nRows,self._nCols),dtype='int8',order='C')
                    tmpMask[startFrame:startFrame + nbFramesRD, \
                            startRow:startRow + nbRowsRD,\
                            startCol:startCol + nbColsRD] = reshapedMask[:,:,:]
                    voi[0] = tmpMask

        #Store all the VOIs as png images
        for voi in self._maskVOIs: 
            if voi[4] == "CLOSED_PLANAR":
                nameROI = ((voi[1].replace(' ','_')).replace('(','_')).replace(')','_')
                obsROI = ((voi[2].replace(' ','_')).replace('(','_')).replace(')','_')
                nameFile = "ROI_%s_%s_%s"%(nameROI,obsROI,str(voi[3]))
                viewer2D.storeMatrixAsImage(self._pathToOutput,voi[0],nameFile, alongAxis = None, cpkl = self._globalVar.storeMasksAsCPKL)
                nbVoxelsVOI = len(np.where(voi[0] == 1)[0])
                print "ROI %s : %i voxels -> %0.3e cm^3."%(nameFile,nbVoxelsVOI,nbVoxelsVOI*self._spacingInfo['frames'] * self._spacingInfo['cols'] * self._spacingInfo['rows'] / 1000.)
                print "\t\t---"
                if voi[2] == 'EXTERNAL':
                    voiExternal = voi[0]
                
            elif voi[4] == "POINT":
                print "ROI point: %s "%str(voi)
            else:
                print "Unknown voi type: %s"%str(voi[4])
                raise ValueError('ROI error : unkown type')

        self.fillIndicesVOIPatient()
        ####End - CReate masks VOIs ############################################################################################
     
        print "############ Patient initialiezd !! ###################"

    
    
    
    
    
    def computeStoppingPowerArray(self,energy = 0):
        ''' Function that computes and fill the relative stopping power 3D array.
        
        :param energy: Energy of the beam being used.
        
        '''
        convTable = self._stopPwrDataRef.getConversionTableStpPwrForEnergy(energy).astype(self._globalVar.typeFloat)
#         tStpPwr0 = time.clock()
        correcFactorTable = self._stopPwrDataRef.getStopPwrCorrecFactor()
        self._stpPwrArray = np.ascontiguousarray(patientTools.fillPatientStopPwr(self._densityArray,convTable,correcFactorTable,self._globalVar.typeFloat),dtype=self._globalVar.typeFloat)
#         tStpPwr1 = time.clock()
#         print "Time Conversion to Stop Pwr Ratio (energy:%f): %s sec"%(energy,str( tStpPwr1 - tStpPwr0 ))    
        if not self._globalVar.skipSaveStpPwrImage:
            for axis in self._globalVar.listSaveStPwrImAlong:
                name = "StoppingPwr_E%0.3fMeV_along%s"%(energy,axis.upper())
                viewer2D.storeMatrixAsImage(self._pathToOutput,self._stpPwrArray,name,alongAxis = axis ,cpkl=self._globalVar.saveStpPwrCPKL)

                    
        
    def setIsocenter(self , isoCoord):
        '''Function used from scanning.deliverySimulation_Static_Dynamic to set the treatment plan isocenter coordinates.
        
        :param isoCoord: treatment plan isocenter coordinates (x,y,z)
        
        '''
        self._Isocenter = isoCoord
        self._IsoInd = None
    
    def getIsocenterInd(self):
        '''Function that computes and return the indices of the voxel representing the isocenter of the treatment plan.
        
        :returns: Indices of the isocenter voxel : [frame, row, columns]
        
        '''
        if self._IsoInd is None :
            if self._Isocenter is None:
                raise ValueError( "Treatment plan isocenter has not been set! Cannot convert to indices")
            self._IsoInd= dcmTools.coordinatesToIndexFromImagePosition(self._Isocenter, self._imagePosition,self._imageOrientation , self._spacingInfo)
        return self._IsoInd

    def shape(self):
        '''
        
        :returns: The patient 3D matrices shape: (number of frames,rows,cols )
        
        '''
        return self._pixelArray.shape

    def getSpacingInfo(self,key = None):
        '''Function to obtain the spacing information used in the patient arrays.
        
        :param key: If None (default), it returns the dictionary with the spacings in mm between the frames ('frames'), rows ('rows') and columns ('cols').\
        if key is not None, it should be 'cols' or 'x_spacing', 'rows' or 'y_spacing', 'frames' or 'z_spacing'.
        
        :returns:  the spacing information corresponding to key. 
        
        '''
        
        if key == 'x_spacing' or key == 'cols':
                return self._spacingInfo['cols']
        elif key == 'y_spacing' or key == 'rows':
                return self._spacingInfo['rows']
        elif key == 'z_spacing' or key == 'frames':
                return self._spacingInfo['frames']
        elif key is None:
            return self._spacingInfo
        else:
            print "Error - getSpacingInfo - no key: "+str(key)+ "found."
            return None
    


    def getUnpaddedPatientDataForAttr(self,attr):
        '''Get the unpadded 3D numpy array of the patient for an attribute. 
        
        :param attr: Attribute for a patient array. See getPatientArrayForAttr(attr) for more information on the attributes.
        
        
        :returns: The unpadded patient array.
        
        '''
        dataPatient = self.getPatientArrayForAttr(attr)
        if self._globalVar.staticDelivery or (not self._globalVar.staticDelivery and not self._globalVar.ctPadding) :
            return dataPatient
        else:
            frames = self._globalVar.padValueFRC[0]
            rows = self._globalVar.padValueFRC[1]
            cols= self._globalVar.padValueFRC[2]
            if attr == 'maskVOIs':
                unpaddedVOIs = []
                for voi in self._maskVOIs:
                    if voi[4] == "CLOSED_PLANAR":
                        shape =  voi[0].shape
                        newMask =voi[0][frames:shape[0]-frames,rows:shape[1]-rows,cols:shape[2]-cols]
                        unpaddedVOIs.append( (newMask , voi[1], voi[2], voi[3], voi[4]) ) 
                    else:
                        unpaddedVOIs.append( voi ) 
                return unpaddedVOIs
            else:     
                shape = dataPatient.shape
            return dataPatient[frames:shape[0]-frames,rows:shape[1]-rows,cols:shape[2]-cols]
    
        

    def getPatientArrayForAttr(self, attr):
        '''Get the 3D numpy array of the static patient for an attribute. 
        
        :param attr: Attribute for a patient array. Attributes are:
            
            * *'x'*: 3D numpy array storing the x values of each voxel in the dicom coordinate system.
            * *'y'*: 3D numpy array storing the y values of each voxel in the dicom coordinate system.
            * *'z'*: 3D numpy array storing the z values of each voxel in the dicom coordinate system.
            * *'density'* : 3D numpy array storing the density of each voxel.
            * *'stpPwrRatio'* : 3D numpy array storing the relative stopping power of each voxel.
            * *'prescribedDose'* : 3D numpy array storing the expected dose distribution obtained from the RD dicom file. \
            If the used did not provide a RD dicom file, it returns None.
            * *'receivedDose'* : 3D numpy array storing the computed dose distribution.
            * *'pixelValues'* : 3D numpy array storing the CT number at each voxel.
            * *'MaskTumor'* : 3D binary array: 1 in the tumor, 0 outside.
            * *'maskVOIs'* :  A list of tuple in which are stored the VOIs information as: ( mask 3D Numpy array , ROI name, ROI observation, ROI index, ROI type)  - Type = 'CLOSED_PLANAR' or 'POINT'
        
        :returns: The patient array.
        
        '''
        
        array = None
        if attr == 'x':
            array = self._XArray
        elif attr == 'y':
            array = self._YArray
        elif attr == 'z':
            array = self._ZArray
        elif attr == 'density':
            array =  self._densityArray
        elif attr == 'stpPwrRatio':
            array = self._stpPwrArray
        elif attr == 'prescribedDose':
            if self._RDFileLoaded:
                array = self._prescribDoseArray
            else:
                array = None
        elif attr == 'receivedDose':
            array = self._receivedDoseArray
        elif attr == 'pixelValues':
            array = self._pixelArray
        elif attr == 'MaskTumor':
            array = self._maskTumor
        elif attr == 'maskVOIs':
            return self._maskVOIs
        else:
            strErr = "Unknown attr: %s, in patient.getPatientArrayForAttr()"% str(attr)
            raise ValueError(strErr)
        
        return array

        
    def getPatientDataForDisplVector(self,attr , dispVec = (0,0,0)):
        '''Get the 3D numpy array of the moving patient for an attribute and a given displacement vector. 
        
        :param attr: Attribute for a patient array. Attributes are:
        
            * *'x'*: 3D numpy array storing the x values of each voxel in the dicom coordinate system. This array is the same if the patient is moving or not.
            * *'y'*: 3D numpy array storing the y values of each voxel in the dicom coordinate system. This array is the same if the patient is moving or not.
            * *'z'*: 3D numpy array storing the z values of each voxel in the dicom coordinate system. This array is the same if the patient is moving or not.
            * *'density'* : 3D numpy array storing the density of each voxel.
            * *'stpPwrRatio'* : 3D numpy array storing the relative stopping power of each voxel.
            * *'pixelValues'* : 3D numpy array storing the CT number at each voxel.
            * *'maskVOIs'* :  A list of tuple in which are stored the VOIs information as: ( mask 3D Numpy array , ROI name, ROI observation, ROI index, ROI type)  - Type = 'CLOSED_PLANAR' or 'POINT'

        
        :param dispVec: Displacement vector in the IEC Fixed coordinate system in cm,  i.e., x_IEC_f = x_Dicom ,y_IEC_f = z_Dicom,z_IEC_f = -y_Dicom 
        
        :returns: The patient array.
        
        '''
            
        if attr == 'x' or attr == 'y' or attr == 'z' :
            return self.getPatientArrayForAttr(attr)
        
        extentStatic = self.getPatientExtentInStaticModel()
        extentDynamic = self.getPatientExtentInDynamicModel(dispVec)
        
                  
        if attr == "density":
            data = self._densityArray
            self._tmpDensityArray[:,:,:] = self._densityArray[:,:,:] 
            newData = self._tmpDensityArray
            valueReplacement = self._densityArray[0][0][0]
        elif attr == "pixelValues":
            data = self._pixelArray
            self._tmpPixelArray[:,:,:] = self._pixelArray[:,:,:]
            newData = self._tmpPixelArray
            valueReplacement = self._pixelArray[0][0][0]
        elif attr == "stpPwrRatio":
            data = self._stpPwrArray
            self._tmpStpPwrArray [:,:,:]= self._stpPwrArray[:,:,:]
            newData = self._tmpStpPwrArray
            valueReplacement = self._stpPwrArray[0][0][0]
        elif  attr == "maskVOIs":
            self._tmpMaskVOIs = []
            for voi in self._maskVOIs:
                self._tmpMaskVOIs.append(voi[:])
                newVOI = self._tmpMaskVOIs[-1]
                if newVOI[4] == "CLOSED_PLANAR":
                    tmpMask = np.zeros((self._nFrames,self._nRows,self._nCols),dtype=self._globalVar.typeFloat,order='C')
                    tmpMask[extentDynamic] = newVOI[0][extentStatic]
                    newVOI[0]= tmpMask
            return self._tmpMaskVOIs
            
        else:
            strErr = "Unknown attr: %s, in patient.getPatientDataForDisplVector()"% str(attr)
            raise ValueError(strErr)            
            
        #Save original values in bounding box
        tmpData = data[extentStatic]
        
        #Fill new values in static bounding box
        newData[extentStatic] = valueReplacement

        #Fill new values in dynamic bounding box (i.e. including the padding)
        newData[extentDynamic] = tmpData
        
        if attr == "density":
            if not self._globalVar.skipSaveDensityImage:
                for axis in self._globalVar.listSaveDensImAlong:
                    name = "density_[%.2f,%.2f,%.2f]_along%s"%(dispVec[0],dispVec[1],dispVec[2],axis.upper())
                    viewer2D.storeMatrixAsImage(self._pathToOutput,newData,name,alongAxis = axis ,cpkl=self._globalVar.saveDensCPKL)
        
        return newData

        
        
        
        
    def getPatientExtentInStaticModel(self):
        '''
        :returns: the patient body indices in the 3D arrays, i.e., the voxels inside the patient body.
        
        '''
        return self._indPatient
        
        
        
        
        
    
    def getPatientExtentInDynamicModel(self, dispVec):
        '''Computes new patient indices in dynamic model, i.e. after a displacement of dispVec.
        
        :param dispVec: Displacement vector in the IEC Fixed coordinate system in cm,  i.e., x_IEC_f = x_Dicom ,y_IEC_f = z_Dicom,z_IEC_f = -y_Dicom 
        
        :returns: the patient body indices in the 3D arrays, i.e., the voxels inside the moving patient body .
        
        '''
        dispIndices = self.getDispIndices(dispVec)
        newIndices = tuple([self._indPatient[0]+dispIndices[0],self._indPatient[1]+dispIndices[1],self._indPatient[2]+dispIndices[2]])
        
        return newIndices
        
        
    def getDispIndices(self, dispVec):
        '''Calculates the displacement vector converted into indices
        
        :param dispVec: Displacement vector in the IEC Fixed coordinate system in cm,  i.e., x_IEC_f = x_Dicom ,y_IEC_f = z_Dicom,z_IEC_f = -y_Dicom 
        
        :returns: Indices in the dicom image (based on the dicom image coordinate system)
        
        '''
        dispVec_mm = [0,0,0]
        for i in range(3):
            dispVec_mm[i] = dispVec[i] * 10.0 + self._imagePosition[i]
        
        (nbFrames,nbRows, nbCols) =  dcmTools.coordinatesToIndexFromImagePosition(dispVec_mm, self._imagePosition, self._imageOrientation, self._spacingInfo)
        
        dispIndices = (int(nbFrames) ,int(nbRows), int(nbCols))
        return dispIndices
 
 
 
 
 
    def computeRadiologicalDepth(self,startVec, incVec, sourceVec, dispVec = (0,0,0) , ener = None):
        '''Computes the radiological depth of the entire patient volume for a given displacement vector.
        
        :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
        :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
        :param sourceVec: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
        :param dispVec: displacement vector in cm in the in the IEC fixed coordinate system.
        :param ener: Beam's energy. Used if the user wants to store the radiological depth every time it is updated.
        
        :returns: 3D numpy array filled with the radiological depth values
        
        '''
        t_start_Deff = time.clock()
        if self._globalVar.staticDelivery:
            stopPwr = self.getPatientArrayForAttr("stpPwrRatio")
        else:
            stopPwr = self.getPatientDataForDisplVector("stpPwrRatio" , dispVec = dispVec)
        reshapedStpPwrGrid = np.ascontiguousarray(np.flipud(np.transpose( stopPwr,(1,0,2))),dtype=self._globalVar.typeFloat) 
        #Note: "Transpose (1,0,2)" changes the slices: slices are not along z dicom axis anymore, but along y dicom axis. "flipud" change the order of the slice (flip upside-down): the 1st becomes last,the last becomes first, the second slice becomes the penultimate and so on... These changes are made such that We can add incVec(positive spacing) to startVec (top left voxel of the first slice) in the IEC FIXED coordinate system to move alonge x,y,z coordinates in the volume until the last voxel(bottom right of the last slice)   
        tmpDeff= patientTools.calcRadiolDepth(reshapedStpPwrGrid, startVec, incVec, sourceVec,self._globalVar.typeFloat)
        radioDepth = np.ascontiguousarray(np.transpose(np.flipud(tmpDeff),(1,0,2)),dtype=self._globalVar.typeFloat)# Undo the initial transformation.
        t_end_Deff = time.clock()
        print "Radiological depth computed in %f sec"%(t_end_Deff - t_start_Deff)
        
        #If the user wants to save the radiological depth as an image
        if not self._globalVar.skipSaveDeff:
            for axis in self._globalVar.listSaveDeffAlong:
                name = "RadiolDepth_[%.2f,%.2f,%.2f]_along%s"%(dispVec[0],dispVec[1],dispVec[2],axis.upper())
                if ener is not None:
                    name = name + "_%iMeV"%ener
                viewer2D.storeMatrixAsImage(self._pathToOutput,radioDepth,name, alongAxis = axis, cpkl = self._globalVar.saveDeffCPKL,bw = False)
        
        return radioDepth
 

    def computeLocalizedRadiologicalDepth(self,startVec, incVec, sourceVec, sourceIECG,rotMatrix,Xr,Yr,Zr,xrshift,yrshift, nbPositionsX = 1,nbPositionsY = 1, dispVec = (0,0,0), ener = None):
        '''Computes the radiological depth only for few rays. All the voxels encountered by the rays will have a value other than -1. The rest of the array will\
        be equal to -1.
        
        :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
        :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
        :param sourceVec: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
        :param sourceIECG: (x,y,z) coordinates of the beam's source in cm in the in the IEC gantry coordinate system.
        :param rotMatrix: Rotation matrix corresponding to  -(gantry angle). It allows to go from the IEC gantry \
        coordinate system to the IEC Fixed coordinate system.
        :param Xr: 3D array of x coordinates in the IEC gantry coordinate system
        :param Yr: 3D array of y coordinates in the IEC gantry coordinate system
        :param Zr: 3D array of z coordinates in the IEC gantry coordinate system
        :param xrshift: beam shift along the x axis of the gantry coordinate system. 
        :param yrshift: beam shift along the y axis of the gantry coordinate system. Note: (xrshift,yrshift) corresponds to the central ray. 
        :param nbPositionsX: Number of positions (i.e. number of rays) along the X axis to the right and to the left of the central ray.
        :param nbPositionsY: Number of positions (i.e. number of rays) along the Y axis to the right and to the left of the central ray.
        :param dispVec: displacement vector in cm in the in the IEC fixed coordinate system.
        :param ener: Beam's energy. Used if the user wants to store the radiological depth every time it is updated.
         
        
        The rays simulated are all the possible combinations of nbPositionsX and nbPositionsY in addition to the central ray. For example, if \
        nbPositionsX = nbPositionsY = 2. All the rays are [xrshift,yrshift], [xrshift + 1 * spacing,yrshift] , \
        [xrshift + 1 * spacing, yrshift + 1 * spacing] ,  [xrshift + 1 * spacing, yrshift - 1 * spacing] , \
        [xrshift , yrshift + 1 * spacing] , [xrshift - 1 * spacing, yrshift + 1 * spacing], [xrshift - 1 * spacing, yrshift - 1 * spacing], \
        [xrshift - 1 * spacing,yrshift] , [xrshift , yrshift - 1 * spacing].....
        
        Spacing = sqrt( incVec[0]^2 + incVec[1]^2)
        
        
        :returns: A tuple:
            
            #. The radiological depth array computed for the rays.
            #. The distance between the source and the patient surface along the central ray
            #. The list of voxels indices encountered in the CT volume along the central ray. The voxels are sorted in increasing distance from the source.
        
        '''
        t_start_beamletRadiolDepth = time.clock()

        minSpacing = np.sqrt( incVec[0] * incVec[0] + incVec[1] * incVec[1])#Note we take only these 2 spacings because the gantry rotates around the 3rd coordinate.
        
        
        if self._globalVar.staticDelivery:
            stopPwr = self.getPatientArrayForAttr("stpPwrRatio")
            density = self.getPatientArrayForAttr("density")
        else:
            stopPwr = self.getPatientDataForDisplVector("stpPwrRatio" , dispVec = dispVec)
            density = self.getPatientDataForDisplVector("density" , dispVec = dispVec)
                
        dynaPatientLoc = self.getPatientExtentInDynamicModel(dispVec)
        patientMask = np.zeros(stopPwr.shape,dtype='int8',order='C')
        patientMask[dynaPatientLoc] = 1
        
        
        auxTargets = []
        for (xPos,yPos) in itertools.product( range(-nbPositionsX,nbPositionsX+1),range(-nbPositionsY,nbPositionsY+1)): 

            if xPos == 0 and yPos == 0:
                targetIECG = np.array([xrshift,yrshift,0],dtype= self._globalVar.typeFloat,order='C')# target point in isocenter frame
                targetIECF =  mathTools.rotateCoordinates(rotMatrix,targetIECG, self._globalVar.typeFloat)      
                ssd0, listIndices,listLenVoxels = self.calculateSSD0(sourceVec,targetIECF,Xr,Yr,Zr,density,startVec,incVec,sourceIECG,patientMask)
                print "SSD0: %f cm"%ssd0

            else:
                pointIECG = np.array([xrshift +  minSpacing * xPos, \
                        yrshift +  minSpacing * yPos,0],dtype= self._globalVar.typeFloat,order='C')# target point in isocenter frame                    
                pointIECF = mathTools.rotateCoordinates(rotMatrix,pointIECG, self._globalVar.typeFloat)
    
                auxTargets.append(pointIECF)
        
        reshapedStpPwrGrid = np.ascontiguousarray(np.flipud(np.transpose( stopPwr,(1,0,2))),dtype=self._globalVar.typeFloat) 
        tmpDeff = patientTools.calcRadiolDepthFewRays(reshapedStpPwrGrid,startVec,incVec,sourceVec, \
                                     np.array(targetIECF,dtype =self._globalVar.typeFloat, order = 'C'),\
                                    np.array(auxTargets,dtype =self._globalVar.typeFloat, order = 'C'), self._globalVar.typeFloat )
        radioDepth = np.ascontiguousarray(np.transpose(np.flipud(tmpDeff),(1,0,2)),dtype=self._globalVar.typeFloat)# Undo the initial transformation.
                                
        t_end_beamletRadiolDepth = time.clock()     
        print "Time to compute  localized radiol depth : %f sec"%(t_end_beamletRadiolDepth - t_start_beamletRadiolDepth)  
               
        if not self._globalVar.skipSaveDeff:
            for axis in self._globalVar.listSaveDeffAlong:
                name = "RadiolDepth_[%.2f,%.2f,%.2f]_along%s"%(dispVec[0],dispVec[1],dispVec[2],axis.upper())
                if ener is not None:
                    name = name + "_%iMeV"%ener
                viewer2D.storeMatrixAsImage(self._pathToOutput,radioDepth,name, alongAxis = axis, cpkl = self._globalVar.saveDeffCPKL,bw = False)
 
        return radioDepth,ssd0,listIndices
        

 
    def calculateSSD0(self,sourceIECF,targetIECF,Xr,Yr,Zr,dens,startVec,incVec,sourceIECG, patientMask):
    
        '''The goal of this function is to find the indices of the voxels crossing the ray starting from source 
        and aiming at target. To obtain the distance source-surface (SSD0), we look along the beam the first voxel
        with a radiological depth > 0.
        
        :param sourceIECF: (x,y,z) coordinates of the beam's source in cm in the in the IEC fixed coordinate system.
        :param targetIECF: (x,y,z) coordinates of the beam's target in cm in the in the IEC fixed coordinate system.
        :param Xr: 3D array of x coordinates in the IEC gantry coordinate system
        :param Yr: 3D array of y coordinates in the IEC gantry coordinate system
        :param Zr: 3D array of z coordinates in the IEC gantry coordinate system
        :param dens: 3D density array. It is used to estimate where the patient body starts in case 'patientMask' is None (i.e. no 'EXTERNAL' structure\
        was present in the RS dicom file)
        :param startVec: Coordinates (x,y,z) in cm of the voxel (0,nb Rows - 1, 0) in the IEC fixed coordinate system.
        :param incVec: Positive increment vector. It corresponds to (x_spacing, y_spacing, z_spacing) in cm.
        :param sourceIECG: (x,y,z) coordinates of the beam's source in cm in the in the IEC gantry coordinate system.
        :param patientMask: 3D binary array : 1 inside the patient body , 0 outside. It can be None if no 'EXTERNAL' structure\
        was present in the RS dicom file)
        
        :returns: A tuple:
            
            #. The source to the patient body surface. If no voxel was found to be considered as the patient surface, the plan isocenter is used and \
            the function returns the distance between the source and the isocenter.
            #. The list of voxels indices encountered in the CT volume along the central ray. The voxels are sorted in increasing distance from the source.
            #. The length of the ray inside each encountered voxel.

        
        '''
        #air_high = 0.001205 # From ICRU report - see stoppingPower.py
        #adipose_high = 0.950
        medDens = 0.01#
        reshapedDensGrid = np.ascontiguousarray(np.flipud(np.transpose( dens,(1,0,2))),dtype=self._globalVar.typeFloat) 
        
        #Note: "Transpose (1,0,2)" changes the slices: slices are not along z dicom axis anymore, but along y dicom axis. 
        #"flipud" change the order of the slice (flip upside-down): the 1st becomes last,the last becomes first, 
        #the second slice becomes the penultimate and so on... 
        #These changes are made such that We can add incVec(positive spacing) to startVec (top left voxel of the first slice) 
        #in the IEC FIXED coordinate system to move alonge x,y,z coordinates in the volume until the last voxel(bottom right of the last slice)   
        
        newTargetIECF = np.array([0,0,0], dtype= self._globalVar.typeFloat, order = 'C')
        newTargetIECF[0] = sourceIECF[0] + ( targetIECF[0] - sourceIECF[0]) * 100.0
        newTargetIECF[1] = sourceIECF[1] + ( targetIECF[1] - sourceIECF[1]) * 100.0
        newTargetIECF[2] = sourceIECF[2] + ( targetIECF[2] - sourceIECF[2]) * 100.0

        
        indices = patientTools.singleRay(reshapedDensGrid,startVec,incVec,sourceIECF, newTargetIECF,self._globalVar.typeFloat)
        [x,y,z] = [None,None,None]
        #Indices: voxel indices from the source to the target
        listIndices = []
        p1 =  sourceIECG
        for idx,(r,f,c,lenVox) in enumerate(indices): # Note: due to the transformations "Transpose" and "flipud" just above on "dens"
            f = int(f)                                   #the order of the indices returned by singleRay() is row,frame,column!
            c - int(c)                                   #And due to flipud row returned must be transformed into dens.shape[1]-row
            r = int((dens.shape[1]-1)-r)
            p2 = [Xr[f,r,c],Yr[f,r,c],Zr[f,r,c]]
            listIndices.append([f,r,c,lenVox,mathTools.distEuclideanP1P2( p1,p2)])
        listIndices = sorted(listIndices,key=itemgetter(4))

        listToReturn = []
        if patientMask is None:
            for idx,(f,r,c,lenVox,dist) in enumerate(listIndices):
                if dens[f,r,c] >= medDens :
                    [x,y,z] = [Xr[f,r,c],Yr[f,r,c],Zr[f,r,c]]
                    ssd0 = dist
                    break
        else:
            for idx,(f,r,c,lenVox,dist) in enumerate(listIndices):
                if patientMask[f,r,c] == 1 :
                    [x,y,z] = [Xr[f,r,c],Yr[f,r,c],Zr[f,r,c]]
                    ssd0 = dist
                    break
        
        if [x,y,z] == [None,None,None]:
            ssd0 = np.sqrt( (sourceIECF[0]-newTargetIECF[0])*(sourceIECF[0]-newTargetIECF[0]) +\
                            (sourceIECF[1]-newTargetIECF[1])*(sourceIECF[1]-newTargetIECF[1]) +\
                            (sourceIECF[2]-newTargetIECF[2])*(sourceIECF[2]-newTargetIECF[2]))
       
                    
        listToReturn = listIndices           
        newListIndices = [ [int(val) for val in item[0:3]] for item in listToReturn ]
        exec "listLenVoxels = [ np.%s(item[3]) for item in listToReturn ]"%self._globalVar.typeFloat
        
        return ssd0, np.array(newListIndices,dtype = 'int', order = 'C'),listLenVoxels


    
    def setReceivedDose(self,computedDose):
        '''Set the "received dose". This should be done, once the delivery simulation ended.
        
        :param computedDose: 3D numpy array - 3D dose distribution
        
        '''
        voiExternal = None
        for voi in self._maskVOIs:
            if voi[4] == "CLOSED_PLANAR":
                if voi[2] == 'EXTERNAL':
                    voiExternal = voi[0]
                    break
        
        if voiExternal is not None:
            tmpDose =  computedDose * np.ascontiguousarray(voiExternal,dtype=self._globalVar.typeFloat)
            
            if not self._RDFileLoaded and np.max(tmpDose) == 0: #If the external is to deep in the image and the radiation did not reach it because of the absence of the RD file.
                self._receivedDoseArray = computedDose
            else:
                self._receivedDoseArray = tmpDose
            
        else:
            self._receivedDoseArray = computedDose
        print "Computed dose distribution set in the patient."
    
        
    def fillIndicesVOIPatient(self):
        '''Function to find the indices of the voxels representing the patient body.\
        It fills the patient attribute '_indPatient'. \
        The Patient body corresponds to the voxels where the value is 1 in the mask of the VOI 'EXTERNAL'. 
        
        
        '''
        voiExternal = None
        for voi in self._maskVOIs:
            if voi[4] == "CLOSED_PLANAR":
                if voi[2] == 'EXTERNAL':
                    voiExternal = voi[0]
                    break
        if voiExternal is not None:
            self._indPatient = np.where(voiExternal == 1)
        else:
            self._indPatient = None
 
         
    def getPatientMask(self):
        '''Get the patient body mask. It corresponds to the VOI 'EXTERNAL' in the RS dicom file.
        
        :returns: 3D binary matrix: 1 in the patient body, 0 outside
        
        '''
        voiExternal = None
        for voi in self._maskVOIs:
            if voi[4] == "CLOSED_PLANAR":
                if voi[2] == 'EXTERNAL':
                    voiExternal = voi[0]
                    break
        return voiExternal    
    
            
        
    def displaySliceOrderingInfo(self):
        '''Function that displays the slice information,i.e. slice number and slice location in mm
        
        '''
        nbSlices = self._pixelArray.shape[0]
        framesSpacing = self.getSpacingInfo('frames')
        print"########### Slice ordering information ##############"
        for idx in xrange(nbSlices):
            print "Slice # %i / %i : %f mm"%(idx+1,nbSlices,self._imagePosition[2]+idx * framesSpacing )
        print"#####################################################"
           
        
        
    def exportSliceOrderingInfo(self):
        '''Export the slice ordering information (i.e. slice number-in first row- and slice location in mm -in second row)to a \
        csv file in the folder /CSV_Data of the current working directory (i.e. 'mainPath'/'nameSimulation') .Note: Rows will become columns if \
        there are more than 255 slices, in order to be able to open the CSV file in excel or numbers.
        '''
        pathToSave = self._pathToOutput + "CSV_Data/"
        
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)     

        filename = "SlicesInformation.csv"
        nbSlices = self._pixelArray.shape[0]
        arrayInfo = np.zeros((2,nbSlices))
        framesSpacing = self.getSpacingInfo('frames')
        for idx in xrange(nbSlices):
            arrayInfo[0,idx] = idx+1
            arrayInfo[1,idx] = self._imagePosition[2]+idx * framesSpacing
            
        if arrayInfo.shape[1] > 255:
            arrayInfo = np.flipud(arrayInfo.transpose((1,0)))
        np.savetxt( pathToSave+filename,arrayInfo,delimiter=',')
        print "\t %s saved."%(pathToSave+filename) 

################################################
#########         Utility functions
################################################

    def hotAndColdSpotsPerStrcuture(self, volRef, volTest, unpadded = False):
        '''Computes the hot and cold spots in the different structures. This function is called from scanning.deliverySimulation_Static_Dynamic for a moving patient. \
        A spot is considered as a hot (cold) spot if its dose is higher (smaller) than 5% of the maximum of volRef.
            
        :param volRef: 3D numpy array. Dose distribution of reference.
        :param volTest: 3D numpy array. Dose distribution to test.
        :param unpadded: True to use unpadded masks, False (default) otherwise.
        
        :returns: Dictionary with the name of the VOIs as keys. The value of each of these key (i.e. for each VOI) is a dictionary whose keys are:
        
            * *'nbVoxels'*: Number of voxels in the VOI
            * *'nbHotSpots'* : Number of hot spots in the VOI
            * *'nbColdSpots'* : Number of cold spots in the VOI
            * *'nbCorrectSpots'* : Number of spots with a 'correct' (i.e. within 5%) dose.
            * *'percHotSpots'* : Percentage of hot spots compared to the entire VOI.
            * *'percColdSpots'* : Percentage of cold spots compared to the entire VOI.
            * *'percCorrectSpots'* : Percentage of correct spots compared to the entire VOI.
        
        '''
        hotAndColdSpots = dict()
        
        if self._RDFileLoaded:
            p =  (5/100)*np.max(volRef) #5% of the maximum reference dose distribution
            diff = simTools(self._globalVar).getMatrixHotAndColdSpot(volRef, volTest)
            
            same = False
            if np.max(diff) > 0:
                same = True
            if unpadded:
                maskVOIs = self.getUnpaddedPatientDataForAttr('maskVOIs')
            else:
                maskVOIs = self._maskVOIs
            for voi in maskVOIs:
                if voi[4] == "CLOSED_PLANAR":
                    nameROI = voi[1]
                    #obsROI = voi[2]
                    mask = voi[0]
                
                    if mask.shape != diff.shape:
                        strErr = "Wrong dimensions in hot and cold spots per structure: diff:%s / mask:%s / pixel array:%s"%(str(diff.shape),str(mask.shape),str(self._pixelArray.shape))
                        raise ValueError(strErr)
                
                    hotAndColdSpots[nameROI] = dict()
                    hotAndColdSpots[nameROI]['nbVoxels'] = len(np.where(mask == 1)[0])
                    if not same:
                        hotAndColdSpots[nameROI]['nbHotSpots'] = len(np.where((mask == 1) & (diff >= p))[0])
                    else:
                        hotAndColdSpots[nameROI]['nbHotSpots'] = 0

                    if not same:                    
                        hotAndColdSpots[nameROI]['nbColdSpots'] = len(np.where((mask == 1) & (diff <= -p))[0])
                    else:
                        hotAndColdSpots[nameROI]['nbColdSpots'] = 0
                    
                    if not same:
                        hotAndColdSpots[nameROI]['nbCorrectSpots'] = len(np.where((mask == 1) & ((diff < p) | (diff > -p)))[0])
                    else:
                        hotAndColdSpots[nameROI]['nbCorrectSpots'] = hotAndColdSpots[nameROI]['nbVoxels']
                        
                    hotAndColdSpots[nameROI]['percHotSpots'] = hotAndColdSpots[nameROI]['nbHotSpots']*100.0/hotAndColdSpots[nameROI]['nbVoxels']
                    hotAndColdSpots[nameROI]['percColdSpots'] = hotAndColdSpots[nameROI]['nbColdSpots']*100.0/hotAndColdSpots[nameROI]['nbVoxels']
                    hotAndColdSpots[nameROI]['percCorrectSpots'] = hotAndColdSpots[nameROI]['nbCorrectSpots']*100.0/hotAndColdSpots[nameROI]['nbVoxels']
                
        return hotAndColdSpots
        
        
    def computeEnergyReceivedByStructures(self,doseData, filename):
        '''Computes the energy in Joules received by each patient structure. The data is exported to a .txt file save in the simulation folder.
        
        :param doseData: 3D numpy array representing the dose distribution.
        :param filename: Name of the file to save. Filename is preceded of 'SummaryEnergy' when the file is being written.
            
        '''
    
        pathToSave = self._pathToOutput 
            
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)
        meV = self._globalVar.MeV #J
        nbProt = self._globalVar.protPerMU
        with open(pathToSave + "SummaryEnergy_"+"%s"%filename+'.txt', 'w') as f:
            f.write('Summary (%s) - energy repartition, per fraction, in patient body:\n'%filename)
            vois = self.getUnpaddedPatientDataForAttr('maskVOIs')
            spacing = self.getSpacingInfo()
            densityMat = self.getUnpaddedPatientDataForAttr('density')  
            for voi in vois:
                if voi[4] == "CLOSED_PLANAR":
                    nameROI = ((voi[1].replace(' ','_')).replace('(','_')).replace(')','_')
                    obsROI = ((voi[2].replace(' ','_')).replace('(','_')).replace(')','_')
                    nameFile = "ROI_%s_%s_%s"%(nameROI,obsROI,str(voi[3]))
                    energyVoiMat = 1e-3 * densityMat * voi[0] * doseData * spacing['frames'] * spacing['rows'] * spacing['cols'] * 1e-3
                    energyVoi = np.sum(energyVoiMat)
                    
                    f.write('\t %s :\t %f \t Joules,  \t %f \t MeV ,  \t %f \t MeV(Fluence: 1e-9 Prot/MU) \n'%( nameROI, energyVoi,energyVoi/meV,energyVoi/(meV*nbProt)))
        
