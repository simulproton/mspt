########################################################################
#
# simulatorMSPT.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# April 2012
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
########################################################################
import sys,os,csv
import traceback
import numpy as np
import scipy as sp
import cProfile
import socket
import time, datetime

from ..dicomReader.dicomManager import DicomManager
from ..patient.patient import Patient
from ..physics.stoppingPower import StoppingPower as stpw
from ..physics.ddCurvesManager import DDCurves
from ..motion.motionManager import Motion
from ..motion.motionMonitor import MotionMonitor
from ..dataViewer import array2DViewer as viewer2D
from ..dataViewer import plotDataGraphs as pltData
from ..dataViewer import dvh
from ..scanning.deliverySimulation_Static_Dynamic import simulatorDelivery
from tools import Tools as simTools

from  ..email import Email
from  ..email.Email import Email as emailManager
from ..zip import zipFiles

from globalVariables import GlobalVariables


if os.name == 'nt': #windows
    cpCmd = 'copy'
    scpCmd = 'copy'
    rmCmd = 'del'
else:
    cpCmd = 'cp'
    scpCmd = 'scp'
    rmCmd = 'rm'


class MSPT(object):
    '''The class "MSPT",  is the main class of the simulator. It controls the global simulation process.
    
    :param params: List of the input files of the simulator:
    
        * The path to a CT dicom directory in which is contained the CT dicom series.
        * The path to the RP dicom file (treatment plan)
        * The path to the RS dicom file (structures defined in the CT)
        * (optional) The path to the RD dicom file (expected dose distribution)
        * (optional) The path to the file configuring MSPT. If no file is provided MSPT will run a simulation on a static patient and \
        MSPT will be configured using default settings (see simulator.globalVariables). If provided, the configuration file is received by MSPT \
        but the configuration can be performed before. In such case 'globalVars' should not be None. 
        
    :param globalVars: Dictionary containing the configuration data of MSPT. If None: 
    
        * If the configuration file is provided: MSPT will be configured using this file. The missing parameters will be set to default.
        * If the configuration file is not provided: SPT will be configured with default settings: static delivery and output data stored to local directory.
        
    .. note::
    
        The configuration file can either be '.txt' or '.py'. In this file do not precede variables with an indentation.\
        Enter 1 variable per line: myVariable = Value . All the following variables can be set in this file, but they are not mandatory since they all have\
        a default value:
        
            * variables common to static and dynamic deliveries:
                
                * **mainPath**: String defining the path where the simulation directory will be created. By default: '../Data_Simulation/'
                * **nameSimulation** : corresponds to the path where the simulation output data will be saved. Default: 'runSimulStatic/'' for a static\
                delivery, 'runSimulDynamic/' for a dynamic simulation.
                * **nameProfilage** : name of the profilage file created using the package 'cProfile'. This profilage data can be after studied using the cProfile package in python or \
                with the profilage viewer `runsnake <http://www.vrplumber.com/programming/runsnakerun/>`_ . Default: 'profilagerunSimulStatic.profile' for a static delivery, \
                'profilagerunSimulDynamic.profile' for a dynamic delivery.
                * **removeProfilageData** : True (default) or False depending if the user wants to delete the profilage data when the simulation ends. 
                * **stdOutErrToFile** : True or False (default). If True the standard output and error are redirected to a text file \
                in the simulation directory.
                * **zipSavedData** : True (default) or False. If True, the simulation directory containing the output data is zipped in 'maintPath'. \
                Note: if the user access directly the class MSPT (not using the '__main__' function) the zip option is not available. 
                * **copyZip** : True (default) or False. If True the zip file is copied to the 'pathToCopyZip' directory.
                * **pathToCopyZip** : String defining where to copy the zip file. Default: '../Data_Zip/'
                * **typeFloat** : The type of numpy float to be used: 'float32' or 'float64'. Default : 'float64'.
                * **emailUser** : If the user wants to send an email saying that the simulation ended with the DVH(s) and the the data printed in the \
                standard output (or if the simulation crashed for any reason), he can fill this variable with his e-mail address en fill the next fields. \
                Note: if the user access directly the class MSPT (not using the '__main__' function) the email option is not available.\
                Default: None
                * **emailPwd** : e-mail password. Default: None
                * **emailSMTP** : e-mail smpt server address. Default: None
                * **emailSMTPPort** : smtp server port. Default: None
                * **emailRecipient** : e-mail address of the e-mail's recipient. Default: None
                * **staticDelivery**: True (default) or False. If True the delivery is performed for a static patient. Otherwise for a moving (dynamic)\
                patient.
                * **ddCurves**: String defining the set of dept-dose curves to use. By default 'RayStation'. It is also possible to use 'MCNPX' or 'Eclipse', \
                but these two datasets have not been tested. The user can also add its own depth dose curves. See physics.ddCurvesManager.DDCurve for more information.
                * **displayScanPathUpdates** : True or False (default). If True it displays the evolution of 2D array representing the scanning path\
                in an energy layer. See dicomReader.scanningPath for more details.
                * **saveCTImage** : True :save the CT dicom image as tiff image files, False otherwise. Default: True
                * **listSaveCTImageAlong** (if 'saveCTImage' is True): list of the axis along which the user wants to save the CT images of the patient. Possible values are\
                'x' (i.e. the images are stored along the column axis of the dicom image), 'y' (i.e. the images are stored along the rows axis of the dicom image), \
                'z' (i.e. the images are stored along the frame axis of the dicom image - default). Examlpes: listSaveCTImageAlong = [] (in this case default value is used) , or listSaveCompDoseAlong = ['x'] or listSaveCompDoseAlong = ['x','z'] ... \
                Default: ['z'].
                * **skipSaveDensityImage** : False to save the 3D density array as a set of png images. True otherwise. Default: True
                * **listSaveDensImAlong** (if 'skipSaveDensityImage' is False): list of the axis along which the user wants to save the density \
                images of the patient. Default: ['z'].
                * **saveDensCPKL** (if 'skipSaveDensityImage' is False): True if the user wants to store the density 3D array in a cPickle structure (this allows to have access to the \
                3D numpy array outdise the simulation). False otherwise. Moreover, if False and the patient is moving\
                it will stored the 3D density array every time the array is updated. Default: False.
                * **skipSaveStpPwrImage** : False to save the 3D relative stopping power array as a set of png images. True otherwise. Default: True
                * **listSaveStPwrImAlong** (if 'skipSaveStpPwrImage' is False): list of the axis along which the user wants to save the rel. stop. power \
                images of the patient. Default: ['z']. 
                * **saveStpPwrCPKL** (if 'skipSaveStpPwrImage' is False):Similar to 'saveDensCPKL' to save the 3D array of the relative stopping power. However, if the patient is moving it only stores it once.
                * **skipSaveDeff** : False to save the 3D radiological depth array as a set of png images. True otherwise. Default: True
                * **listSaveDeffAlong** : Similar to 'listSaveDensImAlong', to save the images of the 3D radiological depth matrix.
                * **saveDeffCPKL** : Similar to 'saveDensCPKL' to save the 3D radiological depth matrix. If the patient is moving it stores the array at every update.
                * **skipSaveBeamlets** : False to save an image of the dose distribution of each beamlet simulated, True otherwise. Default : True
                * **saveBeamletsCPKL** (if skipSaveBeamlets is False): True to save the dose distribution of each beamlet (a 3D numpy array) in a cPickle structure.
                * **listSaveCompDoseAlong** : list of the axis along which the user wants to save the computed dose \
                images (png) of the patient. Default: ['z']. Note: This variable is also used to store the density and the pencil beam dose distribution in a \
                dynamic delivery if the variable 'storeInterDoseDens' (defined below) is True.
                * **saveCompDoseAsDcm** : True to save the computed dose distribution into a RD dicom file. False otherwise. Default : True
                * **nameCompDose** : Name of the output RD dicom file for the computed dose distribution. Default: 'RD_CompDose_MSPT.dcm' 
                * **saveRefDoseAsDcm** : True to save the expect dose distribution into a RD dicom file  (False otherwise). This is useful if the input dose grid\
                must be resampled to the voxel spacing of the input CT dicom series. Moreover, this can only be done if the user provided an inpu RD dicom file. Default : True
                * **nameRefDose** : Name of the output RD dicom file for the expected dose distribution. Default: 'RD_RefDose.dcm'
                * **saveRefDose** : True :save the expected dose distribution image as png image files, False otherwise. Default: True
                * **listSaveRefDoseAlong** (if saveRefDose True) : list of the axis along which the user wants to save the expected dose distribution images.
                * **dvhRelativeVolume** : True to convert the volumes into relative volumes in the dose volume histograms. False otherwise.  Default: True.
                * **displaySliceOrderingInfo** : True to display the order and positions (in mm) of the CT dicom slices. False otherwise. Default: False.
                * **exportSliceOrderingInfo** : True to export the order and positions (in mm) of the CT dicom slices into a cszv file. False otherwise. Default: False.
                * **storeMasksAsCPKL** : True to store binary matrices representing the patient structures (VOIs: Volumes of Interest) into cPickle structures.
                * **protPerMU** : Number of protons per Monitor Units (MU). Default: 1e9 prot.MU-1.
                * **nameStopPwrCorrecFactor** : name of the correction factor for the stopping power to use.
                * **importNewMassStopPwr** : True to import new mass stopping power data, False (default) otherwise.
                * **nameDoseCorrecFactor**: name of the dose correction factor to use.
            
            * variables for dynamic deliveries:
                
                * **unitTimeForDelivery** : time unit to simulate the dynamic delivery (moving patient) expressed in seconds. It corresponds to the time resolution. Default: 1e-3.
                * **typeMotion** : the type of motion applied to the patient: 'motionNoisyCos','motionCos','motionBreathHiccup' or 'motionBreathCough'. See package motion for more information. \
                Default: 'motionCos'.
                * **arrayPeriod** : List of the motion period (in seconds) along x,y, and z axis in the IEC fixed coordinate system. \
                Default : [4.0,  3.5,   0.0].Note: X_IEC fixed= X_Dicom, Y_IEC fixed = Z_Dicom, Z_IEC fixed = -Y_Dicom.               
                * **arrayAmpl** : List of the motion amplitude (in cm) along x,y, and z axis in the IEC fixed coordinate system. \
                Default : [1.0,  1.0,  0.0]
                * **arrayPhase** : List of the motion initial phase (in rad.) along x,y, and z axis in the IEC fixed coordinate system. \
                Default : [0,  0,  0]
                * **arrayStDevPeriod** (if 'typeMotion' is 'motionNoisyCos', 'motionBreathHiccup' or 'motionBreathCough'): List of standard deviations (in seconds) \
                for the gaussian distribution used to simulate changes (noise) in period of the motion along x,y, and z axis in the IEC fixed coordinate system. \
                Default : [0,  0,  0]
                * **arrayStDevAmpl** (if 'typeMotion' is 'motionNoisyCos', 'motionBreathHiccup' or 'motionBreathCough'): List of standard deviations (in cm)\
                for the gaussian distribution used to simulate changes (noise) in amplitude of the motion along x,y, and z axis in the IEC fixed coordinate system. \
                Default : [0,  0,  0]
                * **ctPadding** : True if the user wants to add margins of a specified density (see 'paddedCTValue' ) around\
                the current CT image. This is useful, when the moving patient can be outside of the bounding box formed by the borders of the CT image.  
                * **padValueFRC** (if 'ctPadding' is True) : Numpy array of 3 integer values : [ number of added frames : nbF, number of added rows : nbR, \
                number of added columns : nbC]. For exemple : [ 0,0,10] will add 10 columns to the right of the CT image and 10 columns to the left of the CT image.\
                Default: [10,10,10].
                * **paddedCTValue** (if 'ctPadding' is True): The value in CT number (Hounsfield Units) to use to pad the CT image. For example: air = -1000 HU, water = 0 HU.
                * **storeInterDoseDens** : True to store the density image of the moving patient and the dose distribution of the beam when the \
                patient is moving while being irradiated by a pencil beam. False otherwise. Default: False. Note: It uses the variable \
                'listSaveCompDoseAlong' to save the images.
                * **motionDuringBeam** :  True to allow the patient to move while the beam is irradiating, False otherwise. Default : True.
                * **compensationDynamic** : True to allow motion compensation, False otherwise. Default : False.
                * **measurementStdev** (if compensationDynamic is True):  List of the standard deviation of the motion monitoring accuracy (in cm.) \
                along x,y, and z axis in the IEC fixed coordinate system. Default : [0.1,  0.1,  0]
                * **measurementAvg** (if compensationDynamic is True):List of the average of the motion monitoring accuracy (in cm.) along x,y, and z axis in the IEC fixed\
                coordinate system. Default : [0.1,  0.1,  0].
                * **timeMeasureDisplacement** (if compensationDynamic is True): Time (in seconds) needed to measure the patient displacement with the monitoring system. Default: 1e-3.
                * **findStartPos** (if compensationDynamic is True): 
                * **nameNewRPFile** (if compensationDynamic is True):
                * **addMargins** (if compensationDynamic is True):
                * **marginsParam** (if compensationDynamic is True):
                * **updateMargins** (if compensationDynamic is True):
                * **unlimitedRescanning** : True to allow an unlimited number (the maximum, however, is set to 32767 repaintings to avoid infinite loops) of repainting per energy slice.\
                False otherwise.
                * **repaintingFactor** (if unlimitedRescanning is False): Repainting factor. Default: 1. 
                * **sliceRepaint** : Type of slice repainting method. Types :'IsolayerRepaint','ScaledRepaint' or 'FullWeightRepainting'.\
                Note: 'FullWeightRepainting': we deliver the full weight at once instead of a maximum weight or a scaled weight.
                * **volumeRepaint** : Type volume repainting method. Types : 'Volumetric' and 'NonVolumetric'. Default : 'NonVolumetric'.                
                * **evaluateDiffStatDyna** : True to show and store statistics about the differences between dose distributions  at a beam position when \
                the patient is moving and when he is static. Default: False.
                * **synchrotron** : True to make MSPT behave like a synchrotron (i.e., pulsatile proton delivery, time needed to refill the synchrotron and the \
                time to change the energy is longer than for cyclotrons). Otherwise (False) it behaves like a cyclotron (i.e., continuous proton delivery\
                 the time to change the energy is shorter. Default: True
                * **spillLength** ( if synchrotron is True): length in seconds of a proton spill. Default: 5.
                * **timeEnergyChange** : Time (in seconds) to change the energy. Default: 2.1 ( if synchrotron is True) , 50e-3 ( if synchrotron is False).
                * **spotSettlingTime** : Time (in seconds) needed to stabilize the beam position. Default : 5e-3.
                * **lateralScanningSpeed** : Lateral scanning speed of the beam (in cm/s). Default : 1000.
                * **verticalScanningSpeed** : Vertical scanning speed of the beam (in cm/s). Default : 1000.
                * **speedGantryRotation** : Gantry rotational speed in degrees / min. Default: 360.
                * **beamIntensity** : Beam intensity (or extraction rate), i.e. number of protons per second being delivered. Default: 6.8e10
       
    '''

    def __init__(self, params , globalVars = None):
        '''
      
        '''
        paramsDicom = []
        if len(params) != 5 and len(params) != 4 and len(params) != 3:
            strErr = "Error number of arguments for the simulator..! 3,4 or 5 are expected, %i are given."%len(params)
            raise ValueError(strErr);
        else:            

            print 'Starting MSPT....'
            configFile = None
            for item in params:
                if item.endswith(".py") or item.endswith(".txt"):
                    if globalVars is None:
                        configFile = item
                        globalVars = GlobalVariables(configFile)
                elif item == '--motion':
                    if globalVars is None:
                        with open('tmp_ConfigData.txt', 'w') as myFile:
                            myFile.write("staticDelivery = False \n")
                        globalVar = GlobalVariables('tmp_ConfigData.txt')
                        os.system("%s tmp_ConfigData.txt"%rmCmd)
                else:
                    paramsDicom.append(item)
            
            if configFile is None and globalVars is None:
                globalVars = GlobalVariables(None)
            
            self._globalVar = globalVars 

            print ">>>>>>>>>>> Beam Scanning Proton Machine Settings <<<<<<<<<<"
            print "# of protons/MU: %0.2e"%self._globalVar.protPerMU
            if not self._globalVar.staticDelivery:
                if self._globalVar.synchrotron:
                    print "Machine Type: Synchrotron"
                else:
                    print "Machine Type: Cyclotron"
                print "Extraction rate: %0.2e protons/s"%self._globalVar.beamIntensity
                print "Time to change the energy: %0.2f s"%self._globalVar.timeEnergyChange
                if self._globalVar.synchrotron:
                    print "Spill length: %0.2f s"%self._globalVar.spillDuration
                print "Lateral scanning speed: %0.2f m/s"%(self._globalVar.lateralScanningSpeed/10)
                print "Vertical scanning speed: %0.2f m/s"%(self._globalVar.verticalScanningSpeed/10)
                print "Gantry rotational speed: %0.2f degrees/min"%self._globalVar.speedGantryRotation
            print ">>>>>>>>>>> End Machine Settings  <<<<<<<<<<"
            
            
            if not self._globalVar.staticDelivery:
                print ">>>>>>>>>>> MSPT Settings for Motion <<<<<<<<<<"
                print "Unit time for delivery: %e s"%self._globalVar.unitTimeForDelivery
                if self._globalVar.compensationDynamic:
                    print "Time for motion measurement: %e s"%self._globalVar.timeMeasureDisplacement
                print ">>>>>>>>>>> End MSPT Settings for Motion  <<<<<<<<<<"
            
            pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation
            if not os.path.exists(pathToSave):
                os.makedirs(pathToSave)

            print "Starting simulation %s"%self._globalVar.nameSimulation
            dcmManager = DicomManager(paramsDicom,self._globalVar.typeFloat)
            
            # Create patient from dicom data
            dicomDataForPatient = dcmManager.getDicomReadersForPatient()
            
            #Initialize stopping power ratios data
            self._dataRefStopPwr = stpw(self._globalVar)
            
            
            self._patient = Patient(self._dataRefStopPwr,dicomDataForPatient,self._globalVar)
            
            # get Scanning Path , number of fractions,:
            dicomDataPlan = dcmManager.getDicomReaderForScanningPath()
            numberOfFractionsPlanned = dicomDataPlan.dataForAttribute('NumberOfFractionsPlanned')
            spotPositions = dicomDataPlan.getScanningPath(self._globalVar.settingsScanning)
            dicomDataPlanPath = dicomDataPlan.rtPlanPath
            print "\tTreatment plan loaded !"
            
            #Initialize depth dose curves
            self._ddCurvesData = DDCurves(self._globalVar)
            
            #Simulate treatment
            if self._globalVar.staticDelivery: #Static delivery
                simulator = simulatorDelivery(self._patient,self._ddCurvesData,self._globalVar)
                simulator.simulateDelivery(spotPositions)
                
                
            else: #Dynamic delivery
                motion = Motion(self._globalVar.motionInfo,self._globalVar.typeFloat)
                if self._globalVar.compensationDynamic:
                    monitoringSystm = MotionMonitor(motion, self._globalVar.measurementStdev,self._globalVar.measurementAvg)
                else:
                    monitoringSystm = None
                simulator = simulatorDelivery(self._patient,self._ddCurvesData,self._globalVar, motion = motion , monitor =monitoringSystm )
                compensatedPlan = simulator.simulateDelivery(spotPositions)
                #Note: compensatedPlan is None if there is no compensation
                
            
            #Get computed dose.
            simulatedDose = self._patient.getUnpaddedPatientDataForAttr('receivedDose')
            #Store results - total Dose
            

            if self._globalVar.saveCTImage:
                ctImage = self._patient.getUnpaddedPatientDataForAttr('pixelValues')
                for axis in self._globalVar.listSaveCTImageAlong:
                    viewer2D.storeMatrixAsImage(pathToSave,ctImage,"DcmImage",alongAxis = axis,bw = True)
                spacingList = [['Cols','Rows','Frames'],[self._patient.getSpacingInfo('x_spacing')/10.0,self._patient.getSpacingInfo('y_spacing')/10.0,self._patient.getSpacingInfo('z_spacing')/10.0]]
                self.save2DListToCSV(spacingList,"SpacingData_in_cm")

            if not self._globalVar.skipSaveDensityImage:
                densImage = self._patient.getUnpaddedPatientDataForAttr('density')
                for axis in self._globalVar.listSaveDensImAlong:
                    viewer2D.storeMatrixAsImage(pathToSave,densImage,"DensityImage",alongAxis = axis)


        
            # Get the planned dose matrix
            prescribedDoseLoaded = True
            doseDistrOfRef = self._patient.getUnpaddedPatientDataForAttr('prescribedDose')
            if doseDistrOfRef is None:
                 prescribedDoseLoaded = False
            if prescribedDoseLoaded and numberOfFractionsPlanned > 1:
                doseDistrOfRef = doseDistrOfRef / numberOfFractionsPlanned

            if not self._globalVar.skipSaveStpPwrImage:
                stpPwrRatioImage = self._patient.getUnpaddedPatientDataForAttr('stpPwrRatio')
                for axis in self._globalVar.listSaveStPwrImAlong:
                    viewer2D.storeMatrixAsImage(pathToSave,stpPwrRatioImage,"stpPwrImage")

            if prescribedDoseLoaded:#Clip smallest dose values in computed dose to 0.
                tool = simTools(self._globalVar)
                minRef = np.min( doseDistrOfRef[np.where(doseDistrOfRef > 0)])
                indSmallVlaues = np.where((simulatedDose < minRef) & (simulatedDose > 0))
                simulatedDose[indSmallVlaues] = 0
            
                for idx,frame in enumerate(doseDistrOfRef):
                    if np.max(frame) > 0:
                        self.save2DSliceToCSV(frame,"RefDose_Frame"+str(idx))
            
                if self._globalVar.saveRefDoseAsDcm:
                    self._patient.saveDoseDistrAsRDDcmFile(self._patient.getUnpaddedPatientDataForAttr('prescribedDose'),self._globalVar.nameRefDose)
            
            
            #Save computed dose
            if self._globalVar.saveCompDoseAsDcm:
                self._patient.saveDoseDistrAsRDDcmFile(self._patient.getUnpaddedPatientDataForAttr('receivedDose'),self._globalVar.nameCompDose)

            
            for idx,frame in enumerate(simulatedDose):
                if np.max(frame) > 0:
                    self.save2DSliceToCSV(frame,"simulatedDose_Frame"+str(idx))
            
            if prescribedDoseLoaded:
                #Find min/max values to use to display simulated dose and ref dose in same color scale.
                min = np.min( [np.min(simulatedDose),np.min(doseDistrOfRef)] )
                max = np.max( [np.max(simulatedDose),np.max(doseDistrOfRef)] )
            else:
                min = np.min(simulatedDose)
                max = np.max(simulatedDose)
                
            
            for axis in self._globalVar.listSaveCompDoseAlong:
                viewer2D.storeMatrixAsImage(pathToSave,simulatedDose,"SimulatedDose",alongAxis = axis,minmax=[min,max],bw = False)
                
            if self._globalVar.saveRefDose and prescribedDoseLoaded:
                #Store prescribed dose  
                for axis in self._globalVar.listSaveRefDoseAlong:
                    viewer2D.storeMatrixAsImage(pathToSave,doseDistrOfRef,"RefDose",alongAxis = axis,minmax=[min,max],bw = False)

                
            if prescribedDoseLoaded: 
                self._patient.computeEnergyReceivedByStructures(doseDistrOfRef, "Reference_DoseDistribution") 
            self._patient.computeEnergyReceivedByStructures(simulatedDose, "Simulated_DoseDistribution") 
            
            if prescribedDoseLoaded:
                #Comparision planned dose and computed dose
                print"\n##################################"
                print"Evaluating Differences :  doseDistrOfRef  /  simulatedDose"
                tool.evaluateDifference(doseDistrOfRef,simulatedDose)
                print"##################################\n"
            
            
                #Show Hot and cold spots
                tool.showHotColdSpots(doseDistrOfRef,simulatedDose, title = "Hot and Cold Spots" , pathToSaveImg = pathToSave + 'HotColdSpotsBeamlets/'+'HotAndColdSpots_Final')
            
                diffDose = tool.getMatrixHotAndColdSpot(doseDistrOfRef,simulatedDose)
                #Save Cold Spots
                matColdSpots = np.ones(diffDose.shape)
                indCold = np.where(diffDose<0)
                if len(indCold[0]) > 0:
                    maxCold = np.max(diffDose[indCold])
                    matColdSpots[indCold] = 1 - (diffDose[indCold]/ maxCold)
                    for axis in self._globalVar.listSaveCompDoseAlong:    
                        viewer2D.storeMatrixAsImage(pathToSave,matColdSpots,"ColdSpots-FinalVolume",alongAxis = axis)
                else:
                    print "No differences found between expected dose and computed dose, therefore cold spot will not be plotted" 
                #Save Hot Spots
                matHotSpots = np.zeros(diffDose.shape)
                indHot = np.where(diffDose>0)
                if len(indHot[0]) > 0:
                    maxHot = np.max(diffDose[indHot])
                    matHotSpots[indHot] = (diffDose[indHot]/ maxHot)
                    for axis in self._globalVar.listSaveCompDoseAlong:    
                        viewer2D.storeMatrixAsImage(pathToSave,matHotSpots,"HotSpots-FinalVolume",alongAxis = axis)
                else:
                    print "No differences found between expected dose and computed dose, therefore hot spots will not be plotted" 
                        
                        
                        
                dataHotColdSpots = self._patient.hotAndColdSpotsPerStrcuture(doseDistrOfRef,simulatedDose, unpadded = True)
                
                with open(pathToSave + 'EvalHotColdSpots.txt', 'w') as myFile:
                    myFile.write("Evaluation Hot and Cold spots :\n")
                    myFile.write("Note: Hot (resp. cold) spots are considered as such, if the absolute difference is greater than 5% of the maximum dose in the expected dose distribution.\n\n\n")

                    for roi in dataHotColdSpots:
                        myFile.write("%s :\n\t#of voxels: %i , \n\t#of hot spot:%i ( %0.3f %%) , \n\t#of cold spots:%i (%0.3f %%), \n\t#of correct spots:%i (%0.3f %%)\n"%(roi,dataHotColdSpots[roi]['nbVoxels'],\
                        dataHotColdSpots[roi]['nbHotSpots'],\
                        dataHotColdSpots[roi]['percHotSpots'],\
                        dataHotColdSpots[roi]['nbColdSpots'],\
                        dataHotColdSpots[roi]['percColdSpots'],\
                        dataHotColdSpots[roi]['nbCorrectSpots'],\
                        dataHotColdSpots[roi]['percCorrectSpots']))
                
                
                
                
            if self._globalVar.displaySliceOrderingInfo:
                self._patient.displaySliceOrderingInfo()
                
                
            if self._globalVar.exportSliceOrderingInfo:
                self._patient.exportSliceOrderingInfo()
            
            
            listVOIs = self._patient.getUnpaddedPatientDataForAttr('maskVOIs')
            
            # Build DVH for simulated dose
            listDv = [100,99,98,97,96,95,90,66.66,50,33.33,10,5,4,3,2,1]
            listVd = []
            nameDVH = 'SimulatedDose'
            dvh.processDVHForDoseDistribution(simulatedDose, listVOIs , listDv, listVd, nameDVH,self._globalVar)
            if  prescribedDoseLoaded:
                nameDVH = 'ReferenceDose'
                dvh.processDVHForDoseDistribution(doseDistrOfRef, listVOIs , listDv, listVd, nameDVH,self._globalVar)

             
            #Plot and export motion history
            if not self._globalVar.staticDelivery: 
                historyMotion = motion.getHistory()
                listX = [np.array(historyMotion['time']),np.array(historyMotion['time']),np.array(historyMotion['time'])]
                listY = [np.array(historyMotion['x']),np.array(historyMotion['y']),np.array(historyMotion['z'])]
                listNames = []
    
                for axis in ['X','Y','Z']:
                    nameCurve = "Motion along %s axis"%axis
                    listNames.append(nameCurve)
                namePlot = "Motions along X Y Z axis vs time"
                xLabel = "Time (s)"
                yLabel = "Displacement (cm)"
                (figure,subplot) = pltData.drawData(listX,listY,listNames = listNames, display = False, title =namePlot, xLabel = xLabel, yLabel =  yLabel,posLegend =None )
                pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation + '/'
                pltData.savePlotImage(figure, namePlot.replace(" ","_") , dirPath = pathToSave, ext = 'png', subplot = None)
                pltData.closePlot()
                listData = [np.array(historyMotion['time']),np.array(historyMotion['x']),np.array(historyMotion['y']),np.array(historyMotion['z']) ]
                names = ['time(s)','X(cm)','Y(cm)','Z(cm)']
                exportListDataToCSV(listData,names,pathToSave,namePlot.replace(" ","_"))
                print "Motion history exported"
                
                if self._globalVar.compensationDynamic:
                    newRTPlanPath = pathToSave + self._globalVar.nameNewRPFile
                    if compensatedPlan is not None:
                        tool.saveCompensatedRTPlan(dicomDataPlanPath, newRTPlanPath, compensatedPlan)
                    historyMonitoring =  monitoringSystm.getHistory()
                    listX = [np.array(historyMonitoring['time']),np.array(historyMonitoring['time']),np.array(historyMonitoring['time'])]
                    listYMeasure = [np.array(historyMonitoring['x']['MeasuredMotion']),np.array(historyMonitoring['y']['MeasuredMotion']),np.array(historyMonitoring['z']['MeasuredMotion'])]
                    listYMotion = [np.array(historyMonitoring['x']['PatientMotion']),np.array(historyMonitoring['y']['PatientMotion']),np.array(historyMonitoring['z']['PatientMotion'])]
                    listData = [np.array(historyMonitoring['time'])]
                    listData.extend(listYMeasure)
                    listData.extend(listYMotion)
                    names = ['time(s)','X_Measure(cm)','Y_Measure(cm)','Z_Measure(cm)']
                    names.extend(['X_Motion(cm)','Y_Motion(cm)','Z_Motion(cm)'])
                    exportListDataToCSV(listData,names,pathToSave,"TimeVsMeasureVsMotion")
                    
            return

################################################
#########         Utility functions
################################################


            
    def save2DSliceToCSV(self,data,filename):
        '''Save a 2D numpy array to CSV file. 
        
        :param data: 2D numpy array
        :param filename: Name of the output file.
        
        If the number of columns is greater than 255 and and the number of rows is smaller than 255, we trasnpose the rows and columns\
        in order to be able to open the file in Excel or Numbers. If both dimensions are greater than 255, we remove the same number of columns on the left \
        and on the right, such that the remaining number of columns is 255.
        
        The csv file is stored in the CSV_Data/ folder of the simulation directory.
        
        '''
        
        test=np.zeros((2,2))
        typetest= type(test)
        if type(data) != typetest:
            print "type data: %s"%str(type(data))
            raise ValueError('In save2DSliceToCSV, data is not *NumPy* array')
        if len(np.shape(data)) != 2:
            print len(np.shape(data))
            print "shape data: %s"%str(data.shape)
            raise ValueError('In save2DSliceToCSV, data is not NumPy *2D array*')
        
        pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation + "CSV_Data/"
        
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)     
            
        if len(filename) > 4:    
            if filename.find(".csv",len(filename)-1-4) == -1:
                filename = filename + ".csv"
        else:
                filename = filename + ".csv"
        if data.shape[1] >= 255 and data.shape[0] < 255:
            newData = np.flipud(data.transpose((1,0)))
        elif data.shape[0] >= 255 and data.shape[1] < 255:
            newData = data
        elif data.shape[0] >= 255 and data.shape[1] >= 255:
            nbColsToRemove = (data.shape[1] - 255)
            if nbColsToRemove%2 == 0:
                nbLeft = nbColsToRemove / 2
                nbRight = nbLeft
            else:
                nbLeft = int(nbColsToRemove / 2) + 1
                nbRight = int(nbColsToRemove / 2)
            newData = np.zeros((data.shape[0],255))
            newData[:,:] = data[:,nbLeft:data.shape[1] - nbRight ]
        elif data.shape[0] <= 255 and data.shape[1] <= 255:
            newData = data
        else:
            print "Shape data csv, unexpected, : %s"%str(data.shape)

        np.savetxt( pathToSave+filename,newData,delimiter=',')
        
    def save2DListToCSV(self,data,filename):
        '''Save a 2D list into a CSV file.

        :param data: 2D list
        :param filename: Name of the output file.

         The csv file is stored in the CSV_Data/ folder of the simulation directory.
        
        '''
        pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation + "CSV_Data/"
        
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)             
        if len(filename) > 4:    
            if filename.find(".csv",len(filename)-1-4) == -1:
                filename = filename + ".csv"
        else:
                filename = filename + ".csv"
        myFile = open(pathToSave+filename, 'wb')
        wr = csv.writer(myFile, quoting=csv.QUOTE_ALL)
        for row in data:
            wr.writerow(row)
        myFile.close()
        



def exportListDataToCSV(listsValues,names,pathToSave,filename):
    '''Export a list to a CSV file.
    
    :param listsValues: List of values to store.
    :param names: List of headers.
    :param pathToSave: Path where to save the CSV file.
    :param filename: Name of the output file. 
    
    '''
   
    dataToExport = list()
   
    for idx, values in enumerate(listsValues):
        listVal = list(values)
        listVal.insert(0,names[idx])
        dataToExport.append(listVal)

    
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)
        
    if not filename.endswith(".csv"):
        filename = filename + ".csv"
        
    np.savetxt(pathToSave +filename ,np.array( dataToExport),fmt='%s', delimiter=',',newline='\r\n')
    print "Data exported: %s"%(pathToSave +filename)


if __name__ == '__main__':
    print "Simulator PID: %s on %s"%(str( os.getpid()) , socket.gethostname())
    if len(sys.argv) == 6 or len(sys.argv) == 5 or len(sys.argv) == 4:
        configFile = None
        for item in sys.argv[1:]:
            if item.endswith(".py") or item.endswith(".txt"):
                configFile = item
                globalVar = GlobalVariables(item)
                break
            elif item == '--motion':
                print item
                with open('tmp_ConfigData.txt', 'w') as myFile:
                    myFile.write("staticDelivery = False \n")
#                     myFile.write("arrayAmpl =  [0.0,  0,  0.0] \n")
                configFile = 'tmp_ConfigData.txt'
                globalVar = GlobalVariables(configFile)
                os.system("%s tmp_ConfigData.txt"%rmCmd)
                break
        if configFile is None:
            globalVar = GlobalVariables(None)
        
        pathToTest = globalVar.mainPath + globalVar.nameSimulation
        print pathToTest
        if not os.path.exists(pathToTest):
                os.makedirs(pathToTest)
        pathToSave = pathToTest +  globalVar.nameProfilage

        if globalVar.stdOutErrToFile:
            bkpStdout = sys.stdout
            bkpStderr = sys.stderr    
            fStdout = open(pathToTest + "stdout.txt","w")
            fStdErr = open(pathToTest + "stderr.txt","w")
            sys.stdout = fStdout
            sys.stderr = fStdErr
        print "Simulator PID: %s"%str( os.getpid())
        simulStartDate = datetime.datetime.now()
        print "Simulation began on %s"%simulStartDate.ctime()
        try:
            cProfile.run('MSPT(sys.argv[1:], globalVar)',pathToSave)
            
            simulEndDate = datetime.datetime.now()
            timeSimul = simulEndDate - simulStartDate
            print "Simulation ended on %s in %s"%(simulEndDate.ctime(),str(timeSimul))
            
            #Copy config file
#             if configFile is not None:
#                 os.system("%s %s %s"%(cpCmd,configFile, pathToTest))
#             else:
            with open(pathToTest + 'MSPT_%s_ConfigData.txt'%globalVar.nameSimulation[:-1], 'w') as myFile:
                myFile.write("#############################\n")
                myFile.write("## Configuration file for simulation %s \n"%(globalVar.nameSimulation[:-1]))
                myFile.write("#############################\n")
                for key in sorted(globalVar.keys()):
                    if key != 'motionInfo' and key != 'settingsScanning':
                        if key in globalVar.authorizedKeys:
                            if not isinstance(globalVar[key],str):
                                myFile.write("%s = %s \n"%(key,str(globalVar[key])))
                            else:
                                myFile.write("%s = '%s' \n"%(key,globalVar[key]))
                    else :
                        for subkey in  sorted(globalVar[key].keys()):
                            if key in globalVar.authorizedKeys:
                                if not isinstance(globalVar[key][subkey],str):
                                    myFile.write("%s = %s \n"%(subkey,str(globalVar[key][subkey])))
                                else:
                                    myFile.write("%s = '%s' \n"%(subkey,globalVar[key][subkey]))

            
            if os.path.exists(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]):
                sys.stdout.flush()
                os.system("%s %s %s"%(cpCmd,globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1], globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt"))
            
            textEmail = "<p>Simulation %s started on %s (server: %s), ended properly on %s in %s.</p><p>Good job!</p>"%( Email.htmlTextBoldRed(globalVar.nameSimulation),Email.htmlTextBoldRed(simulStartDate.ctime()),Email.htmlTextBoldRed(str(socket.gethostname())),simulEndDate.ctime(),str(timeSimul))
            attachment = [globalVar.mainPath + globalVar.nameSimulation + "DVH_Data/"+"DVH_SimulatedDose.png",\
            globalVar.mainPath + globalVar.nameSimulation + "DVH_Data/"+"DVH_ReferenceDose.png"]
            dataToZip = [globalVar.mainPath + globalVar.nameSimulation]
            if os.path.exists(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt"):
                time.sleep(40) # Allow to finish writting in the file
                attachment.append(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt")
                dataToZip.append(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt")
            zipname = globalVar.mainPath + globalVar.nameSimulation[0:-1] + '.zip'
            zipFiles.archiveDataFromListPath(dataToZip , zipname, goToSubDir = True ,listExt = None, rmCommonPath = True)
            if globalVar.copyZip:
                if not os.path.isdir(globalVar.pathToSendZip):
                    os.makedirs(globalVar.pathToSendZip)
                os.system("%s %s %s"%(scpCmd,zipname, globalVar.pathToSendZip))
                textEmail = textEmail + "<p>Zipfile %s sent to %s.</p>"%(zipname, globalVar.pathToSendZip)
            emailSent = True
            try:
                myEmail = emailManager( globalVar.emailUser,globalVar.emailPwd,globalVar.emailSMTP,globalVar.emailSMTPPort)
                myEmail.mail(globalVar.emailRecipient,"Simulation done",textEmail,attach=attachment)
            except:
                attachment = [globalVar.mainPath + globalVar.nameSimulation + "DVH_Data/"+"DVH_SimulatedDose.png",\
            globalVar.mainPath + globalVar.nameSimulation + "DVH_Data/"+"DVH_ReferenceDose.png"]
                textEmail = textEmail + "<p>Only DVHs attached.</p>"
                try:
                    myEmail = emailManager( globalVar.emailUser,globalVar.emailPwd,globalVar.emailSMTP,globalVar.emailSMTPPort)
                    myEmail.mail(globalVar.emailRecipient,"Simulation done",textEmail,attach=attachment)
                except:
                    print "Sorry ... something went wrong. The final email has not been sent. Check emailUser, emailPwd, emailSMTP,emailSMTPPort or emailRecipient in the configuration file"
                    emailSent = False
            if emailSent:
                print "Email sent to %s"%globalVar.emailRecipient
            if os.path.exists(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt"):
                os.remove(globalVar.mainPath+"output"+globalVar.nameSimulation[3:-1]+".txt")
        except:
            simulEndDate = datetime.datetime.now()
            print "Simulation ended on %s"%simulEndDate.ctime()
            textEmail = "<p>The simulation %s (PID: %s - server: %s) started on %s crashed...</p><p>An unexpected error happened on %s : </p><p>----------------</p><p>%s  </p><p>----------------</p>"%(Email.htmlTextBoldRed(globalVar.nameSimulation),str( os.getpid()),Email.htmlTextBoldRed(str(socket.gethostname())),Email.htmlTextBoldRed(simulStartDate.ctime()),simulEndDate.ctime(),Email.htmlTextBoldRed(str(sys.exc_info()[0]).replace('<',"").replace('>',"")))
            strInfo = str(traceback.format_exc())
            newStrInfo =''
            for line in strInfo.split('\n'):
                newStrInfo =newStrInfo +'<br>'+line+'</br>'
            textEmail = textEmail +"<p>"+ newStrInfo + "</p><p>----------------</p><p>Try again...!</p>"
            print "Unexpected error:\n", sys.exc_info()[0]
            
            try:
                myEmail = emailManager( globalVar.emailUser,globalVar.emailPwd,globalVar.emailSMTP,globalVar.emailSMTPPort)
                myEmail.mail(globalVar.emailRecipient,"Simulation crashed",textEmail,attach=None)
                print "Crash report sent to %s"%globalVar.emailRecipient
            except:
                print "Sorry ... something went wrong. The crash report email has not been sent. Check emailUser, emailPwd, emailSMTP,emailSMTPPort or emailRecipient in the configuration file"
            raise
        
        finally:
            if globalVar.removeProfilageData:
                os.remove(pathToSave)
            if globalVar.stdOutErrToFile:
                sys.stdout = bkpStdout
                sys.stderr = bkpStderr  
            
                fStdout.close()
                fStdErr.close()
    else:
        print "MSPT terminated because it doesn't have the right number of input parameters. Expected: 3 (CT,RP,RS), 4 (CT,RP,RS,RD/config/--motion) or 5(CT,RP,RS,RD,config/--motion), given: %i"%(len(sys.argv)-1)
