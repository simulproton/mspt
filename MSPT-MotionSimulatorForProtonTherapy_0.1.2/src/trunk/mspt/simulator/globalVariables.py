########################################################################
#
# globalVariables.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# July 2013 - Modified Jan 13 2014
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
########################################################################
import numpy as np
import traceback


listAuthorizedKeys = ['MeV', 'typeFloat', 'removeProfilageData', 'stdOutErrToFile', 'staticDelivery', 'ddCurves', 'displayScanPathUpdates', \
'mainPath', 'zipSavedData', 'copyZip', 'pathToCopyZip', 'listSaveCompDoseAlong', 'saveRefDose', 'listSaveRefDoseAlong' , 'skipSaveDeff',\
'listSaveDeffAlong', 'saveDeffCPKL', 'saveCTImage', 'listSaveCTImageAlong', 'skipSaveDensityImage', 'listSaveDensImAlong', 'saveDensCPKL', \
'skipSaveStpPwrImage', 'listSaveStPwrImAlong', 'saveStpPwrCPKL', 'skipSaveBeamlets', 'saveBeamletsCPKL', 'displaySliceOrderingInfo',\
'exportSliceOrderingInfo', 'saveCompDoseAsDcm', 'nameCompDose', 'saveRefDoseAsDcm', 'nameRefDose', 'storeMasksAsCPKL' , 'emailUser', \
'emailPwd', 'emailSMTP', 'emailSMTPPort', 'emailRecipient', 'protPerMU', 'dvhRelativeVolume', 'importNewMassStopPwr','nameStopPwrCorrecFactor', 'nameSimulation',\
'nameDoseCorrecFactor','nameCTToDensTable','nameProfilage', 'sliceRepaint', 'volumeRepaint', 'ctPadding', 'padValueFRC', 'paddedCTValue', 'storeInterDoseDens', 'motionDuringBeam', \
'compensationDynamic', 'unlimitedRescanning', 'repaintingFactor', 'sliceRepaint', 'volumeRepaint', 'findStartPos', 'NameNewRPFile', \
'measurementStdev', 'measurementAvg', 'addMargins', 'marginsParam', 'updateMargins', 'evaluateDiffStatDyna', 'arrayPeriod', \
'arrayAmpl', 'arrayPhase', 'typeMotion','fileRecordedMotion', 'breathingPeriod', 'magnitude', 'initialPhase', 'arrayStDevPeriod', 'arrayStDevAmpl', \
'distributionMag', 'distributionPeriod', 'variationsMag', 'variationsPeriod', 'motion3D', 'spotSettlingTime', 'synchrotron' ,\
'timeEnergyChange', 'lateralScanningSpeed', 'verticalScanningSpeed', 'spillDuration', 'speedGantryRotation', 'beamIntensity', \
'timeMeasureDisplacement', 'unitTimeForDelivery' , "TypeScanning", "Repainting", "3DMargins"]


class GlobalVariables(dict):

    '''The class GlobalVariables is a subclass of dictionary to manage the configuration parameters from the input configuration file provided.
    
     :param configFilename: Path to the configuration file. (.txt or .py)
     
    '''
    
    
    def __init__(self, configFilename):
        super(GlobalVariables, self).__init__()

        self.__dict__ = self
        if configFilename is not None:
            try:
                execfile(configFilename, self)
            except:
                print "Error when reading the configuration file."
                print "Correct it based on the following information:"
                print traceback.format_exc()
                
        print ">>>>>>>>>> Configuring MSPT.... <<<<<<<<<<<"
        self.checkConfigData()
        if self.staticDelivery:
            self.processConfigStaticSimul()
        else:
            self. processConfigMotionSimul()
        
        self.buildSettingsForScanningDelivery()
        
        print ">>>>>>>>>> Configuration done! <<<<<<<<<<<"
        
    @property
    def authorizedKeys(self):
        '''Get a list of the authorized keys in the configuration file.
        
        '''
        return listAuthorizedKeys
    
    def checkConfigData(self):
        '''Check the configuration keys to ensure that all the needed keys are present. If they are missing they are added \
        to the GlobalVariables instance with the default settings. The parameters verified in this function correspond to \
        variables common between static and dynamic delivery.
        
        '''
    
        self['MeV'] = 1.602e-13 #1MeV = 1.602e-13J
    
        if 'typeFloat' not in self:
            self['typeFloat'] = 'float64'
            print "'typeFloat' set to default:  'float64'"
            print "\t\t-------------------"
            
        if 'removeProfilageData' not in self:
            self['removeProfilageData'] = True
            print "'removeProfilageData' set to default:  True"
            print "\t\t-------------------"
            
        if 'stdOutErrToFile' not in self:
            self['stdOutErrToFile'] = False
            print "'stdOutErrToFile' set to default:  False"
            print "\t\t-------------------"
        
        if 'staticDelivery' not in self:
            self['staticDelivery'] = True
            print "'staticDelivery' set to default:  True"
            print "\t\t-------------------"
            
        if 'ddCurves' not in self:
            self['ddCurves'] = 'RayStation'
            print "'staticDelivery' set to default:  'RayStation'"
            print "\t\t-------------------"
            
        if 'displayScanPathUpdates' not in self:
            self['displayScanPathUpdates'] = False
            print "'displayScanPathUpdates' set to default:  False"
            print "\t\t-------------------"
            
        if 'mainPath' not in self:
            self['mainPath'] = '../Data_Simulation/'
            print "'mainPath' set to default:  '../Data_Simulation/'"
            print "\t\t-------------------"
             
        if 'zipSavedData' not in self:
            self['zipSavedData'] = True
            print "'zipSavedData' set to default:  True"
            print "\t\t-------------------"
            
        if 'copyZip' not in self:
            self['copyZip'] = False
            print "'copyZip' set to default:  False"
            print "\t\t-------------------"
            
        if self.copyZip:
            if 'pathToCopyZip' not in self:
                self['pathToCopyZip'] = '../Data_Zip/'
                print "'pathToCopyZip' set to default:  '../Data_Zip/'"
                print "\t\t-------------------"        
                
                
        if 'listSaveCompDoseAlong' not in self:
            self['listSaveCompDoseAlong'] = ['z']
            print "'listSaveCompDoseAlong' set to default:  ['z']"
            print "\t\t-------------------"
            
        if 'saveRefDose' not in self:
            self['saveRefDose'] = True
            print "'saveRefDose' set to default:  True"
            print "\t\t-------------------"
            
      
        if self.saveRefDose and 'listSaveRefDoseAlong' not in self:
            self['listSaveRefDoseAlong'] = ['z']
            print "'listSaveRefDoseAlong' set to default:  ['z']"
            print "\t\t-------------------"        
        
    
        if 'skipSaveDeff' not in self:
            self['skipSaveDeff'] = True
            print "'skipSaveDeff' set to default:  True"
            print "\t\t-------------------"    
    
        if not self.skipSaveDeff and 'listSaveDeffAlong' not in self:
            self['listSaveDeffAlong'] = ['z']
            print "'listSaveDeffAlong' set to default:  ['z']"
            print "\t\t-------------------"   
            
        if not self.skipSaveDeff and 'saveDeffCPKL' not in self:
            self['saveDeffCPKL'] = False
            print "'saveDeffCPKL' set to default:  True"
            print "\t\t-------------------"   
            
        if 'saveCTImage' not in self:
            self['saveCTImage'] = True
            print "'saveCTImage' set to default:  True"
            print "\t\t-------------------"    
            
        if  self.saveCTImage and 'listSaveCTImageAlong' not in self:
            self['listSaveCTImageAlong'] = ['z']
            print "'listSaveCTImageAlong' set to default:  ['z']"
            print "\t\t-------------------"   
                    
        if 'skipSaveDensityImage' not in self:
            self['skipSaveDensityImage'] = True
            print "'skipSaveDensityImage' set to default:  True"
            print "\t\t-------------------"    
           
        if not self.skipSaveDensityImage and 'listSaveDensImAlong' not in self:
            self['listSaveDensImAlong'] = ['z']
            print "'listSaveDensImAlong' set to default:  ['z']"
            print "\t\t-------------------"  
        
        if not self.skipSaveDensityImage and 'saveDensCPKL' not in self:
            self['saveDensCPKL'] = False
            print "'saveDensCPKL' set to default:  False"
            print "\t\t-------------------"   
           
        if 'skipSaveStpPwrImage' not in self:
            self['skipSaveStpPwrImage'] = True
            print "'skipSaveStpPwrImage' set to default:  True"
            print "\t\t-------------------"    
           
        if not self.skipSaveStpPwrImage and 'listSaveStPwrImAlong' not in self:
            self['listSaveStPwrImAlong'] = ['z']
            print "'listSaveStPwrImAlong' set to default:  ['z']"
            print "\t\t-------------------"  
        
        if not self.skipSaveStpPwrImage and 'saveStpPwrCPKL' not in self:
            self['saveStpPwrCPKL'] = False
            print "'saveStpPwrCPKL' set to default:  False"
            print "\t\t-------------------"   
            
        if 'skipSaveBeamlets' not in self:
            self['skipSaveBeamlets'] = True
            print "'skipSaveBeamlets' set to default:  True"
            print "\t\t-------------------"    
 
        if not self.skipSaveBeamlets and 'saveBeamletsCPKL' not in self:
            self['saveBeamletsCPKL'] = False
            print "'saveBeamletsCPKL' set to default:  False"
            print "\t\t-------------------"   

        
        if 'displaySliceOrderingInfo' not in self:
            self['displaySliceOrderingInfo'] = False
            print "'displaySliceOrderingInfo' set to default:  False"
            print "\t\t-------------------"             
        
        if 'exportSliceOrderingInfo' not in self:
            self['exportSliceOrderingInfo'] = False
            print "'exportSliceOrderingInfo' set to default:  False"
            print "\t\t-------------------"              
        
        if 'saveCompDoseAsDcm' not in self:
            self['saveCompDoseAsDcm'] = True
            print "'saveCompDoseAsDcm' set to default:  True"
            print "\t\t-------------------"              
       
        if self.saveCompDoseAsDcm and 'nameCompDose' not in self:
            self['nameCompDose'] = 'RD_CompDose_MSPT.dcm'
            print "'nameCompDose' set to default:  'RD_CompDose_MSPT.dcm'"
            print "\t\t-------------------"              


        if 'saveRefDoseAsDcm' not in self:
            self['saveRefDoseAsDcm'] = True
            print "'saveRefDoseAsDcm' set to default:  True"
            print "\t\t-------------------"              
       
        if self.saveRefDoseAsDcm and 'nameRefDose' not in self:
            self['nameRefDose'] = 'RD_RefDose.dcm'
            print "'nameRefDose' set to default:  'RD_RefDose.dcm'"
            print "\t\t-------------------"              

        if 'storeMasksAsCPKL' not in self:
            self['storeMasksAsCPKL'] = False
            print "'storeMasksAsCPKL' set to default:  False"
            print "\t\t-------------------"              


        if 'emailUser' not in self:
            self['emailUser'] = None
            print "'emailUser' set to default:  None"
            print "\t\t-------------------"              


        if 'emailPwd' not in self:
            self['emailPwd'] = None
            print "'emailPwd' set to default:  None"
            print "\t\t-------------------"              

        if 'emailSMTP' not in self:
            self['emailSMTP'] = None
            print "'emailSMTP' set to default:  None"
            print "\t\t-------------------"              
        
        if 'emailSMTPPort' not in self:
            self['emailSMTPPort'] = None
            print "'emailSMTPPort' set to default:  None"
            print "\t\t-------------------"              

        if 'emailRecipient' not in self:
            self['emailRecipient'] = None
            print "'emailRecipient' set to default:  None"
            print "\t\t-------------------"              
    
        if 'protPerMU' not in self:
            self['protPerMU']= 1e9
            print "'protPerMU' set to default:  1e9"
            print "\t\t-------------------"              

        
        if 'dvhRelativeVolume' not in self:
            self['dvhRelativeVolume'] = True
            print "'dvhRelativeVolume' set to default:  True"
            print "\t\t-------------------"              

        
        if 'importNewMassStopPwr' not in self:
            self['importNewMassStopPwr'] = False
            print "'importNewMassStopPwr' set to default:  False"
            print "\t\t-------------------"              
            
        if 'nameStopPwrCorrecFactor' not in self:
            self['nameStopPwrCorrecFactor'] = 'MSPT'
            print "'nameStopPwrCorrecFactor' set to default:  'MSPT'"
            print "\t\t-------------------"                        
            
        if 'nameDoseCorrecFactor' not in self:
            self['nameDoseCorrecFactor'] = 'MSPT'
            print "'nameDoseCorrecFactor' set to default:  'MSPT'"
            print "\t\t-------------------"                        


        if 'nameCTToDensTable' not in self:
            self['nameCTToDensTable'] = 'MSPT'
            print "'nameCTToDensTable' set to default:  'MSPT'"
            print "\t\t-------------------"                        
    
    
    def processConfigStaticSimul(self):
        '''Process the keys:values specifically needed for the static delivery.
        
        '''
        
        if 'nameSimulation' not in self:
            self['nameSimulation'] = 'runSimulStatic/'
            print "'nameSimulation' set to default:  'runSimulStatic/'"
            print "\t\t-------------------"
        else:
            if self.nameSimulation[-1] != "/":
                self.nameSimulation = self.nameSimulation + '/'
            if self.nameSimulation[0:3] != 'run':
                self.nameSimulation = 'run'+self.nameSimulation

        if 'nameProfilage' not in self:
            self['nameProfilage'] = 'profilagerunSimulStatic.profile'
            print "'nameProfilage' set to default:  'profilagerunSimulStatic.profile'"
            print "\t\t-------------------"              
    
        self['sliceRepaint'] = 'FullWeightRepainting'
        self['volumeRepaint'] = 'NonVolumetric'

    
    def processConfigMotionSimul(self):
        '''Process the keys:values specifically needed for the dynamic delivery.
        
        '''
        if 'nameSimulation' not in self:
            self['nameSimulation'] = 'runSimulDynamic/'
            print "'nameSimulation' set to default:  'runSimulDynamic/'"
            print "\t\t-------------------"     
        else:
            if self.nameSimulation[-1] != "/":
                self.nameSimulation = self.nameSimulation + '/'
            if self.nameSimulation[0:3] != 'run':
                self.nameSimulation = 'run'+self.nameSimulation

        if 'nameProfilage' not in self:
            self['nameProfilage'] = 'profilagerunSimulDynamic.profile'
            print "'nameProfilage' set to default:  'profilagerunSimulDynamic.profile'"
            print "\t\t-------------------"              
            
        if 'ctPadding' not in self:
            self['ctPadding'] = True
            print "'ctPadding' set to default:  True"
            print "\t\t-------------------"              
            
        if 'padValueFRC' not in self:
            self['padValueFRC'] = [10,0,10]
            print "'padValueFRC' set to default:  [10,0,10]"
            print "\t\t-------------------"              
        
        if 'paddedCTValue' not in self:
            self['paddedCTValue'] = -1000
            print "'paddedCTValue' set to default:  -1000"
            print "\t\t-------------------"              
        
        
        if 'storeInterDoseDens' not in self:
            self['storeInterDoseDens'] = False
            print "'storeInterDoseDens' set to default:  False"
            print "\t\t-------------------"              
            
 
        if 'motionDuringBeam' not in self:
            self['motionDuringBeam'] = True
            print "'motionDuringBeam' set to default:  True"
            print "\t\t-------------------"              
 
        
        if 'compensationDynamic' not in self:
            self['compensationDynamic'] = False
            print "'compensationDynamic' set to default:  False"
            print "\t\t-------------------"              
           
        
        if 'unlimitedRescanning' not in self:
            self['unlimitedRescanning'] = False
            print "'unlimitedRescanning' set to default:  False"
            print "\t\t-------------------"              
        
        if not self.unlimitedRescanning:
            if 'repaintingFactor' not in self:
                self['repaintingFactor'] = 1
                print "'repaintingFactor' set to default:  1"
                print "\t\t-------------------"              
        
        if 'sliceRepaint' not in self:
            self['sliceRepaint'] = 'ScaledRepaint'#['IsolayerRepaint','ScaledRepaint', 'FullWeightRepainting']
                                                #Note: 'FullWeightRepainting': we deliver full weight at once instead of a maximum weight or a scaled weight
            print "'sliceRepaint' set to default:  'ScaledRepaint'"
            print "\t\t-------------------"              
 
        if 'volumeRepaint' not in self:
            self['volumeRepaint'] = 'NonVolumetric' #['Volumetric','NonVolumetric']
            print "' volumeRepaint' set to default:  'NonVolumetric'"
            print "\t\t-------------------"              

        if self.compensationDynamic:
            if 'findStartPos' not in self:
                self['findStartPos'] = False
                print "'findStartPos' set to default:  False"
                print "\t\t-------------------"              
        else:
            self['findStartPos'] = False
        
        
        if self.compensationDynamic:
            if 'NameNewRPFile' not in self:
                self['NameNewRPFile'] = 'RPFile_Compensation.dcm'
                print "'NameNewRPFile' set to default:  'RPFile_Compensation.dcm'"
                print "\t\t-------------------"              

            if 'measurementStdev' not in self:
                self['measurementStdev'] = [ 0.1, 0, 0.1] # Standard deviation of the normal distributions modeling noise in the motion monitor
                                                            #It is defined along [X,Y,Z] IEC Fixed in cm
                print "'measurementStdev' set to default:  [ 0.1, 0, 0.1 ]"
                print "\t\t-------------------"              

            if 'measurementAvg' not in self:
                self['measurementAvg'] = [ 0.1,0 , 0.1 ]#It is defined along [X,Y,Z] IEC Fixed in cm
                print "'measurementAvg' set to default:  [ 0.1, 0, 0.1 ]"
                print "\t\t-------------------"              
           
            if 'addMargins' not in self:
                self['addMargins'] = False
                print "'addMargins' set to default:  False"
                print "\t\t-------------------"              
            
            if self.addMargins:
                if 'marginsParam' not in self:
                    self['marginsParam'] = [10,0,10] # extension in mm for [x , y, z] -> this means for ex. add a margin of (x mm) at the right of the map and (x mm) at the left of the map.
                    print "'marginsParam' set to default:  [10,0,10]"
                    print "\t\t-------------------"              
                
                if 'updateMargins' not in self:
                    self['updateMargins'] = True
                    print "'updateMargins' set to default:  True"
                    print "\t\t-------------------"              
       
        if 'evaluateDiffStatDyna' not in self:
            self['evaluateDiffStatDyna'] = False
            print "'evaluateDiffStatDyna' set to default:  False"
            print "\t\t-------------------"              
        
        
                #Motion parameters
                #Motion. Array[0]: X direction,Array[1]: Y direction,Array[2]: Z direction. 
                #Directions refer to axis directions in the IEC fixed coordinate system : X_IEC fixed= X_Dicom, Y_IEC fixed = Z_Dicom, Z_IEC fixed = -Y_Dicom 
        if 'arrayPeriod' not in self:
            self['arrayPeriod'] = [4.0,  0.0,   4.0] #Periods in sec 
            print "'arrayPeriod' set to default:  [4.0,  0.0,   4.0]"
            print "\t\t-------------------"              
  
        if 'arrayAmpl' not in self:
            self['arrayAmpl'] = [1.0,  0.0,  1.0] #Ampl in cm   
            print "'arrayAmpl' set to default:  [1.0,  0.0,  1.0]"
            print "\t\t-------------------"              


        if 'arrayPhase' not in self:
            self['arrayPhase'] = [0,   0,   0] #Phase in rad  
            print "'arrayPhase' set to default:  [0,   0,   0]"
            print "\t\t-------------------"              
  

            #Fill motion Info
        motionInfo = dict()
        if 'typeMotion' not in self:
            motionInfo['typeMotion'] = 'motionCos'  # ['motionNoisyCos','motionCos','motionBreathHiccup','motionBreathCough']
            print "'typeMotion' set to default:  'motionCos'"
            print "\t\t-------------------"              
        else:
            motionInfo['typeMotion'] = self.typeMotion

        
        if  motionInfo["typeMotion"] == 'motionRecorded':
        
            if 'fileRecordedMotion' not in self:
                print "'fileRecordedMotion' not in the configuration data. The motion will be set to motionCos."
                motionInfo['typeMotion'] = 'motionCos'  # ['motionNoisyCos','motionCos','motionBreathHiccup','motionBreathCough']
                print "'typeMotion' set to default:  'motionCos'"
                print "\t\t-------------------"  
            else:
                print "motionRecorded: file used %s"%str(self.fileRecordedMotion)
        
        
        if  motionInfo["typeMotion"] in ['motionCos',"motionNoisyCos",'motionBreathHiccup','motionBreathCough']:
            motionInfo['breathingPeriod'] = np.array(self.arrayPeriod, dtype = self.typeFloat)
            motionInfo['magnitude'] = np.array(self.arrayAmpl, dtype = self.typeFloat)
            motionInfo['initialPhase'] = np.array(self.arrayPhase, dtype = self.typeFloat)
            
        else: 
            print " 'typeMotion' is set to %s"%self.typeMotion
            print "\t\t-------------------"          

        if  motionInfo["typeMotion"] in ["motionNoisyCos",'motionBreathHiccup','motionBreathCough']:
            if 'arrayStDevPeriod' not in self:
                self['arrayStDevPeriod'] = [0,   0,   0]  # sigma Periods in sec
                print "'arrayStDevPeriod' set to default:  [0,   0,   0]"
                print "\t\t-------------------"              
            self.arrayStDevPeriod = np.array(self.arrayStDevPeriod, dtype = self.typeFloat)

            if 'arrayStDevAmpl' not in self:
                self['arrayStDevAmpl'] = [0,   0,   0]  #sigma  Ampl in cm 
                print "'arrayStDevAmpl' set to default:  [0,   0,   0]"
                print "\t\t-------------------"              
            self.arrayStDevAmpl  = np.array(self.arrayStDevAmpl , dtype = self.typeFloat)

            motionInfo['distributionMag'] = 'normal'  # ['normal' , 'lognormal']
            motionInfo['distributionPeriod'] = 'normal' # ['normal' , 'lognormal']
            motionInfo['variationsMag'] = self.arrayStDevAmpl
            motionInfo['variationsPeriod'] = self.arrayStDevPeriod
        
        self['motionInfo'] = motionInfo
        
        if 'arrayAmpl' in self:
            if self.arrayAmpl[0] == 0 or self.arrayAmpl[1] == 0 or self.arrayAmpl[2] == 0:
                self['motion3D'] = False
            else:
                self['motion3D'] = True
                

                #Machine delivery characteristics

        if 'spotSettlingTime' not in self:
            self['spotSettlingTime'] = 5e-3 # seconds :  1ms,5ms,10ms / From Dowdell et al 2013
            print "'spotSettlingTime' set to default:  5e-3 s"
            print "\t\t-------------------"              
       
        if 'synchrotron' not in self:
            self['synchrotron'] = True 
            print "'synchrotron' set to default:  True"
            print "\t\t-------------------"              
    
        if 'timeEnergyChange' not in self:
            if self.synchrotron:
                self['timeEnergyChange'] = 2.1 # s   : from Smith et al. 2009 : time to fill synchrotron included
                print "'timeEnergyChange' synchrotron set to default:  2.1 s"
                print "\t\t-------------------"              
            else:
                self['timeEnergyChange'] = 50e-3 # s in Schippers et al 2007  / 150e-3 # s in Pedronni 2004
                print "'timeEnergyChange' cyclotron set to default:  50e-3 s"
                print "\t\t-------------------"              
                
        if 'lateralScanningSpeed' not in self:
            self['lateralScanningSpeed'] = 1000 # cm.s^-1 : from Smith et al. 2009,Gillin 2010 Kraus et al. 2011, Kang 2007
            print "'lateralScanningSpeed' set to default:  1000 cm/s"
            print "\t\t-------------------"              
            
        if 'verticalScanningSpeed' not in self:
            self['verticalScanningSpeed'] = 1000 # cm.s^-1 : Kraus et al. 2011, Kang 2007
            print "'verticalScanningSpeed' set to default:  1000 cm/s"
            print "\t\t-------------------"              

        if 'spillDuration' not in self:
            self['spillDuration'] = 5 #s : Kraus et al. 2011 // 4.4 s in Smith et al. 2009
            print "'spillDuration' set to default:  5s"
            print "\t\t-------------------"              
        
        if 'speedGantryRotation' not in self:
            self['speedGantryRotation'] = 360  # deg/min
            print "'speedGantryRotation' set to default:  360 deg/min"
            print "\t\t-------------------"              
        
        if 'beamIntensity' not in self:
            self['beamIntensity'] = 6.8e10 # protons.s^ -1 : no information regarding this, but this seems more realistic than 4e7 knowing that we use 9e8 protons/MU
                                    #Beam intensity (extraction rate): from test on plan Full CT: average weight: 0.073186 MU. If we assume an average time spent per spot of 10ms (found in paper Grozinger 2008:"The irradiation time for an individual raster point is typically in the order of 5-10 ms."), then average intensity is 0.073186 * 9.303078e8  / 10e-3 = 6.8e10
            print "'beamIntensity' set to default:  6.8e10  protons/s"
            print "\t\t-------------------"              
        
 
        if 'timeMeasureDisplacement' not in self:
            self['timeMeasureDisplacement'] = 1e-3 #sec : time needed to measure the displacement - in Spadea 2011 the sampling rate is 30Hz (1 sampling every 33ms), therfore we assume the time to measure the displacement negligeable compared to 30ms -> 1ms
            print "'timeMeasureDisplacement' set to default:  1e-3 s"
            print "\t\t-------------------"              
 
        if 'unitTimeForDelivery' not in self:
            self['unitTimeForDelivery'] = 1e-3 #sec: unit time used during the delivery
            print "'unitTimeForDelivery' set to default:  1e-3 s"
            print "\t\t-------------------"              
 
 
 
 
 
 

    def buildSettingsForScanningDelivery(self):
        '''Builds a dictionary for the initialization of the ScanningPath (in dicomReader package) .
        
        '''
        settings = dict()
        if self.staticDelivery:
            settings["TypeScanning"] = 'Regular' # ['Random', 'Regular']
            settings["Repainting"] = False
        else:
            settings["TypeScanning"] = 'Regular' # ['Random', 'Regular']
            settings["Repainting"] = True
            settings["SliceRepaint"] = self.sliceRepaint #['IsolayerRepaint','ScaledRepaint', 'FullWeightRepainting']#Note: 'FullWeightRepainting': used for compensation: we deliver full weight instead of o maximum weight or a scaled weight.
            settings["VolumeRepainting"] = self.volumeRepaint #['Volumetric','NonVolumetric']
            settings["UnlimitedRescan"] = self.unlimitedRescanning
            if not settings["UnlimitedRescan"]:
                settings["RepaintFactor"] = self.repaintingFactor
            settings["FindStartPos"] = self.findStartPos
            if self.compensationDynamic:
                settings["AddMapMargins"] = self.addMargins
                if self.addMargins:
                    settings["3DMargins"] = self.motion3D
                    settings["MarginParams"] = self.marginsParam
                    settings["UpdateMargins"] = self.updateMargins
        self['settingsScanning']  = settings        
        
 
 
 
 
    