########################################################################
#
# deliverySimulation_Static_Dynamic.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# June 2013
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

import numpy as np
from ..dataViewer import array2DViewer as viewer2D 
from ..mathTools import tools as mathTools
from ..simulator.tools import Tools as simTools
from ..dicomReader import scanningPath as scanPlan
import sys, os 
import datetime
import cPickle
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
from scipy.interpolate import interp1d
import itertools as iter

#C extension to simulate the beamlet
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBeamletCalc import _beamletSimulation_PythonExt as _beamletSimulation 
from ..extensionsC.PythonExtensionAllPlatforms.PythonExtensionBeamletCalcDouble import _beamletSimulation_PythonExtDouble as _beamletSimulationDouble 

limitSmallWeights = 1e-4
useFastRadiolDepth_DynaDeliv = False 

class simulatorDelivery(object):
    '''Class simulatorDelivery is the core of the simulation. This where the delivery is actually simulated.
    
    :param patient: Patient object (in patient.patient) containing all the 3D arrays for the CT, density, ...
    :param ddCurvesData: DDcurve object (in physics.ddCurveManager) used to store all the depth-dose curves.
    :param configData: dictionary (globalVariables) storing the data configuring the simulator. \
    Mandatory keys:
        
        * *'mainPath'* : corresponds to the path where any simulation will be saved.
        * *'nameSimulation'* : corresponds to the path where the current simulation will be saved.
        * *'typeFloat'* : The type of numpy float to be used: 'float32' or 'float64'
        * *'protPerMU'* : Number of protons per MU to use.
        * *'MeV'*: Constant to convert 1MeV to Joules: 1MeV = 1.602e-13J
        * *'displayScanPathUpdates'* : Display the evolution of the current scanning path of the current energy layer. Scanning Paths are\
        stored in 3D numpy array ( see dicomReader.scanningPath for more information) where each frame is an energy layer. Therefore, if True \
        the frame corresponding the current energy layer will be displayed to show the remaining beam positions.
        * *'skipSaveBeamlets'*: False to save an image of the dose distribution of each beamlet simulated, True otherwise. 
        * *'saveBeamletsCPKL'*: True to save the dose distribution of each beamlet (a 3D numpy array) in a cPickle structure.
        * *'staticDelivery'* : True if the patient is not moving, or False otherwise
        * *'motionDuringBeam'* (if 'staticDelivery' is False): True to allow the patient to move while the beam is irradiating, False otherwise. 
        * *'motion3D'* (if 'staticDelivery' is False): True if the patient motion is in 3D.
        * *'evaluateDiffStatDyna'* (if 'staticDelivery' is False) : True to evaluate the difference for each beam position between the static delivery and the dynamic delivery (i.e., moving patient).\
        False otherwise.
        * *'storeInterDoseDens'* (if 'staticDelivery' is False) : True to store the density image of the moving patient and the dose distribution of the beam when the \
        patient is moving while being irradiated by a pencil beam. False otherwise. 
        * *'listSaveCompDoseAlong'* (if 'staticDelivery' is False):  list of the axis along which the user wants to save the density and dose images of the moving patient. Possible values are\
        'x' (i.e. the images are stored along the column axis of the dicom image), 'y' (i.e. the images are stored along the rows axis of the dicom image), \
        'z' (i.e. the images are stored along the frame axis of the dicom image). Examlpes: listSaveCompDoseAlong = [] , or listSaveCompDoseAlong = ['x'] or listSaveCompDoseAlong = ['x','z'] ... 
        * *'timeEnergyChange'* (if 'staticDelivery' is False) : Time to change the energy of the beam in s.
        * *'lateralScanningSpeed'* (if 'staticDelivery' is False): Beam lateral scanning speed in cm/s. 
        * *'spotSettlingTime'* (if 'staticDelivery' is False): Beam vertical scanning speed in cm/s. 
        * *'speedGantryRotation'* (if 'staticDelivery' is False): Gantry rotational speed in degrees/min.
        * *'beamIntensity'* (if 'staticDelivery' is False): Number of protons per seconds being delivered.
        * *'unitTimeForDelivery'* (if 'staticDelivery' is False): Unit time for motion simulation in seconde. Every unit time the patient position and radiological depth is updated.
        * *'synchrotron'* (if 'staticDelivery' is False): True to simulate a synchrotron ( time energy change, refill the synchrotron). False for a cyclotron.
        * *'spillLength'* (if 'staticDelivery' and synchrotron' are True): Duration of a proton spill in s.
        * *'compensationDynamic'* (if 'staticDelivery' is False): True to allow motion compensation, False otherwise.
        * *'findStartPos'* (if 'compensationDynamic' is True): Find a good starting position in the energy layer( see dicomReader.scanningPath for more information).
        * *'timeMeasureDisplacement'* (if 'compensationDynamic' is True): Time to measure the patient motion is s. 
       
    
    
    '''
    
    
    def __init__(self, patient,ddCurvesData, configData, motion = None ,monitor = None):
        self._globalVar = configData
        self._pathToOutputImages = self._globalVar.mainPath + self._globalVar.nameSimulation 
        ############# If Dynamic delivery ############# 
        if not self._globalVar.staticDelivery:
            if motion == None:
                raise ValueError('In simulator delivery - dynamic mode - no motion has been provided for the patient...')
            else:
                self._motion = motion
                if self._globalVar.compensationDynamic:
                    self._monitoringSystm = monitor
                    if monitor is None:
                        raise ValueError("Monitor is None in delivery with compensation")
                self._time = 0
                self._deffStatic = None # To be used when comparing computed dose for motion to computed dose in static

        ############# End - If Dynamic delivery #############       
        if patient is None or ddCurvesData is None:
            raise ValueError ("Input data (patient or ddCurves) is None")
        self._patient = patient
        self._ddCurvesData = ddCurvesData
        
        
        

    def simulateDelivery(self, spotPositions ): 
        '''Simulates the delivery of a proton therapy treatment plan to a patient based on the simulator configuration. \
        It uses the beam scanning method.
        
        
        :param spotPosition: Treatment plan stored as a scanning path (dicomReader.scanningPath).
        
        :returns: None if the delivery is performed for a static patient or a moving patient without compensation. It returns a compensated treatment plan if the patient \
        is moving and 'compensationDynamic' is True: dictionary[beam angle][Energy Rounded][ "ScanSpotPositionMap" / "ScanSpotMetersetWeights"]

        '''
        self._tool = simTools(self._globalVar)
        self._receivedDose = np.zeros(self._patient.shape(), dtype = self._globalVar.typeFloat, order ='C')
        
        #Get spacing between spots in cm
        self._xSpacing = self._patient.getSpacingInfo('x_spacing') / 10.0   # from mm to cm
        self._ySpacing = self._patient.getSpacingInfo('y_spacing') / 10.0 # from mm to cm
        self._zSpacing = self._patient.getSpacingInfo('z_spacing') / 10.0 # from mm to cm
        ############# If Dynamic delivery #############
        if not self._globalVar.staticDelivery:
            print "Starting dynamic delivery simulation..."
            print "\t------------------"
            print "Motion Information:\n\t%s"%str(self._motion)
            print "Motion allowed during a beamlet delivery: %s"%(self._globalVar.motionDuringBeam)
            print "\t------------------"
            if self._globalVar.compensationDynamic:
                    print self._monitoringSystm
                    print "\t------------------"
            self._staticPatientLoc = self._patient.getPatientExtentInStaticModel()
            self.updateTime('Reset') 
            self.funcGetDataPatient = self._patient.getPatientDataForDisplVector
            #self.computeDeff = self.updateDeff
            self._prevRotationAngle = 0
            if self._globalVar.compensationDynamic:
                nbCompensatedSpots = 0
                nbCorrectSpots = 0
                nbMissedSpots = 0
                dictRepaintStat = dict()
                compensatedPlan = dict()
        ############# End - If Dynamic delivery ############# 
        else: 
            print "Starting static delivery simulation..."
            self.funcGetDataPatient = self._patient.getPatientArrayForAttr
            #self.computeDeff = self.computeDeffStatic
        nbBeams = len(spotPositions)
        t_simul_start = datetime.datetime.now()
        nbTargets = 0
        sumWeights = 0
        
        for idxBeam,dataBeam in enumerate(spotPositions):
            #####
            #dataBeam : (0 : scanningPathSingleBeam ,1: angle , 2: isocenter , 
            #         3: vsad, 4: string 'NewAngle' or 'SameAngle' depending if the angle change between beams) 
            #####              
            t_start_beam = datetime.datetime.now()
            print "Beam %i / %i"%(idxBeam+1,nbBeams)
            ### Get beam angle: In proton therapy, 1 angle/beam 
            angle = dataBeam[1]
            #angle = 45
            print "Beam angle(degrees): %f"%angle
            if not self._globalVar.staticDelivery and self._globalVar.compensationDynamic:
                dictRepaintStat[angle] = dict()
                compensatedPlan[angle] = dict()
            
            ############# If Dynamic delivery #############
            if not self._globalVar.staticDelivery and dataBeam[4] == 'NewAngle' :
                self.updateTime('GantryRotation',info = angle)
            ############# End - If Dynamic delivery ############# 
            
            #Set the patient and the source in the 3D space:
            ###Get isocenter
            isocenterPos = [(value / 10.0) for value in dataBeam[2]] #spotPositions: data provided in mm! it must be converted to cm!
            
            self._patient.setIsocenter(dataBeam[2])
            indIsocenter = self._patient.getIsocenterInd()
            print "Index isocenter: %s"%str(indIsocenter)
            
            
            # Get coordinate system in the patient
            X = np.ascontiguousarray(self.funcGetDataPatient('x'))
            X = X / 10.0 # from mm to cm
            Y = np.ascontiguousarray(self.funcGetDataPatient('y'))
            Y = Y / 10.0 # from mm to cm
            Z = np.ascontiguousarray(self.funcGetDataPatient('z') )
            Z = Z / 10.0 # from mm to cm 
            
            
            
            
            #Convert coordinate system in the IEC Fixed coordinate system
                #Get isocenter
            isocenterPos = [(value / 10.0) for value in dataBeam[2]] #spotPositions: data provided in mm! it muste be converted to cm!

            #Convert X,Y,Z  coordinates from the patient dicom coordinate system to the IEC fixed gantry coordinate system
            (X_iec_f,Y_iec_f,Z_iec_f) = mathTools.fromCTToIECFixedCoord(X,Y,Z,isocenterPos)
            
            #Starting vector for radiological depth computation
            self._startVec = np.array( [X_iec_f[0,-1,0],Y_iec_f[0,-1,0],Z_iec_f[0,-1,0]],dtype = self._globalVar.typeFloat, order = 'C')
            #Incremental vector for radiological depth computation
            self._incVec = np.array([self._xSpacing,self._zSpacing,self._ySpacing],dtype = self._globalVar.typeFloat, order = 'C')
            #First voxel coordinates used in the beamlet computation
            self._x0y0z0 = np.array([X_iec_f[0,0,0],Y_iec_f[0,0,0],Z_iec_f[0,0,0]], dtype = self._globalVar.typeFloat, order = 'C')
            self._IECFSpacing = np.array([ X_iec_f[0,0,1]-X_iec_f[0,0,0],\
                                             Y_iec_f[1,0,0]-Y_iec_f[0,0,0],\
                                              Z_iec_f[0,1,0]-Z_iec_f[0,0,0]], dtype = self._globalVar.typeFloat, order = 'C')# Note: the y IEC F axis is along the frames axis
    
            

    
            sourceLocation = [0,0] # in the IEC fixed coordinate system
            distanceSourceIsocenter = [(value / 10.0) for value in dataBeam[3]][0] #We assume that the beam is composed of 2 gaussian "curves" 1 in X direction, 1 in Y direction. Data provided in  'VirtualSourceAxisDistances' is composed by 2 coordinates (X,Y), we think that correspond to the distance with the virtual source generating the X gaussian and th Y gaussian.
            sourceLocation.append(distanceSourceIsocenter)
            
            print "Source location: %s"%str(sourceLocation)
            
            
            ###Compute rotation (around the patient - around y axis of the IEC fixed coord systm) matrix.
            if dataBeam[4] == 'NewAngle':
                self._rotMatrix = np.array(mathTools.rotationMatrixForAngle(angle),dtype=self._globalVar.typeFloat,order='C')
                nVec = np.zeros(3) #nVec corresponds to the beam trajectory("depth" vector, on-axis): normal vector to the BEV plane. This is along the z axis rotated.
                nVec[0] = self._rotMatrix[0][2] 
                nVec[1] = self._rotMatrix[1][2]
                nVec[2] = self._rotMatrix[2][2]
                
                
                self._sourceVec = (sourceLocation[2] * nVec).astype(self._globalVar.typeFloat)
                #Compute rotated coordinates (cm) 
                (Xr,Yr,Zr) = mathTools.rot3D(X_iec_f,Y_iec_f,Z_iec_f,np.identity(3),self._rotMatrix,self._globalVar.typeFloat)
            nbSpots = 0
            
                        
            listSpotsVisited = None
            listSpotsDelivered = None
            spotsVisitedButWeightNull = 0
            unexplainedWeightNull = 0
            countSmallWeights = 0
            logSpotsSmallWeights = []
            countNotSmallWeights = 0
            sumNotSmallWeights = 0
            
            for idxSpot,dataSpot in enumerate(dataBeam[0]): # dataBeam[0] : class ScanningPathSingleBeam
                #####
                #dataSpot : (0 : weight ,1: index in scanning path structure (3D) , 2: x shift in IEC Gantry , 
                #         3: y shift in IEC Gantry, 4: energy , 5:spot size , 
                #         6:string 'SameEnergy' or 'NewEnergy' depending on what was 
                #         the energy of the previous spot - it will be 'NewEnergy' for 1st spot ) 
                #####  
                ### Get beam size (width) 
                sig0 = [((value/ 2.3548200) / 10.0) for value in dataSpot[5]]
                T0 = dataSpot[4]
                
                
                T0_round = np.floor(T0)
                print "Start new spot #%i"%idxSpot
                print "Energy rounded: %f"%T0_round
                nbSpots = nbSpots + 1
                
                compensated = False
                
                
                if self._globalVar.displayScanPathUpdates:
                    #Display the current scanning path slice
                    currentEnergyLayer = dataBeam[0].getEnergySliceWeights(T0)
                    viewer2D.displayImageFromArray( currentEnergyLayer, "Scan path - energy layer %i MeV"%T0_round)
                    filenameScanPath = "ScanningPath%0.9d_%0.3dMeV"%(idxSpot,T0_round)
                    viewer2D.saveFigure( self._pathToOutputImages + "ScanningPathUpdates/",filenameScanPath , ext = 'png', subplotIdx = 111)
                    viewer2D.closePlot()
            

                if dataSpot[6] == 'NewEnergy':
                    ############# If Dynamic delivery #############
                    if not self._globalVar.staticDelivery:
                        self.updateTime('Energy')
                    ############# End - If Dynamic delivery ############# 
                    self._patient.computeStoppingPowerArray(energy = T0_round)
                    self._protPerMU = self._globalVar.protPerMU
                    dmax = self._ddCurvesData.dataForEnergy('BraggPeak' , T0)
                    print "For energy %0.2f Bragg Peak at:%s cm in water." %(T0,str(dmax))
                    #Depth dose curve data (Gy.cm^2)
                    z_dd = self._ddCurvesData.dataForEnergy('DD_Curve' , T0)[0,:].astype(self._globalVar.typeFloat) # depths at which dose was computed
                    dd =  self._ddCurvesData.dataForEnergy('DD_Curve' , T0)[1,:].astype(self._globalVar.typeFloat)# computed dose for depth dose curve
                    
                    
                    #Get number of protons used for the Monte Carlo simulation providing the depth dose curve
                    nbProtons = self._ddCurvesData.dataForEnergy('NbProtons', T0)

                    if self._globalVar.staticDelivery:
                        self._StPwrRatioGrid =np.ascontiguousarray(self.funcGetDataPatient('stpPwrRatio'))
                        self._DensityGrid = self.funcGetDataPatient('density')
                        self._radiolDepth = self._patient.computeRadiologicalDepth(self._startVec, self._incVec, self._sourceVec, dispVec = (0,0,0) , ener = T0_round)
                


                                
                if not self._globalVar.staticDelivery and self._globalVar.compensationDynamic:
                ####################"
                #### Compensation : update the table for the repainting statistics
                ####################"
                    if T0_round not in dictRepaintStat[angle].keys():
                        dictRepaintStat[angle][T0_round] = dict()
                        compensatedPlan[angle][T0_round] = dict()
                        compensatedPlan[angle][T0_round]["ScanSpotPositionMap"]=list()
                        compensatedPlan[angle][T0_round]["ScanSpotMetersetWeights"]=list()                        
                    if dataSpot[8] not in dictRepaintStat[angle][T0_round].keys():
                        if  self._globalVar.findStartPos:
                        ####################"
                        #### Compensation : find a good starting position at each repainting
                        ####################"
                            startSpot = self.compensationFindStartSpotEnerLayer(dataBeam[0] , self._ddCurvesData ,dataSpot[6] )  
                            if startSpot is not None:
                                print "Starting Point found: %s"%str(startSpot[2])
                                dataSpot = startSpot                          
                        
                        dictRepaintStat[angle][T0_round][dataSpot[8]] = 0
                        #If margins added to scan spot position map, display it:
                        if self._globalVar.displayScanPathUpdates:
                            mapWithMargin = dataBeam[0].getMapWithMargin()
                            if mapWithMargin is not None:
                                viewer2D.displayImageFromArray( mapWithMargin, "Current Scan Spot Map with Margins - energy layer %i MeV"%T0_round)
                                filenameScanMap = "ScanSpotMapWithMargin%0.3dMeV_Repainting%0.3d"%(T0_round,dataSpot[8])
                                viewer2D.saveFigure( self._pathToOutputImages + "ScanSpotMapWithMargin/",filenameScanMap , ext = 'png', subplotIdx = 111)
                                viewer2D.closePlot()
                     
                    print "Current repainting: %i - nb comp spot so far: %i"%(dataSpot[8],dictRepaintStat[angle][T0_round][dataSpot[8]] )                
                    
                    
                if self._globalVar.displayScanPathUpdates:
                    mapWithCurrTarget = dataBeam[0].getMapWithEnhancedSpot(dataSpot[1])
                    viewer2D.displayImageFromArray( mapWithCurrTarget, "Current Scan Spot Position - energy layer %i MeV"%T0_round)
                    filenameScanPos = "ScanSpotPosition%0.9d_%0.3dMeV"%(idxSpot,T0_round)
                    viewer2D.saveFigure( self._pathToOutputImages + "ScanSpotPosition/",filenameScanPos , ext = 'png', subplotIdx = 111)
                    viewer2D.closePlot()                   
                    
                ###Spot delivery
                t_start_spot = datetime.datetime.now()
                weight = dataSpot[0]
                if weight < 1:
                    print "\t\tSpot %i , weight: %e"%(idxSpot+1,weight)
                else:
                    print "\t\tSpot %i , weight: %f"%(idxSpot+1,weight)
                
                xrshift = dataSpot[2] / 10.0# lateral shift of spot: from mm to cm
                yrshift = dataSpot[3] / 10.0#on-axis shift of spot: from mm to cm
                print "\t(xrshift,yrshift) : (%f,%f)"%(xrshift,yrshift)
                if listSpotsVisited is None:
                    listSpotsVisited = dict()
                if (xrshift,yrshift,T0) in listSpotsVisited.keys():
                    listSpotsVisited[(xrshift,yrshift,T0)] = listSpotsVisited[(xrshift,yrshift,T0)] + 1
                    print "spot (%f,%f) : number of visits : %i times"%( xrshift,yrshift, listSpotsVisited[(xrshift,yrshift,T0)])
                else:
                    listSpotsVisited[(xrshift,yrshift,T0)] = 1
                    
                
                if weight <= limitSmallWeights:
                    countSmallWeights += 1
                    logSpotsSmallWeights.append((xrshift,yrshift,T0,weight))
                else:
                    countNotSmallWeights += 1 
                    sumNotSmallWeights += weight
                    
                if weight <= sys.float_info.epsilon:
                     
                    if listSpotsDelivered is not None:
                        flagXYInList = 0
                        for indSpot in listSpotsDelivered:
                            if tuple(dataSpot[1]) == tuple(indSpot):
                                flagXYInList = 1
                                print "Spot with weight 0. Original weight %f."%dataBeam[0].getPlannedWeightForSpotInd(indSpot)
                                break
                        if flagXYInList == 1:
                            spotsVisitedButWeightNull = spotsVisitedButWeightNull + 1
                            print "Weight null - already visited"
                        else:
                            unexplainedWeightNull = unexplainedWeightNull + 1
                            print "Weight null - don't know why"
                    if (not self._globalVar.staticDelivery and not self._globalVar.compensationDynamic) or self._globalVar.staticDelivery:
                        continue
                ############# If Dynamic delivery #############
                if not self._globalVar.staticDelivery:
                    if dataSpot[6] == 'NewEnergy' :
                        prevSpot = [xrshift,yrshift]
                    if (prevSpot[1]) != yrshift and (prevSpot[0]) != xrshift:
                            deltaX = np.abs((prevSpot[0])-xrshift)
                            deltaY = np.abs((prevSpot[1])-yrshift)
                            dist = np.max( np.array([deltaX, deltaY]) )
                            self.updateTime('MoveToSpot',info=dist)                
                ############# If Dynamic delivery #############
                if not self._globalVar.staticDelivery:
                    ####################"
                    #### Compensation
                    ####################"
                
                    if self._globalVar.compensationDynamic:
                        print"Start compensation"
                        (dispVec_IECFixed ,dispVec_IECG,dispVec_CT,indDispCT) = self.computeDisplacementVectors(useMotionMonitor = True)['MeasuredMotion']
                        self.updateTime('Displacement Measure')
                        pointX = (xrshift - dispVec_IECG[0])*10.0 #convert from cm to mm
                        pointY = (yrshift - dispVec_IECG[1])*10.0 #convert from cm to mm
                        if self._globalVar.motion3D:
                            energyForCompens = self._ddCurvesData.findEnergyForBraggPeak(dmax + dispVec_IECG[2]) ## Double check depth
                            if energyForCompens is None:
                                print "Spot cannot be compensated: no energy achievable for BP %f cm!\n"%(dmax + dispVec_IECG[2])
                                dataBeam[0].updateSpot(dataSpot[1], 0)
                                nbMissedSpots = nbMissedSpots + 1
                                print "End Slice compensation"
                                continue
                            energyForCompens = dataBeam[0].findCorrespondingEnergy(energyForCompens)
                        else:
                            energyForCompens = T0
                        if energyForCompens is not None:
                            neighbor = dataBeam[0].isPointNearOtherSpot(energyForCompens,pointX, pointY)
                        else:
                            neighbor = [False]
                        if neighbor[0]:
                            locationNeighbors = dataBeam[0].locateSpots(([xrshift*10.0,pointX,neighbor[2][0]],[yrshift*10.0,pointY,neighbor[2][1]]))
                            viewer2D.displayImageFromArray( locationNeighbors, "Neighbours in energy slice")
                            filenameNeighbors = "Neighbors%0.9d_%0.3dMeV"%(idxSpot,T0_round)
                            viewer2D.saveFigure( self._pathToOutputImages + "Neighbors/",filenameNeighbors , ext = 'png', subplotIdx = 111)
                            viewer2D.closePlot()
                            
                            newWeight = dataBeam[0].getWeightForSpotInd(neighbor[1])
                            compensated = True
                            nbCompensatedSpots = nbCompensatedSpots + 1
                            dictRepaintStat[angle][T0_round][dataSpot[8]] += 1
                            print "Compensation is possible!"
                        else:
                            print "Spot cannot be compensated: no spot position found\n"
                            dataBeam[0].updateSpot(dataSpot[1], 0)
                            nbMissedSpots = nbMissedSpots + 1
                            print "End Slice compensation"
                            continue
                        
                        
                        print "End Slice compensation"
                        if compensated:
                            weight = newWeight
                            
                        
                        exec "timeToSpendAtSpot = self.timeDeliverSpotWeight(np.%s(weight))"%self._globalVar.typeFloat
                        if  timeToSpendAtSpot > 1e-3:
                            print "Time to spent at current spot is: %e s"%timeToSpendAtSpot


                    ####################"
                    #### End compensation
                    ####################"
                                
                ############# End - If Dynamic delivery #######
                ####################"
                ####    End  code for scanning path optimization
                ####################"                
                
                #Get number of loops to deliver each unit time
                if not self._globalVar.staticDelivery and self._globalVar.motionDuringBeam:
                    loopsDeliv = self.numberOfUnitTimeToDeliverWeight(weight)
                else:
                    loopsDeliv = (1 , 0 , weight)

                if compensated:
                    tupleIndices = tuple(neighbor[1])
                    tupleCoord = (neighbor[2][0]/10.0,neighbor[2][1]/10.0)
                    verifyCompensation = True
                    
                else: 
                    tupleIndices = tuple(dataSpot[1])
                    tupleCoord = (xrshift,yrshift)
                if listSpotsDelivered is not None:
                    if tupleIndices in listSpotsDelivered.keys():
                        listSpotsDelivered[tupleIndices]['visits'] = listSpotsDelivered[tupleIndices]['visits'] + 1
                        listSpotsDelivered[tupleIndices]['sumWeight'] = listSpotsDelivered[tupleIndices]['sumWeight'] + weight
                    else:
                        listSpotsDelivered[tupleIndices] = { 'visits':1 , 'sumWeight':weight ,'coord':tupleCoord}
                else:
                    listSpotsDelivered = dict()
                    listSpotsDelivered[tupleIndices] = { 'visits':1 , 'sumWeight':weight ,'coord':tupleCoord}
                    
                print "Nb of time units for current weight: %i"%loopsDeliv[0]
                if not self._globalVar.staticDelivery:
                    prevDispIndices = None
                
                if not self._globalVar.skipSaveBeamlets:
                    beamletDoseDistrib = np.zeros(self._DensityGrid.shape, dtype = self._globalVar.typeFloat, order = 'C')
                
                
                for loopIdx in xrange(loopsDeliv[0]):   
                    print "\t----- # of time units: %i / %i"%(loopIdx + 1, loopsDeliv[0])
                    
                    if loopIdx < loopsDeliv[0] - 1:
                        weightLoop = loopsDeliv[1]
                    else:
                        weightLoop = loopsDeliv[2]
                    ############## If Dynamic delivery #############
                    if not self._globalVar.staticDelivery:
                        
                        (dispVec_IECFixed ,dispVec_IECG,dispVec_CT,indDispCT) = self.computeDisplacementVectors()
                         
                        self._dynaPatientLoc = self._patient.getPatientExtentInDynamicModel(dispVec_CT)
                        if not  self.isPatientInTheVolume(self._patient.shape(),self._dynaPatientLoc):
                            strErr = "Shifted patient is outside the volume defined "
                            raise ValueError(strErr)
                        self._StPwrRatioGrid = self.funcGetDataPatient('stpPwrRatio',dispVec_CT)
                        self._DensityGrid = self.funcGetDataPatient('density',dispVec_CT) 
        
                        if prevDispIndices is None or (prevDispIndices is not None and (prevDispIndices[0] != indDispCT[0] or \
                                                                                      prevDispIndices[1] != indDispCT[1] or \
                                                                                      prevDispIndices[2] != indDispCT[2])):
                            targetIECG = np.array([xrshift,yrshift,0],dtype= self._globalVar.typeFloat,order='C')# target point in isocenter frame
                            targetIECF =  mathTools.rotateCoordinates(self._rotMatrix,targetIECG, self._globalVar.typeFloat)
                            sourceIECG = np.array([0,0,sourceLocation[2]],dtype= self._globalVar.typeFloat,order='C')
                            if useFastRadiolDepth_DynaDeliv:
                                nbPositionsX = 2
                                nbPositionsY = 2
                                self._radiolDepth,ssd0,listIndices = self._patient.computeLocalizedRadiologicalDepth(self._startVec, self._incVec, self._sourceVec,\
                                                        sourceIECG,self._rotMatrix,Xr,Yr,Zr,xrshift,yrshift, \
                                                        nbPositionsX = nbPositionsX,nbPositionsY = nbPositionsY, dispVec = dispVec_CT, ener = T0_round)
                            else:
                                self._radiolDepth = self._patient.computeRadiologicalDepth(self._startVec, self._incVec, self._sourceVec, dispVec = dispVec_CT , ener = T0_round)
                                patientMask = np.zeros(self._DensityGrid.shape,dtype= 'int8',order='C') 
                                patientMask[self._dynaPatientLoc] = 1
                                ssd0, listIndices,listLenVoxels = self._patient.calculateSSD0(self._sourceVec,targetIECF,Xr,Yr,Zr,self._DensityGrid,self._startVec,self._incVec,sourceIECG,patientMask)
                            
                            prevDispIndices = indDispCT
                        else:
                            print "No need to update radiol depth: indices of the displacement vector didn't change."
                        
                        #In dynamic the first computation of self._subVol returned None and X,Y,Z, X_iec, Y_iec, Z_iec and Xr, Yr,Zr have been computed for the entire volume.
                        Xr_tmp = Xr 
                        Yr_tmp = Yr 
                        Zr_tmp = Zr 
                        print "Displacement vector (f,r,c) in CT image:%s"%str(indDispCT)
                        
                    ############ End - If Dynamic delivery ############# 
                    else: #Static
                        Xr_tmp = Xr
                        Yr_tmp = Yr
                        Zr_tmp = Zr
                        dispVec = None
                        patientMask = self._patient.getPatientMask()
                        targetIECG = np.array([xrshift,yrshift,0],dtype= self._globalVar.typeFloat,order='C')# target point in isocenter frame
                        targetIECF =  mathTools.rotateCoordinates(self._rotMatrix,targetIECG, self._globalVar.typeFloat)
                        sourceIECG = np.array([0,0,sourceLocation[2]],dtype= self._globalVar.typeFloat,order='C')
                        ssd0, listIndices,listLenVoxels = self._patient.calculateSSD0(self._sourceVec,targetIECF,Xr_tmp,Yr_tmp,Zr_tmp,self._DensityGrid,self._startVec,self._incVec,sourceIECG,patientMask)
                    
                    t_start_beamlet = datetime.datetime.now()
                
                    
                    computedDoseRaw = beamletSimul(self._DensityGrid,self._StPwrRatioGrid,\
                                     Xr_tmp,Yr_tmp,Zr_tmp,np.array(dmax,dtype=self._globalVar.typeFloat),\
                                     z_dd,dd,np.array(sig0,dtype=self._globalVar.typeFloat),\
                                     np.array(xrshift ,dtype=self._globalVar.typeFloat),\
                                     np.array(yrshift ,dtype=self._globalVar.typeFloat),np.array(sourceLocation[2],\
                                     dtype=self._globalVar.typeFloat),np.array(ssd0 ,dtype=self._globalVar.typeFloat),self._radiolDepth,listIndices,self._x0y0z0,self._IECFSpacing,self._rotMatrix,self._globalVar.typeFloat)


                    t_end_beamlet = datetime.datetime.now()
                    print "Beamlet simul in %s"%str(t_end_beamlet - t_start_beamlet)
                    
                
                    if not self._globalVar.staticDelivery:
                        if not self._globalVar.compensationDynamic:
                            indexForVerification = idxSpot
                        else:
                            indexForVerification = nbCompensatedSpots
                    
                    if not self._globalVar.staticDelivery and self._globalVar.evaluateDiffStatDyna and loopIdx == loopsDeliv[0]-1 and (indexForVerification == 1 or indexForVerification%10 == 0):
                        if compensated:
                            if verifyCompensation: 
                                (xStat , yStat) = (neighbor[2][0]/10.0,neighbor[2][1]/10.0)
                                nameEvaluationText  = "EvaluationCompensation.txt"
                                outputEvalImages = 'CompensatedBeamlet'
                                titleFig = 'Compensation '
                                if self._globalVar.motion3D:
                                    dmax_stat = self._ddCurvesData.dataForEnergy('BraggPeak' , energyForCompens)
                                    z_dd_stat = self._ddCurvesData.dataForEnergy('DD_Curve' ,  energyForCompens)[0,:].astype(self._globalVar.typeFloat)
                                    dd_stat = self._ddCurvesData.dataForEnergy('DD_Curve' , energyForCompens)[1,:].astype(self._globalVar.typeFloat)
                                else:
                                    dmax_stat = dmax
                                    z_dd_stat = z_dd
                                    dd_stat = dd
                        elif not self._globalVar.compensationDynamic:
                            (xStat , yStat) = (xrshift,yrshift)
                            nameEvaluationText  = "EvaluationNoCompensation.txt"
                            outputEvalImages = 'NotCompensatedBeamlet'
                            titleFig = 'Motion '
                            dmax_stat = dmax
                            z_dd_stat = z_dd
                            dd_stat = dd
                        else:
                            raise ValueError("in dynamic delivery - wrong choice for result evaluation!")
                        
                        
                        
                        beamletStatic = self.computeBeamletStatic(Xr,Yr,Zr,dmax_stat,z_dd_stat,dd_stat,sig0,xStat  ,yStat,sourceLocation)
                        
                        if not self._globalVar.motion3D:

                            tmpComputedDoseRaw = np.zeros(computedDoseRaw.shape) 
                            indicesCompDose = np.where(computedDoseRaw > 0)
                            newIndCompDose = tuple((indicesCompDose[0]-indDispCT[0],indicesCompDose[1]-indDispCT[1],indicesCompDose[2]-indDispCT[2]))
                            tmpComputedDoseRaw[newIndCompDose] = computedDoseRaw[indicesCompDose]
                        
                            dataStatic = beamletStatic
                            dataDyna = tmpComputedDoseRaw
                        else:
                            
                            
                            staticLoc = self._patient.getPatientExtentInStaticModel()
                                                 
                            dataStatic = np.zeros(computedDoseRaw.shape)
                            dataDyna = np.zeros(computedDoseRaw.shape)
                            
                            patientLoc = (np.min(staticLoc[0]) + self._dynaPatientLoc[0]-np.min(self._dynaPatientLoc[0]) ,\
                                        np.min(staticLoc[1]) + self._dynaPatientLoc[1]-np.min(self._dynaPatientLoc[1]) ,\
                                        np.min(staticLoc[2]) + self._dynaPatientLoc[2]-np.min(self._dynaPatientLoc[2]) )
                            
                            dataStatic[patientLoc] = beamletStatic[staticLoc]
                            dataDyna[patientLoc] = computedDoseRaw[self._dynaPatientLoc]
                        
                        print "\t*********Beamlet evaluation*********"
                        valuesEval = self._tool.evaluateDifference(dataStatic, dataDyna)
                        print "\t******************"
                        
                        self._tool.showHotColdSpots(dataStatic, dataDyna, title = 'Hot and cold spots' , pathToSaveImg = self._globalVar.mainPath + self._globalVar.nameSimulation + 'HotColdSpotsBeamlets/'+'HotAndColdSpots_Spot%0.9d_loop%0.9d'%(idxSpot,loopIdx))
                        
                        dataHotColdSpots = self._patient.hotAndColdSpotsPerStrcuture(dataStatic, dataDyna, unpadded = False)
                        for roi in dataHotColdSpots:
                            print "%s :\n\t#of voxels: %i , \n\t#of hot spot:%i ( %0.3f %%) , \n\t#of cold spots:%i (%0.3f %%), \n\t#of correct spots:%i (%0.3f %%)"%(roi,dataHotColdSpots[roi]['nbVoxels'],\
                            dataHotColdSpots[roi]['nbHotSpots'],\
                            dataHotColdSpots[roi]['percHotSpots'],\
                            dataHotColdSpots[roi]['nbColdSpots'],\
                            dataHotColdSpots[roi]['percColdSpots'],\
                            dataHotColdSpots[roi]['nbCorrectSpots'],\
                            dataHotColdSpots[roi]['percCorrectSpots'])
                        
                        
                        if (compensated and nbCompensatedSpots == 1 )or (not self._globalVar.compensationDynamic and idxSpot == 0): 
                            fileOption = 'wb'
                            textToWrite = 'Ener\tSpot\tShiftX\tShiftY\tShiftZ\tSens\tSpec\tJac'+\
                            '\tMinStat\tMaxStat\tAvgStat\tStdStat'+\
                            '\tMinDyna\tMaxDyna\tAvgDyna\tStdDyna\tDiffMax\tDiffAvg\tDiffStd\n'
                        else:
                            fileOption = 'a'
                            textToWrite = ''
                        fileOutput = self._globalVar.mainPath + self._globalVar.nameSimulation + nameEvaluationText
                        
                        textHeader ='%0.2f\t%i\t%0.3f\t%0.3f\t%0.3f'%(T0,idxSpot,dispVec_CT[0],dispVec_CT[1],dispVec_CT[2])
                        resEvalValues = '\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n'%(valuesEval['vol1']['min'],valuesEval['vol1']['max'],valuesEval['vol1']['avg'],valuesEval['vol1']['std'],valuesEval['vol2']['min'],valuesEval['vol2']['max'],valuesEval['vol2']['avg'],valuesEval['vol2']['std'],np.abs(valuesEval['vol1']['max']- valuesEval['vol2']['max']),np.abs(valuesEval['vol1']['avg']- valuesEval['vol2']['avg']),np.abs(valuesEval['vol1']['std']- valuesEval['vol2']['std']))
                        
                        textToWrite =  textToWrite + textHeader + resEvalValues
                        with open(fileOutput, fileOption) as myFile:
                            myFile.write(textToWrite)
                        
                        
                        indStat = np.where(dataStatic > 0)
                        frameStat = int((np.max(indStat[0])+np.min(indStat[0]))/2)
                        if self._globalVar.motion3D:
                            indDyna = np.where(dataDyna > 0)
                            frameDyna = int((np.max(indDyna[0])+np.min(indDyna[0]))/2)
                        else:
                            frameDyna = frameStat
                            
                        viewer2D.display2ImagesFromArray( dataStatic[frameStat], dataDyna[frameDyna], "Beamlet Static: target(%0.3f,%0.3f) - frame:%i"%(xStat , yStat,frameStat), "Beamlet Dynamic: target(%0.3f,%0.3f) - frame:%i"%(xrshift,yrshift,frameDyna),titleFig+" observation")
                        filenameBeamlets = outputEvalImages+"%0.9d_frames%0.4d_and_%0.4d"%(idxSpot,frameStat,frameDyna)
                        viewer2D.saveFigure( self._pathToOutputImages + outputEvalImages+"s/",filenameBeamlets , ext = 'png')
                        viewer2D.closePlot()
                        if compensated:  
                            if verifyCompensation:
                                verifyCompensation = False

                    ############# If Dynamic delivery #############
                    if not self._globalVar.staticDelivery:
                        self.updateTime('Delivery',info = weightLoop)
                        if self._globalVar.storeInterDoseDens:
                            dVec = dispVec_IECG
                        else:
                            dVec = None
                        computedDose = self.processDynamicDoseDistrib(computedDoseRaw, T0_round, dispVec = dVec)
                    ############# End - If Dynamic delivery #######
                    else:
                        
                        computedDose = computedDoseRaw
    
                    #Find Bragg Peak
                    self.findBraggPeak(computedDose, self._radiolDepth)
        
                    correcFactor =  self._ddCurvesData.doseCorrectionFactor(T0)

                    weightedDose = (weightLoop * computedDose * (self._protPerMU / (nbProtons)) *correcFactor )
                    
                    if not self._globalVar.skipSaveBeamlets:
                        beamletDoseDistrib = beamletDoseDistrib + weightedDose

                
                    self._receivedDose = self._receivedDose +  weightedDose
                    
                #For info:
                nbTargets = nbTargets + 1
                sumWeights = sumWeights + weight

                if compensated:
                    dataBeam[0].updateSpot(neighbor[1], weight, decrRepaint = False)
                    dataBeam[0].updateSpot(dataSpot[1], 0)
                    compensatedPlan[angle][T0_round]["ScanSpotPositionMap"].append(dataSpot[2])
                    compensatedPlan[angle][T0_round]["ScanSpotPositionMap"].append(dataSpot[3])
                    compensatedPlan[angle][T0_round]["ScanSpotMetersetWeights"].append(weight)
 
                else:
                    dataBeam[0].updateSpot(dataSpot[1], weight)
                    print "Update done.."
                if not self._globalVar.skipSaveBeamlets:
                    indImg = np.where(beamletDoseDistrib > 0)
                    try:
                        frameImg = int((np.max(indImg[0])+np.min(indImg[0]))/2)
#                     print "Dose not empty"
                    except:
                        print "Computed dose empty"
                        raise ValueError("Computed dose empty")
                
                    viewer2D.displayImageFromArray( beamletDoseDistrib[frameImg], "Beamlet %i"%idxSpot)
                    filenameBeamlets = "Beam_%0.2d_Beamlet%0.9d_frame%0.4d"%(idxBeam,idxSpot,frameImg)
                    viewer2D.saveFigure( self._pathToOutputImages + "Beamlets/",filenameBeamlets , ext = 'png', subplotIdx = 111)
                    viewer2D.closePlot()
                    name = "Beam_%0.2d_Beamlet%0.9d_E%iMeV_frame%0.4d"%(idxBeam,idxSpot,T0_round,frameImg)
                    viewer2D.storeMatrixAsImage(self._pathToOutputImages,beamletDoseDistrib[frameImg].reshape(1,computedDose[frameImg].shape[0],computedDose[frameImg].shape[1]),name, alongAxis = None,cpkl = self._globalVar.saveBeamletsCPKL, bw = False)
            
                
                t_end_spot = datetime.datetime.now()
                print "Spot # %i processed in %s"%(idxSpot+1,str(t_end_spot - t_start_spot))            
                
            
            #Display resulting scanning path "matrix"
            if self._globalVar.displayScanPathUpdates:
                    listEnergies = dataBeam[0].getListEnergies()
                    #Display the current scanning path slice
                    for T0 in listEnergies:
                        currentEnergyLayer = dataBeam[0].getEnergySliceWeights(T0)
                        viewer2D.displayImageFromArray( currentEnergyLayer, "Scan path - energy layer %i MeV"%np.round(T0))
                        filenameScanPath = "FinalScanningPath_%0.3dMeV"%(np.round(T0))
                        viewer2D.saveFigure( self._pathToOutputImages + "ScanningPathUpdates/",filenameScanPath , ext = 'png', subplotIdx = 111)
                        viewer2D.closePlot()
            
            
            t_end_beam = datetime.datetime.now()
            print "Beam # %i processed in %s "%(idxBeam+1, str(t_end_beam-t_start_beam))
            if listSpotsVisited is not None and dataBeam[0].repaintingMethodSlice == 'ScaledRepaint':
                print "#### Checking repainting ####"
                flagRepaint = 0
                repaintFactor = dataBeam[0].getRepaintingFactor()
                print "Repainting factor used: %i"%repaintFactor
                print "Number of unique spots: %i"%len(listSpotsVisited.keys())
                for visitedSpots in  listSpotsVisited:
                    if listSpotsVisited[visitedSpots] > repaintFactor:
                        flagRepaint = 1
                        print "WARNING : spot %s visited %i times instead of %i "%(str(visitedSpots),listSpotsVisited[visitedSpots],repaintFactor)
                if flagRepaint == 0:
                    print "Nothing wrong has been detected in the repaintings."

            if listSpotsDelivered is not None and dataBeam[0].repaintingMethodSlice == 'ScaledRepaint':
                print "### checking weights deposited ###"
                flagHotSpot = 0
                flagVisit = 0
                for indSpot in listSpotsDelivered:
                    expectedWeight = dataBeam[0].getPlannedWeightForSpotInd(indSpot)
                    if np.abs(expectedWeight-listSpotsDelivered[indSpot]['sumWeight']) > sys.float_info.epsilon:
                        flagHotSpot = 1
                        print "WARNING : spot %s : expected weight: %f - delivered weight: %f"%(str(listSpotsDelivered[indSpot]['coord']),expectedWeight,listSpotsDelivered[indSpot]['sumWeight'])
                    if listSpotsDelivered[indSpot]['visits'] > repaintFactor:
                        flagVisit = 1
                        print "WARNING: spot %s : delivered %i times"%(str(listSpotsDelivered[indSpot]['coord']),listSpotsDelivered[indSpot]['visits'])
                if flagVisit == 0:
                    print "Nothing wrong has been detected in the number of times spots have been delivered."
                if flagHotSpot == 0:
                    print "Nothing wrong has been detected in the spots weight."
            print "#### Max Number of repaintings per energy slice ####"
            for indSlice,repaintNb in enumerate(dataBeam[0].getNumberOfRepaintingsPerSlice()):
                print "slice %i: %i repainting(s)"%(indSlice,repaintNb)
            print "#### Checking spots with weight Null ####"
            print "Spots with weight Null because already visited: %i"%spotsVisitedButWeightNull 
            print "Spots with unexplained weight null: %i"%unexplainedWeightNull
            print "################################"
                
                
                
                
        print "End beamlet calculation..."
        t_simul_end = datetime.datetime.now()
        totalTimeSimul = t_simul_end-t_simul_start
        print "Total delivery simulation in %s"%(str(totalTimeSimul) )
        self._patient.setReceivedDose(self._receivedDose)
        
        ############# If Dynamic delivery #############
        if not self._globalVar.staticDelivery:
            print "The delivery lasted (simulated time): %f min"%(self.getCurrentTime()/60.0)
            print "The time spent irradiating is: %f min"%(self.getCumulDelivDuration()/60.0)
        ############# End - If Dynamic delivery #######
        
        if countSmallWeights > 0:
            print ">>>>>>>>>>>>>>>>>>>>>>>>>>>"
            print "List of spots (%i) Null ( <= %e ) that shouldn't:"%(countSmallWeights,limitSmallWeights)
            sumSmallWeights = 0
            for idx,tup in enumerate(logSpotsSmallWeights):
                sumSmallWeights+= tup[3]
                print "Spot %i: %s"%(idx,str(tup))
            print "Sum weights: %f"%sumSmallWeights
            print ">>>>>>>>>>>>>>>>>>>>>>>>>>>"
            print "Not small weights: %i spots : %f MU"%(countNotSmallWeights,sumNotSmallWeights)


        
        print "Average weight : %f "%(sumWeights/nbTargets)       
        print "Sum weights : %f"%sumWeights
        print "Total number of spots: %i"%nbSpots
        if not self._globalVar.staticDelivery and self._globalVar.compensationDynamic:
            textToWrite = "*********************************\nSummary compensation:\nNumber of missed spots:%i\nNumber of compensated spots:%i"%(nbMissedSpots,nbCompensatedSpots)
            print textToWrite
            print "Statistics on repaintings:\n%s"%str(dictRepaintStat)
            self.storeDictStatisticsRepaint(dictRepaintStat)
            return compensatedPlan
        else:
            return None
            
        
       
#############################################################################
###         Functions related to dynamic delivery
#############################################################################

    def updateTime(self,event,info = None):
        '''For dynamic delivery: Update the delivery timer (expressed in seconds).
        
        :param event: The type of event that requires a timer update:
            
            * 'Reset': Reset the timer
            * 'Energy': Energy change
            * 'MoveToSpot': Change beam position
            * 'GantryRotation': Gantry rotation
            * 'Delivery': Time spent at a beam position for the delivery
            * 'Refill': Time to refill the synchrotron ( if the machine setting is set to Synchrotron)
            * 'Displacement Measure': Time to measure the patient displacement (if Compensation is activated)
            
        :param info: Information regarding the event
            
            * If event 'MoveToSpot', 'info' is the distance between the spots in cm.
            * If event 'GantryRotation', 'info' is the rotation angle in degrees.
            * If event 'Delivery', 'info' is the weight in MU.
            * Otherwise info is None  
        
            
        '''
        print "\t\t\t---> Update time for event '%s'"%event
        prevTime = self._time
        if event == 'Reset': #Reset delivery timer
            self._time = 0 #Time of the entire delivery
            self._cumulDelivSpot = 0 # Cumulative delivery time at each spot (time spent over a spot)
            self._spillTimer = 0 # Time incremented every time a spot is irradiated and every time the beam moves from ones spot to another. This time must always stay below the maximal spill length, otherwise (with a synchrotron) it means we have to refill the tank.
            #print"reset Done"
            
        elif event == 'Energy':#Energy change
            self._time = self._time + self._globalVar.timeEnergyChange
        elif event == 'MoveToSpot': #Move to another spot: info shouldn't be None, it should be a number representing the distance in cm.
            if info is None:
                raise ValueError("In updateTime for delivery, information expected to update time in order to move between 2 spots: distance in cm !")
            if not isinstance( info, ( int, long, float , np.float64,np.float32) ):
                strErr = "In updateTime for delivery, move between 2 spots: information expected is distance (number) in cm, not a '%s' - %s! "%(type(info),str(info))
                raise ValueError(strErr) 
            exec "timeMove = (np.%s(info) / self._globalVar.lateralScanningSpeed) + self._globalVar.spotSettlingTime"%self._globalVar.typeFloat
            self._time = self._time +  timeMove 
            self._spillTimer = self._spillTimer + timeMove
            print "\tTime to move to spot:%f ms"%(timeMove * 1e3)
        elif event == 'GantryRotation': #Beam rotation
            if info is None:
                raise ValueError("Update time for gantry rotations needs the angle")
            if not isinstance( info, ( int, long, float , np.float64,np.float32) ):
                strErr = "In updateTime for gantry rotation, information expected is an angle (number) indegrees, not a '%s' - %s! "%(type(info),str(info))
                raise ValueError(strErr)
            exec "testInfo = np.%s(info)"%self._globalVar.typeFloat
            if testInfo < 0:
                exec "newAngle = 360 - np.%s(info)"%self._globalVar.typeFloat
            else:
                exec "newAngle = np.%s(info)"%self._globalVar.typeFloat
            res = newAngle /  360.0
            if newAngle >  360.0 :
                newAngle = (res - np.floor(res)) * 360.0
            
            angleRot = np.abs(newAngle - self._prevRotationAngle)
            self._time = self._time + (angleRot/self._globalVar.speedGantryRotation)
            exec "self._prevRotationAngle = np.%s(newAngle)"%self._globalVar.typeFloat
        elif event == 'Delivery': # Time spent at one spot: info shouldn't be None, it should be a number representing the weight in MU.
            if info is None:
                raise ValueError("In updateTime for delivery, information expected to update time in order to deliver at 1 spot: weight in MU!")
            if not isinstance( info, ( int, long, float, np.ndarray , np.float32, np.float64) ):
                print "info: %s"%str(info)
                strErr = "In updateTime for delivery, deliver at 1 spot: information expected is weight (number) in MU, not a '%s'! "%type(info)
                raise ValueError(strErr)
            else:
                if isinstance(info, np.ndarray) and len(np.atleast_1d(info)) != 1:
                    print "info: %s"%str(info)
                    strErr = "In updateTime for delivery, deliver at 1 spot: information expected is weight (number or numpy  array of 1 element) in MU, not a '%i' elements! "%len(np.atleast_1d(info))
                    raise ValueError(strErr)                   
            exec "timeSpent = self.timeDeliverSpotWeight(np.%s(info))"%self._globalVar.typeFloat
            self._time = self._time + timeSpent
            self._cumulDelivSpot = self._cumulDelivSpot + timeSpent
            self._spillTimer = self._spillTimer + timeSpent
            print "\tTime spent over spot:%f ms"%(timeSpent * 1e3)
        elif event == 'Refill':
            self._time = self._time + self._globalVar.timeEnergyChange
        elif event == 'Displacement Measure':
            self._time = self._time + self._globalVar.timeMeasureDisplacement
        
        
        else:
            strErr = "In updateTime for delivery, unknowkn event '%s'! "%str(info)
        
        if self._spillTimer >= self._globalVar.spillDuration and self._globalVar.synchrotron:
            self._spillTimer = self._spillTimer - self._globalVar.spillDuration
            self.updateTime('Refill') # in order to simulate the refill
        print "\t\t\t\tUpdated timer: %f ms"%(self._time)
            
            
    def getCurrentTime(self):
        '''Get the current value of the timer (dynamic delivery)
        
        :returns: Current timer value
        
        '''
        return self._time
    
    def getCumulDelivDuration(self):
        '''Get the cumulative time spent at each beam position
        
        :returns: The cumulative irradiation time.
        
        '''
        return self._cumulDelivSpot
    
    def timeDeliverSpotWeight(self,weight):
        '''Return the time spent at one spot to deliver a given weight (MU), knowing the beam intensity and the number_of_protons/MU. 
        
        Formula : eq. 1 in "Interplay effects in proton scanning for lung: a 4D Monte Carlo study assessing the impact of tumor and beam delivery parameters." Dowdell et al. 2013
        
        :param weight: Weight in MU of the beam position.
        
        :returns: Time spent at beam position
        
        '''
        return (weight * self._protPerMU / self._globalVar.beamIntensity)


    def processDynamicDoseDistrib(self, computedDose,  T0_round ,dispVec = None):
        '''Process the computed dose in dynamic mode:
        
        The dose received by the patient (located at its dynamic position: self._dynaPatientLoc) is translated toward the "static position" (self._staticPatientLoc). \
        The dose values located at indices self._dynaPatientLoc (body dynamic position) are moved to indices self._staticPatientLoc (body static position). \
        Then everything outside the patient body is set to zero.
       
        :param computedDose: 3D numpy array representing the 3D dose distribution of the pencilbeam.
        :param T0_round: Rounded beam energy
        :param dispVec: The displacement vector. By default it is None. It is used if the global variable 'storeInterDoseDens' is set to True when configuring MSPT.
        
        :returns: The processed dose distribution
        
        '''
        tmpComputedDose =  computedDose
        tmpDoseToStore = np.zeros(computedDose.shape,dtype=self._globalVar.typeFloat)
        tmpDoseToStore[self._dynaPatientLoc] = tmpComputedDose[self._dynaPatientLoc]
        if self._globalVar.storeInterDoseDens and dispVec is not None:
            nameOutputDose = "interCompDose_[%.2f,%.2f,%.2f]_%iMeV"%(dispVec[0],dispVec[1],dispVec[2],T0_round)
            nameOutputDensity = "interDensity_[%.2f,%.2f,%.2f]_%iMeV"%(dispVec[0],dispVec[1],dispVec[2],T0_round)
            for axis in self._globalVar.listSaveCompDoseAlong:
                viewer2D.storeMatrixAsImage(self._pathToOutputImages,tmpDoseToStore,nameOutputDose,alongAxis = axis ,cpkl=False,bw=False)
                viewer2D.storeMatrixAsImage(self._pathToOutputImages,self._DensityGrid,nameOutputDensity,alongAxis = axis ,cpkl=False,bw=True)
        saveDoseTumor = tmpComputedDose[self._dynaPatientLoc]
                     
        #Otherwise:
        newCompDose = np.zeros(computedDose.shape,dtype=self._globalVar.typeFloat)
        newCompDose[self._staticPatientLoc] = np.asarray(saveDoseTumor,dtype=self._globalVar.typeFloat,order='C')                        
                
        return newCompDose
    
    
    def computeBeamletStatic(self,Xr,Yr,Zr,dmax,z_dd,dd,sig0,xrshift,yrshift,sourceLocation):
        '''This function is only used if MSPT global variables 'staticDelivery' is set to False and 'evaluateDiffStatDyna' is set to False.
        
        It simulates the irradiation at a beam position on the static body. This allows to highlight the differences between the beams delivered to \
        the static and the dynamic body for a same beam position.
        
        :param Xr: 3D matrix representing x coordinates in the IEC Gantry coordinate system (cm).
        :param Yr: 3D matrix representing y coordinates in the IEC Gantry coordinate system (cm).
        :param Zr: 3D matrix representing z coordinates in the IEC Gantry coordinate system (cm).
        :param dmax: Bragg Peak depth (cm)
        :param z_dd: "depth" coordinates for the depth dose curve (cm).
        :param dd: "dose" coordinates for the depth dose curve. (Gy.cm ^2 )
        :param sig0: beam size (along x and y IEC gantry axis) at patient entrance (cm): [sig0X,sig0Y]
        :param xrshift: x beam shift (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
        :apram yrshift: y beam shift (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
        :param sourceLocation: Source location in the IEC Fixed coordinate system(cm).

        :returns: The processed dose distribution
        
        '''
        if not self._globalVar.staticDelivery:
            stpPwrStatic = np.ascontiguousarray(self._patient.getPatientArrayForAttr('stpPwrRatio'))
            densityStatic = np.ascontiguousarray(self._patient.getPatientArrayForAttr('density'))
            
            if self._deffStatic is None:
                self._deffStatic = self._patient.computeRadiologicalDepth(self._startVec, self._incVec, self._sourceVec, dispVec = (0,0,0) )
            
            patientMask = self._patient.getPatientMask()
            targetIECG = np.array([xrshift,yrshift,0],dtype= self._globalVar.typeFloat,order='C')# target point in isocenter frame
            targetIECF =  mathTools.rotateCoordinates(self._rotMatrix,targetIECG, self._globalVar.typeFloat)
            sourceIECG = np.array([0,0,sourceLocation[2]],dtype= self._globalVar.typeFloat,order='C')
            ssd0, listIndices,listLenVoxels = self._patient.calculateSSD0(self._sourceVec,targetIECF,Xr,Yr,Zr,self._DensityGrid,self._startVec,self._incVec,sourceIECG,patientMask)
 
            computedDose = beamletSimul(densityStatic,stpPwrStatic,\
                                Xr,Yr,Zr,np.array(dmax,dtype=self._globalVar.typeFloat),\
                                z_dd,dd,np.array(sig0,dtype=self._globalVar.typeFloat),\
                                np.array(xrshift,dtype=self._globalVar.typeFloat),\
                                np.array(yrshift,dtype=self._globalVar.typeFloat),np.array(sourceLocation[2],\
                                dtype=self._globalVar.typeFloat),np.array(ssd0,\
                                dtype=self._globalVar.typeFloat),self._deffStatic,listIndices,self._x0y0z0,self._IECFSpacing,self._rotMatrix, self._globalVar.typeFloat)
            return computedDose
        else:
            raise Exception("Trying to compute beamlet static designed for dynamic mode, in static mode... ")
        return None
        
        
    def isPatientInTheVolume(self,shapeVolume,patientIndices):
        '''Function to verify that the patient is in the field of view when moving.
        
        :param shapeVolume: Shape of the field of view volume (nb frames, nb rows, nb columns).
        :param patientIndices: Patient indices in the field of view 3D array.
        
        :returns: True if the patient is inside, False otherwise
        
        '''

        nbFrames = shapeVolume[0]
        nbRows = shapeVolume[1]
        nbCols = shapeVolume[2]
        
        #Check if patient in field of view:
        if np.all((patientIndices[0]>= 0) & (patientIndices[0]< nbFrames)) and \
        np.all((patientIndices[1]>= 0) & (patientIndices[1]< nbRows)) and \
        np.all((patientIndices[2]>= 0) & (patientIndices[2]< nbCols)):
            return True
        else:
            return False


    def numberOfUnitTimeToDeliverWeight(self,weight):
        '''Calculate the number of time units needed to deliver a given weight:
        
        :param weight: Beam weight in MU
        
        :returns: A list [ number of time units , weight to be delivered per time unit, remainder to deliver at last time unit]
        '''
        
        weightPerTimeUnit = self.weightDeliveredByUnitTime()
        if weight < weightPerTimeUnit:
            return (1 , 0 , weight)
        print "Weight per time unit: %f"%weightPerTimeUnit
        nbLoops = np.ceil(weight / weightPerTimeUnit)
        remainWeight = weight - ((nbLoops -1 ) *  weightPerTimeUnit)
        return (int(nbLoops),weightPerTimeUnit,remainWeight)
        
        
    def weightDeliveredByUnitTime(self):
        '''Calculates the weight delivered per time unit
        
        :returns: The calculated weight
        
        '''
        return (self._globalVar.beamIntensity * self._globalVar.unitTimeForDelivery) / self._protPerMU 



    def computeDisplacementVectors(self , useMotionMonitor = False):
        '''Compute the displacement vector in the IEC fixed coordinate system from the simulated timer. \
        Then convert the coordinated into the IEC Gantry coordinate system and into the CT coordinate system.\
        
        :param useMotionMonitor: True if MSPT is configured to use motion monitoring, False otherwise. 
        
        
        :returns: If useMotionMonitor is True:
        
            * A tuple T with:
                
                * tuple[0] : disp. vector in the IEC Fixed system
                * tuple[1] : disp vector in the IEC Gantry
                * tuple[2] : disp vector in the CT system
                * tuple[3] : disp indices in the CT image (matrix)
            
            * Else, a dictionary with keys:
                
                * 'PatientMotion' : the tuple T
                * 'MeasuredMotion' : a tuple equivalent to T corresponding to the monitor measurements.
            
        '''
        currTime = self.getCurrentTime()
        if not useMotionMonitor:
            dispVec_IECFixed = self._motion.getDisplacementVectorAtTime(currTime)
            return self.getVectorInDiffCoordSystm(dispVec_IECFixed)

        else:
            
            motion = self._monitoringSystm.monitorPatientMotionAtTime(currTime)
            vecMotion = motion['PatientMotion']
            vecMeasure = motion['MeasuredMotion']
            
            vecInCoordSystmMotion = self.getVectorInDiffCoordSystm(vecMotion,' real motion ')
            vecInCoordSystmMeasure = self.getVectorInDiffCoordSystm(vecMeasure,' measured motion ')
            return {'PatientMotion': vecInCoordSystmMotion, 'MeasuredMotion': vecInCoordSystmMeasure}
            
            
    def getVectorInDiffCoordSystm(self,dispVec_IECFixed, info = ''):
        '''Converts a vector defined in the IEC Fixed coordinate system into the IEC grantry coordinate system and the dicom coordinate system.
        
       :param dispVec_IECFixed: Displacement vector in the IEC fixed coordinate system (x,y,z)
       :param info: String for information, to display, about the motion.
       
       :returns: A tuple with the displacement vector in different coordinate systems: \
       ( dispVec_IECFixed ,dispVec_IECG,dispVec_CT,indDispCT ). Note: indDispCT corresponds to the displacement vector as indices \
       in the CT volume: (frame, row, colum).
       
        '''
        
        dx_IEC_F = np.array( dispVec_IECFixed[0]).reshape((1,1,1))
        dy_IEC_F = np.array( dispVec_IECFixed[1]).reshape((1,1,1))
        dz_IEC_F = np.array( dispVec_IECFixed[2]).reshape((1,1,1))
        #Convert in the IEC Gantry
        (dx,dy,dz) = mathTools.rot3D(dx_IEC_F,dy_IEC_F,dz_IEC_F,np.identity(3),self._rotMatrix,self._globalVar.typeFloat)
        dispVec_IECG = np.zeros((3),dtype=self._globalVar.typeFloat)
        dispVec_IECG[0] = dx[0,0,0]
        dispVec_IECG[1] = dy[0,0,0]
        dispVec_IECG[2] = dz[0,0,0]
        #Convert in the CT system
        isocenterPos = (0,0,0) # : we just need the displacement vector, centered on zero, i.e. not at the isocenter
        (X_ct,Y_ct,Z_ct) = mathTools.fromIECFixedToCTCoord(dx_IEC_F,dy_IEC_F,dz_IEC_F,isocenterPos)
        dispVec_CT = np.zeros((3),dtype=self._globalVar.typeFloat)
        dispVec_CT[0] = X_ct[0,0,0]
        dispVec_CT[1] = Y_ct[0,0,0]
        dispVec_CT[2] = Z_ct[0,0,0]        
        #Convert CT coordinates into matrix indice
        indDispCT = self._patient.getDispIndices(dispVec_CT)
        
        print "\t***** Displacement info %s *****"%str(info)
        print "Displacement in IEC fixed (x,y,z) (cm): (%f,%f,%f)"%(dispVec_IECFixed[0],dispVec_IECFixed[1],dispVec_IECFixed[2])
        print "Displacement in IEC gantry (x,y,z) (cm): (%f,%f,%f)"%(dispVec_IECG[0],dispVec_IECG[1],dispVec_IECG[2])
        print "Displacement in CT (x,y,z) (cm): (%f,%f,%f)"%(dispVec_CT[0],dispVec_CT[1],dispVec_CT[2])
        print "Displacement in CT matrix (f,r,c): (%i,%i,%i)"%(indDispCT[0],indDispCT[1],indDispCT[2])
        print "\t*****************************"
        return (dispVec_IECFixed ,dispVec_IECG,dispVec_CT,indDispCT )
        
        
#############################################################################
###         Functions for compensation methods
#############################################################################           
    def storeDictStatisticsRepaint(self,dictRep):
        '''Store statistics about the repainting: the number of compensated beam positions for each repainting index.\
        It is stored to a csv file in the simulation directory in 'StatisticsRepaintings/'.
        
        :param dictRep: Dictionary storing for each repainting index the number of compensated positions: \
        dictRep[gantry angle][energy][repainting Index].
        
        '''
        
        pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation +'StatisticsRepaintings/'
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)
        for angle in dictRep.keys():
            for ener in dictRep[angle].keys():
                enerData = list()
                enerData.append(['RepIdx','NbComp'])
                for repIndex in dictRep[angle][ener].keys():
                    enerData.append([repIndex,dictRep[angle][ener][repIndex]])
                array = np.array(enerData)
                filename = pathToSave+"StatRepaint-Angle%0.3d-%0.3dMeV.csv"%(int(angle),int(ener))
                np.savetxt(filename, array, delimiter=",", fmt="%s")
                
    def compensationFindStartSpotEnerLayer(self, scanPath , ddCurvesData , strEner):
        '''If MSPT global variables compensationDynamic and findStartPos are True, it look for a 'good starting position'\
        for delivery of the current the energy layer. See dicomReader.scanningPath method findStartingPointForMotion() for more information.
        
        :param scanPath: scanning path object (dicomReader.scanningPath.ScanningPathSingleBeam)
        :param ddCurvesData: depth dose curves data : physics.ddcurveManager.DDCurves
        :param strEner: String indicating if the energy being used has just been set ('NewEnergy') or has already been set for at least one \
        beam position ('SameEnergy').
            
        :returns: None if no starting position has been found, otherwise a tuple:
            
            * index 0: spot weight
            * index 1: spot in the 3D matrices : frame, row, column 
            * index 2: x position
            * index 3: y position
            * index 4: energy
            * index 5: list of 2 elements : spot size along x and along y 
            * index 6: 'NewEnergy' if the energy is changing or 'SameEnergy' otherwise 
            * index 7: spot type = 1 (native spot) or 2 (spot from added margin) 
            * index 8: repainting index 
            
        '''
        (dispVec_IECFixed ,dispVec_IECG,dispVec_CT,indDispCT) = self.computeDisplacementVectors()
        self.updateTime('Displacement Measure')
        return scanPath.findStartingPointForMotion( dispVec_IECG * 10.0 ,ddCurvesData, strEner) # dispVec_IECG * 10.0 : convert from cm to mm

    
#############################################################################
###         Functions shared by static and dynamic
#############################################################################    
                 
    def findBraggPeak(self, compDose, deff):
        '''Finds and display the Bragg Peak depth in the beam dose distribution. 
    
        :param compDose: Computed 3D dose distribution
        :param deff: 3D numpy array of the radiological depth
    
        '''
        #We assume beam along z dcm axis
        indBraggPeak = np.unravel_index(compDose.argmax(), compDose.shape)
#         print "ind Bragg Peak: %s"%(str(indBraggPeak))
        print "Bragg Peak at :%f in patient - spacing = %f "%(float(self._ySpacing * (compDose.shape[1] - indBraggPeak[1])) ,self._ySpacing)
        #print "Bragg Peak at :%f in water"%(deff[indBraggPeak])
        

def beamletSimul(rhomw,Smw,Xr,Yr,Zr,dmax,z_dd,dd,sig0,xrshift,yrshift,sad,ssd0,Deff, indicesRay,x0y0z0_IECF,IECFSpacing,rotMatrixInv,typeFloat):
    '''Calls the C extension: beamletSimulation_PythonExt and returns its result.
        
    The goal of this function is to simulate a beamlet and get the dose deposited by this beamlet.
    
    :param rhomw: 3D numpy array density ratio to water (i.e. density matrix since water density = 1.0 g.cm-3)
    :param Smw: 3D numpy array of the relative stopping power ratio to water
    :param Xr: 3D matrix representing x coordinates in the IEC Gantry coordinate system (cm).
    :param Yr: 3D matrix representing y coordinates in the IEC Gantry coordinate system (cm).
    :param Zr: 3D matrix representing z coordinates in the IEC Gantry coordinate system (cm).
    :param dmax: Bragg Peak depth (cm).
    :param z_dd: "depth" coordinates for the depth dose curve (cm).
    :param dd: "dose" coordinates for the depth dose curve. (Gy.cm ^2 )
    :param sig0: beam size (along x and y IEC gantry axis) at patient entrance (cm): [sig0X,sig0Y]
    :param xrshift: x beam shift (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
    :param yrshift: y beam shift (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
    :param sad: source-axis distance: distance between the source and the isocenter.(cm)
    :param ssd0: source to patient surface distance (cm)
    :param Deff: 3D matrix representing the radiological depth.(cm)
    :param indicesRay: list of indices (frame,row, column) of the voxels encountered by the ray from the source to the target spot (xrshift,yrshift)\
    located in the isocenter plane. They are sorted by increasing distance from the source.
    :param x0y0z0_IECF: beam target point in the IEC fixed coordinate system.
    :param spacingIECF: voxel spacing in the IEC fixed coord system
    :param rotMatrixInv: 3x3 array (rotation matrix) to convert coordinates from IECG to IECF
    :param typeFloat: the type of numpy float to be used: 'float32' or 'float64'

    :returns: The resulting 3D dose distribution in a numpy array. 
 
    '''
    test=np.zeros((2,2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    for idx,name in enumerate(['rhomw','Smw','Xr','Yr','Zr','Deff']):
        mat = eval(name)
        if type(mat) != typetest:
            print "type '%s' : %s"%(name,type(mat))
            raise Exception('In beamletSimulation, mat is not *NumPy* array')
        if len(np.shape(mat)) != 3:
            print "len shape '%s' : %i"%(name,len(np.shape(mat)))
            raise Exception('In beamletSimulation, mat is not NumPy *3D array*')
        if mat.dtype != datatest:
            print "type '%s' : %s"%(name,mat.dtype)
            raise Exception('In beamletSimulation, mat is not *Float* NumPy array')
    
    test=np.zeros((2,2),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    for idx,name in enumerate(['z_dd','dd']):
        mat = eval(name)
        if type(mat) != typetest:
            print "type '%s' : %s"%(name,type(mat))
            raise Exception('In beamletSimulation, mat is not *NumPy* array')
        if len(np.shape(mat)) != 1:
            print "len shape '%s' : %i"%(name,len(np.shape(mat)))
            raise Exception('In beamletSimulation, mat is not NumPy *3D array*')
        if mat.dtype != datatest:
            print "type '%s' : %s"%(name,mat.dtype)
            raise Exception('In beamletSimulation, mat is not *Float* NumPy array')
        
    test=np.zeros((1),dtype=typeFloat)
    typetest= type(test)
    datatest=test.dtype
    for idx,name in enumerate(['dmax','xrshift','yrshift','sad','ssd0']):
        mat = eval(name)
        if type(mat) != typetest:
            print "type '%s' : %s"%(name,type(mat))
            raise Exception('In beamletSimulation, mat is not *NumPy* array')
        if mat.size != 1:
            print "size '%s' : %i"%(name,mat.size)
            raise Exception('In beamletSimulation, mat is not NumPy *3D array*')
        if mat.dtype != datatest:
            print "type '%s' : %s"%(name,mat.dtype)
            raise Exception('In beamletSimulation, mat is not *Float* NumPy array')

    test=np.zeros((1),dtype='int')
    typetest= type(test)
    datatest=test.dtype
    for idx,name in enumerate(['indicesRay']):
        mat = eval(name)
        if type(mat) != typetest:
            print "type '%s' : %s"%(name,type(mat))
            raise Exception('In beamletSimulation, mat is not *NumPy* array')

        if mat.dtype != datatest:
            print "type '%s' : %s"%(name,mat.dtype)
            raise Exception('In beamletSimulation, mat is not *int* NumPy array')
            
    for idx,name in enumerate(['x0y0z0_IECF','IECFSpacing']):
        item = eval(name)
        test=np.zeros((3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
 
 
    for idx,name in enumerate(['sig0']):
        item = eval(name)
        test=np.zeros((2),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 1:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)
 
 
    for idx,name in enumerate(['rotMatrixInv']):
        item = eval(name)
        test=np.zeros((3,3),dtype=typeFloat)
        typetest= type(test)
        datatest=test.dtype
        
        if type(item) != typetest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *NumPy* array'
            raise Exception(stringError)
        if len(np.shape(item)) != 2:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not 2D NumPy *matrix*'
            raise Exception(stringError)
        if item.dtype != datatest:
            print item
            stringError = 'In beamletSimulation, '+str(name) +' is not *Float* NumPy array'
            raise Exception(stringError)


    for idx,name in enumerate(['rhomw','Smw','Xr','Yr','Zr','dmax','z_dd','dd','sig0','xrshift','yrshift','sad','ssd0','Deff','indicesRay','x0y0z0_IECF','IECFSpacing','rotMatrixInv']):
        item = eval(name)
        if not item.flags['C_CONTIGUOUS']:
            strErr = 'In beamletSimulation, '+str(name) +' is not *C-Contiguous* NumPy array'
            raise ValueError(strErr)


    if typeFloat == 'float32':
        return _beamletSimulation.beamletSimul(rhomw,Smw,Xr,Yr,Zr,dmax,z_dd,dd,sig0,xrshift,yrshift,sad,ssd0,Deff,indicesRay,x0y0z0_IECF,IECFSpacing,rotMatrixInv)
    else:
        return _beamletSimulationDouble.beamletSimulDouble(rhomw,Smw,Xr,Yr,Zr,dmax,z_dd,dd,sig0,xrshift,yrshift,sad,ssd0,Deff,indicesRay,x0y0z0_IECF,IECFSpacing,rotMatrixInv)


    