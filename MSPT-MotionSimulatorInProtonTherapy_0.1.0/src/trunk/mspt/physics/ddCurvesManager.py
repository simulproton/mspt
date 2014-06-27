########################################################################
#
# ddCurvesData.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# Jan 13 2014
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
import numpy as np
import os,sys
import fnmatch
from scipy.interpolate import interp1d
from ..mathTools import tools as mathTool

patternDDCurvesFilenames = 'depthDoseData-???_??MeV*.csv'
patternEnerDigit0 = 14
patternEnerDigit2 = 16
patternEnerDigit3 = 18
patternEnerDigit4 = 19

patternDoseCorrecFile = 'doseCorrecFactor_MSPT.txt'


class DDCurves(object):
    '''Class managing the depth dose curves used in MSPT.
    
    :param configData: dictionary (globalVariables) storing the data configuring the simulator. \
    Mandatory keys for the depth dose curves are:
        
        * 'ddCurves' : The type of Depth-Dose curves to use. In MSPT only RayStation DD curves are available.
        * 'typeFloat': type of numpy float to use: 'float32' or 'float64'.
        * 'nameDoseCorrecFactor': name of the dose correction factor to use.
    
    .. note::
    
        If one wants to add new depth-dose curves, please follow the procedure:
        
            #. Create a new directory *"NAME_NEWFOLDER"* in the folder:
    
                *MSPT-MotionSimulatorInProtonTherapy/src/trunk/RefData/DepthDoseData/*
                
            #.  Add the depth dose curves in this new directory. Please, note that the depth-dose curves must be stored in .csv files:
                
                #. 1 file per energy
                #.  The name of each file must respect this pattern: 'depthDoseData-xxx_xxMeV*.csv' where the x corresponds to digits of the \
                energy and '*' can be any string. For example: depth-dose curve files for 10.20MeV: "depthDoseData-010_20MeV.csv" \
                or  "depthDoseData-010_20MeV_SomeInformation.csv" or "depthDoseData-010_20MeV_NAME_NEWFOLDER.csv"
                #. Each file must be organized as follows: 1st row: the depths in cm, 2nd row the corresponding dose in Gy.cm2.
            #. To use these depth-dose curves set the variable *ddCurves*, of the configuration file, to 'NAME_NEWFOLDER'.
            
    
    
    .. note::
    
        If one wants to modify the dose corrections factors:
        
            #. Create a table with 2 columns into a tab delimited '.txt' file:
                
                #. Then fist line must contains the headers: 1st column: 'Energy', 2nd column: 'Factor'
                #. Then, for each new line:  1st column: the energy in MeV , 2nd column: the corresponding correction factor unitless.
                #. The filename must respect this pattern: 'doseCorrecFactor__NAMEDATA.txt' where 'NAMEDATA' will be the name \
                you want to give to this data. For example, by default the name for the correction factors used in MSPT is 'MSPT', so the file name is 'doseCorrecFactor_MSPT.txt'
            
            #. Place the file in the folder *MSPT-MotionSimulatorInProtonTherapy/src/trunk/RefData/DoseCorrecFactor/*
            #. Set the MSPT configuration variable *nameDoseCorrecFactor* to 'NAMEDATA'
            
        And you're done!
    
    '''

    def __init__(self,configData):
        self._globalVar = configData
        if self._globalVar.ddCurves == 'RayStation': 
            self._depthDoseDataDirectory = (os.getcwd()+'/RefData/DepthDoseData/RayStation/')
#         elif self._globalVar.ddCurves == 'MCNPX':
#             self._depthDoseDataDirectory = (os.getcwd()+'/RefData/DepthDoseData/MCNPX/')
#         elif self._globalVar.ddCurves == 'Eclipse': 
#             self._depthDoseDataDirectory = (os.getcwd()+'/RefData/DepthDoseData/Eclipse/')
        else:
            if os.path.exists(os.getcwd()+'/RefData/DepthDoseData/%s/'%self._globalVar.ddCurves):
                self._depthDoseDataDirectory = os.getcwd()+'/RefData/DepthDoseData/%s/'%self._globalVar.ddCurves
            else:
                strErr = 'Unknown ddCurve type %s'%(self._globalVar.ddCurves)
                raise ValueError(strErr)
            
        self._DDDataFilename = self._depthDoseDataDirectory
        self._lookupTable = self.readCSVFilesInFolder(self._DDDataFilename)

        self._luTableMissingEnergies = None  
        self._funcBPvsEnergy = None #Function BP = f(energy)
        self._funcEnergyvsBP = None #Function Energy = f(BP)
        self._funcDoseAtBPVsEnergy = None # Function DoseAtBP = f(Energy)
        self._listEnergies = np.array(self._lookupTable.keys())
        self.loadDoseCorrectionFactor()
        
        
        
    
    def dataForEnergy(self, dataType , energy):
        '''Returns the data for a given energy. Data types are:
            
            * 'DD_Curve': the depth dose curve stored in an array (2 rows, n columns).
            * 'BraggPeak': Depth of the Bragg Peak (BP) in cm.
            * 'MaxDose': Maximum dose value (at Bragg Peak).
            * 'CumulDose': Array of the cumulative dose (2 rows, n columns).
            * 'NbProtons': Number of protons estimated from the depth dose curve.
            * 'Ratio-DoseEntrance-BP': Ratio between the entrance dose and the dose at the BP.
            * 'DoseAtEntrance': Dose at patient entrance
            * 'EnergyLoss': Array of the energy loss along the depth (2 rows, n columns). 
            * 'EnergyRemaining': Array of the remaining energy along the depth (2 rows, n columns).        
        
        :param dataType: string specifying the dataType
        :param energy: energy 
        
        :returns: The data.
        
        '''
        maxEner = max(self._lookupTable.keys())
        minEnergy =  min(self._lookupTable.keys())
        
        energy = float("%0.3f"%energy)
        if 0 < (energy - max(self._lookupTable.keys())) < 1:
            energy = float(max(self._lookupTable.keys()))
#         print "Getting data for energy: %s"%(energy)
        diffEner = np.abs(self._listEnergies - energy)
        nearestEner = np.where(diffEner <= 0.2)
        if len(nearestEner[0]) != 0:
            enerToUse = self._lookupTable.keys()[nearestEner[0][0]]
            return self._lookupTable[ enerToUse][dataType]
        elif self._luTableMissingEnergies is None or (self._luTableMissingEnergies is not None \
                and energy not in self._luTableMissingEnergies.keys()) :
            if energy <= max(self._lookupTable.keys()):
                self.generateApproxDDCurveForEnergy(energy)
                return self._luTableMissingEnergies[energy][dataType]
            else:
                
                if minEnergy <= energy and (energy - maxEner) < 6:
                    return self._lookupTable[maxEner][dataType]
                elif energy <= max(self._lookupTable.keys()) and (minEnergy - energy) < 6:
                    return self._lookupTable[minEnergy][dataType]
                else:
                    strErr = "Energy required : %i MeV"%energy + "\n" + "Depth Dose curves can be generated only if energy is lower than %i en greater than %i. Current energy: %s"%(int(min(self._lookupTable.keys())),int(max(self._lookupTable.keys())) ,str(energy))
                    raise ValueError (strErr)
            
        else :
            return self._luTableMissingEnergies[energy][dataType]
 
        




    def generateApproxDDCurveForEnergy(self,energy):
        '''Generate Depth Dose Curve for a given energy if this energy is outside the lookup table energy range
        
        :param energy: Energy for which to generate a depth-dose curve.
        
        The depth dose curves is generated using a linear interpolation between encompassing energies to obtain the Bragg Peak (BP) location\
        Then it uses the DD curves of the closest and higher energy available E2: the depths position are shifted by the difference \
        between the BP of the generated curve and the BP of E2. The part of the curve is kept where the depth is positive.
        
        '''
        print "Generating DD Curve %f MeV"%energy
        listEnergies = self._lookupTable.keys()
        listEnergies.sort()
        
        prevEner = None
        nextEner = None
        
        if energy < listEnergies[0]:
            dataEner0 = self._lookupTable[listEnergies[0]]['DD_Curve']
            idxBP0 =  dataEner0[1].argmax()
            bp0 = dataEner0[0,idxBP0]
            dataEner1 = self._lookupTable[listEnergies[1]]['DD_Curve']
            idxBP1 =  dataEner1[1].argmax()
            bp1 = dataEner1[0,idxBP1]   
        
            newBP = mathTool.interp(energy, listEnergies[0],listEnergies[1],bp0,bp1)
            deltaZ = bp0 - newBP

            
            doseAtBP0 = self.findDoseAtBPForEner(listEnergies[0])
            doseAtBP1 = self.findDoseAtBPForEner(listEnergies[1])
            doseAtBP = mathTool.interp(energy, listEnergies[0],listEnergies[1],doseAtBP0,doseAtBP1)

            ddCurve = np.zeros(dataEner0.shape)
            
            ddCurve[0,:]= dataEner0[0,:] - deltaZ
            positivZ = np.where(ddCurve[0,:] > 0)
            
            waterDens = 1e-3 #kg/cm3
            MeV = 1.602e-13
            nbProt = self._lookupTable[listEnergies[0]]['NbProtons']
            removedArea = 0
            startz = 0
            for idxCol in range(positivZ[0][0]):
                removedArea = removedArea + dataEner0[1,idxCol] * np.abs(startz - dataEner0[0,idxCol])
                startz = dataEner0[0,idxCol]
            factorDose = ( nbProt / waterDens) * energy * MeV / np.abs(self._lookupTable[listEnergies[0]]['CumulDose'][-1] - removedArea)
        
            ddCurve[1,:] = dataEner0[1,:] *  factorDose
        
            newDDCurve = np.zeros((2, len(positivZ[0])))
            newDDCurve[0,:] = ddCurve[0,positivZ]
            newDDCurve[1,:] = ddCurve[1,positivZ]
            
        else:
            for ener0,ener1 in zip(listEnergies[0:-1],listEnergies[1:]):
                if ener0 < energy < ener1:
                    prevEner = ener0
                    nextEner = ener1
                    break
        
            dataEner0 = self._lookupTable[prevEner]['DD_Curve']
            dataEner1 = self._lookupTable[nextEner]['DD_Curve']
        
            #Bragg Peak
            idxBP0 =  dataEner0[1].argmax()
            bp0 = dataEner0[0,idxBP0]
            idxBP1 =  dataEner1[1].argmax()
            bp1 = dataEner1[0,idxBP1]
        
            print "Ener0: %f , BP: %f"%(prevEner,bp0)
            print "Ener1: %f , BP: %f"%(nextEner,bp1)
        
            if self._funcBPvsEnergy is None:
                self.buildFuncBPvsEner()
        
            newBP = mathTool.interp(energy, prevEner,nextEner,bp0,bp1)
            print "Ener: %f , BP: %f"%(energy,newBP)
            print "test 0: %f"%self._funcBPvsEnergy(prevEner)
            print "test 1: %f"%self._funcBPvsEnergy(nextEner)
        
            deltaZ = bp1 - newBP
        
            doseAtBP = self.findDoseAtBPForEner(energy)
            doseAtBP1 = self.findDoseAtBPForEner(nextEner)
        
        
            ddCurve = np.zeros(dataEner1.shape)
            ddCurve[0,:]=  dataEner1[0,:] - deltaZ
            positivZ = np.where(ddCurve[0,:] > 0)
            
            waterDens = 1e-3 #kg/cm3
            MeV = 1.602e-13
            nbProt = self._lookupTable[nextEner]['NbProtons']
            removedArea = 0
            startz = 0
            for idxCol in range(positivZ[0][0]):
                removedArea = removedArea + dataEner1[1,idxCol] * np.abs(startz - dataEner1[0,idxCol])
                startz = dataEner1[0,idxCol]
            
            
            factorDose = (nbProt/waterDens)*energy * MeV / np.abs(self._lookupTable[nextEner]['CumulDose'][-1] - removedArea)
            ddCurve[1,:] = dataEner1[1,:] *  factorDose
              
            newDDCurve = np.zeros((2, len(positivZ[0])))
            newDDCurve[0,:] = ddCurve[0,positivZ]
            newDDCurve[1,:] = ddCurve[1,positivZ]
            idxBPNew =  newDDCurve[1].argmax()
            bpNew = newDDCurve[0,idxBPNew]
            
        pathToSave = self._DDDataFilename
        pathToSave = pathToSave + "GeneratedDDCurves/"
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)  
        name = "depthDoseData-%0.3d_%0.2dMeV"%(int(energy),int((energy - int(energy))*100))
        filename = pathToSave + name +'.csv'           
        np.savetxt( filename,newDDCurve,delimiter=',')            
        print "Saved %s"%filename
        self.generateTableMissingEnergy(energy,newDDCurve)

                       
            
    def generateTableMissingEnergy(self,energy,newDDCurve):
        '''Function that generates the dictionary in which is stored all the information about the depth dose curve.
        
        :param energy: Energy in MeV
        :param newDDCurve: depth-dose curve generated: 1st row depths in cm, second row dose.
        
        The new dictionary created is added to the look-up table storing the missing energies. The keys of this dictionary are:
        
            * *'DD_Curve'* : the depth dose curve
            * *'BraggPeak'* : Brag Peak depth in cm
            * *'MaxDose'* : Maximum dose in the DD curve.
            * *'CumulDose'* : Cumulative dose along the depth
            * *'NbProtons'* : Number of protons used in the beam corresponding to this DD Curve
            * *'Ratio-DoseEntrance-BP'* : Ratio between the dose at initial depth and at BP depth.
            * *'DoseAtEntrance'* : Dose at initial depth
            * *'EnergyLoss'* : Energy loss along the depth
            * *'EnergyRemaining'* : Remaining energy along the depth
        
        '''
        if self._luTableMissingEnergies is None:
            self._luTableMissingEnergies = dict()
        if energy in self._luTableMissingEnergies.keys():
            print "Energy %i MeV already in table."%energy
            return
        else:
            self._luTableMissingEnergies[energy] = dict()
            self._luTableMissingEnergies[energy]['DD_Curve'] =  newDDCurve
            bpData = findBraggPeak(newDDCurve)
            self._luTableMissingEnergies[energy]['BraggPeak'] = bpData[0]
            self._luTableMissingEnergies[energy]['MaxDose'] = bpData[1]
            cumDose = cumulDose(newDDCurve)
            self._luTableMissingEnergies[energy]['CumulDose'] = cumDose
            self._luTableMissingEnergies[energy]['NbProtons'] = computeNbProtonsMCSimul( cumDose[-1],energy)
            self._luTableMissingEnergies[energy]['Ratio-DoseEntrance-BP'] =newDDCurve[1,0] / bpData[1]
            self._luTableMissingEnergies[energy]['DoseAtEntrance'] = newDDCurve[1,0]
            self._luTableMissingEnergies[energy]['EnergyLoss'] = computeEnergyLoss(cumDose , energy)
            self._luTableMissingEnergies[energy]['EnergyRemaining'] = computeEnergyRemaining( self._luTableMissingEnergies[energy]['EnergyLoss'], energy)
            print "Lookup table updated for energy %f:\n %s"%(energy,str(self._luTableMissingEnergies.keys()))
            
                
    
    def findEnergyForBraggPeak(self, desiredBP):
        '''Find an energy for a given Bragg Peak
        
        :param desiredBP: Desired BP depth
        
        :returns: Corresponding energy estimated using a linear interpolation.
        
        '''
        if self._funcEnergyvsBP is None:
            self.buildFuncEnervsBP()
        try:
            ener = self._funcEnergyvsBP(desiredBP)
        except:
            print "findEnergyForBraggPeak: energy outside the possible range..."
            return None
        return ener
 

    def findDoseAtBPForEner(self,energy):
        '''Obtain the dose at BP for a given energy
        
        :param energy: Desired energy in MeV
        
        :returns: Corresponding BP depth
        
        '''
        if self._funcDoseAtBPVsEnergy is None:
            self.buildFuncDoseAtBPVsEnergy()
        try:
            dose = self._funcDoseAtBPVsEnergy(energy)
        except:
            print "findDoseAtBPForEner: energy outside the possible range..."
            return None
        return dose        

    def buildFuncBPvsEner(self):
        '''Creates an interpolation function to provide the BP for a given energy
        
        '''
        dataEnergies = np.array(self._lookupTable.keys())
        dataEnergies = np.sort(dataEnergies)
        dataBP = []
        for ener in dataEnergies:
             dataBP.append(self.dataForEnergy('BraggPeak' , ener))
        
        dataBP = np.array(dataBP)
        self._funcBPvsEnergy = interp1d(dataEnergies, dataBP, kind = 'linear') #kind='linear')   kind = 'quadratic') #kind='linear')


    def buildFuncEnervsBP(self):
        '''Creates an interpolation function to provide the energy for a given BP
        
        '''
        dataEnergies = np.array(self._lookupTable.keys())
        dataEnergies = np.sort(dataEnergies)
        dataBP = []
        for ener in dataEnergies:
             dataBP.append(self.dataForEnergy('BraggPeak' , ener))
        
        dataBP = np.array(dataBP)
        self._funcEnergyvsBP  = interp1d(dataBP, dataEnergies,  kind='linear')


    def buildFuncDoseAtBPVsEnergy(self):   
        '''Creates the interpolation function providing the dose at BP for a specific energy
        
        '''
        dataEnergies = np.array(self._lookupTable.keys())
        dataEnergies = np.sort(dataEnergies)
        dataDoseBP = []
        for ener in dataEnergies:
             dataDoseBP.append(self.dataForEnergy('MaxDose' , ener))
        
        dataDoseBP = np.array(dataDoseBP)
        self._funcDoseAtBPVsEnergy = interp1d(dataEnergies, dataDoseBP, kind = 'linear') #kind='linear')   kind = 'quadratic') #kind='linear')
        

    def readFile(self,filename):
        '''Load depth dose data and computes additional data (data at BP, cumulative dose,number of protons, ratio dose entrance/BP\
        , energy loss, energy remaining) and returns it.
        
        :param filename: Path to file to load
        
        :returns: A list with data:
            
            #. energy
            #. DD curve stored in a (2,n) array
            #. BP data : depth and dose at the BP
            #. Cumulative dose
            #. Number of protons
            #. Ratio dose
            #. Energy loss
            #. Remaining energy
        
        '''
        name = os.path.split(filename)[1]
        if fnmatch.fnmatch(name, patternDDCurvesFilenames): 
            dataDDCurve = np.loadtxt(filename,dtype=self._globalVar.typeFloat,delimiter=',')
            subFilename = os.path.basename(filename)
            energy = np.float('%i.%i%i'%(int(name[patternEnerDigit0:patternEnerDigit2+1]),\
                                            int(name[patternEnerDigit3]),\
                                            int(name[patternEnerDigit4])))

            bpData = findBraggPeak(dataDDCurve)
            cumDose = cumulDose(dataDDCurve)
            nbProt = computeNbProtonsMCSimul(cumDose[-1],energy)
            ratioDose =  dataDDCurve[1,0] / bpData[1]
            enerLoss = computeEnergyLoss(cumDose , energy)
            enerRemain = computeEnergyRemaining( enerLoss, energy)
    
            data = [energy , dataDDCurve ,bpData , cumDose, nbProt, ratioDose, enerLoss,enerRemain]
            return data
        else:  
            print "Does not mach pattern"


    
    def readCSVFilesInFolder(self, dirPath):
        '''Function that loads all the Depth-Dose curves contained in a directory
        
        :param dirPath: Path to the directory
        
        The function creates a lookup table (dictionary) whose keys are the energies of the depth dose curves and the value \
        is a dictionary storing the data computed by 'readFile(filename)'. The keys for each energy:
        
            * *'DD_Curve'* : the depth dose curve
            * *'BraggPeak'* : Brag Peak depth in cm
            * *'MaxDose'* : Maximum dose in the DD curve.
            * *'CumulDose'* : Cumulative dose along the depth
            * *'NbProtons'* : Number of protons used in the beam corresponding to this DD Curve
            * *'Ratio-DoseEntrance-BP'* : Ratio between the dose at initial depth and at BP depth.
            * *'DoseAtEntrance'* : Dose at initial depth
            * *'EnergyLoss'* : Energy loss along the depth
            * *'EnergyRemaining'* : Remaining energy along the depth
        
        :returns: the look up table.
        
        '''
        lookupTable = dict()
        for r,d,f in os.walk(dirPath):
            for files in f:
                if fnmatch.fnmatch(files, patternDDCurvesFilenames):
                    filename = os.path.join(r,files)
                    data = self.readFile(filename)
                    lookupTable[data[0]] = dict()
                    lookupTable[data[0]]['DD_Curve'] = data[1]
                    lookupTable[data[0]]['BraggPeak'] = data[2][0]
                    lookupTable[data[0]]['MaxDose'] = data[2][1]
                    lookupTable[data[0]]['CumulDose'] = data[3]
                    lookupTable[data[0]]['NbProtons'] = data[4]
                    lookupTable[data[0]]['Ratio-DoseEntrance-BP'] = data[5]
                    lookupTable[data[0]]['DoseAtEntrance'] = data[1][1,0]
                    lookupTable[data[0]]['EnergyLoss'] = data[6]
                    lookupTable[data[0]]['EnergyRemaining'] = data[7]
            return lookupTable# stop after the first level of files has been processed
               

    def loadDoseCorrectionFactor(self):
    
        '''Function used to import dose correction factor data. It scans the folder */RefData/DoseCorrecFactor/*\. 
        
        Each file matching the pattern 'doseCorrecFactor_NAMEDATA.txt'.\
        NAMEDATA must also be provided by the mspt global variable : nameDoseCorrecFactor = 'NAMEDATA'
        
        '''
        dirPath = os.getcwd()+'/RefData/DoseCorrecFactor/'
        print "#### Importing dose correction factor ####"
        self._DoseCorrecFactor = None
        self._funcInterpDoseCorrec = None
        for r,d,f in os.walk(dirPath):
            for files in f:
                if fnmatch.fnmatch(files, patternDoseCorrecFile):
                    nameData = os.path.splitext(files)[0].split('doseCorrecFactor_')[1]
                    if nameData == self._globalVar.nameDoseCorrecFactor:
                        filename = os.path.join(r,files)
                        loadedArray = np.loadtxt(filename,dtype="|S",delimiter='\t',skiprows = 1,usecols=(0,1))
                        loadedArray = loadedArray.astype(self._globalVar.typeFloat)
                        loadedArray = loadedArray.transpose()
                        loadedArray.view('%s,%s'%(self._globalVar.typeFloat,self._globalVar.typeFloat)).sort(order=['f1'], axis=0)
                        print "Dose correction factor %s imported"%nameData
                        self._DoseCorrecFactor = np.array(loadedArray,dtype = self._globalVar.typeFloat, order = 'C')
            print "YOYO"
            if self._DoseCorrecFactor is None:
                print "Correction factor data (dose) for name : %s , has not been found"%self._globalVar.nameMassStopPwrCorrecFactor
                raise ValueError('Error - No dose correction factor data has been found in /RefData/DoseCorrecFactor/. Please check the data and try again..')
            print "#### End Importing dose correction factor ####"
            return # stop after the first level of files has been processed


    def doseCorrectionFactor(self, energy):
        '''Function created to match RayStation absolute dose. It uses a linear interpolation to obtain a multiplicative dose correction factor for an energy.
    
        The function has been obtained by measuring the ratio between the dose obtained in RayStation and in our simulator
        for single beamlets in water for energies ranging from 30MeV to 245MeV.
    
        :param energy: energy in MeV
    
        :returns: The correction factor
    
        '''
    
#         listCoeff = [1.926,	1.346,	1.215,	1.155,	1.117,	1.090,	1.070,	1.055,	1.044,	1.034,	1.026,	1.019,	1.014,	1.008,	1.004,	0.999,	0.995,	0.992,	0.989,	0.986,	0.983,	0.980,	0.977,	0.974,	0.972,	0.969,	0.966,	0.964,	0.961,	0.959,	0.956,	0.954,	0.951,	0.949,	0.947,	0.944,	0.942,	0.939,	0.937,	0.934,	0.932,	0.930,	0.927,	0.923,	0.920,	0.918]
#         listEner =   [20.1    ,	25.1    ,	30.1    ,	35.1    ,	40.1    ,	45.1    ,	50.1    ,	55.1    ,	60.1    ,	65.09   ,	70.09   ,	75.09   ,	80.09   ,	85.09   ,	90.09   ,	95.09   ,	100.09  ,	105.09  ,	110.09  ,	115.09  ,	120.09  ,	125.09  ,	130.09  ,	135.09  ,	140.09  ,	145.09  ,	150.09  ,	155.09  ,	160.09  ,	165.09  ,	170.09  ,	175.09  ,	180.09  ,	185.09  ,	190.09  ,	195.09  ,	200.09  ,	205.09  ,	210.09  ,	215.09  ,	220.09  ,	225.09  ,	230.09  ,	235.09  ,	240.09  ,	245.09]
        if self._funcInterpDoseCorrec is None:
            listCoeff = self._DoseCorrecFactor[1]
            listEner = self._DoseCorrecFactor[0]
            self._funcInterpDoseCorrec = interp1d(listEner, listCoeff, kind='linear',bounds_error=False,fill_value = 1.0)
        factor = self._funcInterpDoseCorrec(energy)
        return factor



def findBraggPeak( dataDDCurve):
    '''Finds the Bragg Peak (cm) and return the BP with the deposited dose associated in [Gy.cm^2]
    
    :param dataDDCurve: stored in a (2,n) array
    
    :returns: list [ BP depth (cm), dose at BP]
    
    '''
    indMaxDose = np.argmax(dataDDCurve[1,:])
    maxDose = dataDDCurve[1,indMaxDose]
    braggPeak =  dataDDCurve[0,indMaxDose]
    return [braggPeak ,maxDose]
    
def cumulDose( dataDDCurve):
    '''Computes the cumulative dose along the depth:
    
        Cumul dose = sum ( dose * dz ), where dose in [Gy.cm^2] * [cm]
    
    :param dataDDCurve: stored in a (2,n) array
    
    :returns: Cumulative dose array
    
    '''
    cumDoseArray = np.zeros((dataDDCurve.shape[1]))
    
    startDepth = 0
    cumDose = 0
    for idx in range(dataDDCurve.shape[1]):
        cumDose = cumDose + (dataDDCurve[0,idx] - startDepth) * dataDDCurve[1,idx]
        startDepth = dataDDCurve[0,idx]
        cumDoseArray[idx] = cumDose 
    return cumDoseArray
    
def computeNbProtonsMCSimul(cumDoseVal,energy):
    '''For given energy , estimate the number of protons used in the corresponding beam in water:
    
        nbProtons = sum_dose[Gy.cm^2] * z_spacing [cm] * rhoWater [kg.cm ^-3] / Energy_Joules
        
    
    .. note::
        
        Gy = J / Kg   and  1 MeV = 1.602e-13 Joules.
    
    :param cumDoseVal: cumulative dose value at the end of the DD curve range (after the BP)
    :param energy: energy in MeV
    
    :returns: The number of protons estimated.
    
    '''
    meV = 1.602e-13
    rhoWater = 1e-3 #rhoWater = 1e-3 kg/cm^3
    nbProtonsMCSimul = cumDoseVal * rhoWater / ( energy * meV)
    return nbProtonsMCSimul

    
    

def computeEnergyLoss(cumDose , energy):
    '''Computes energy loss along depth:
    
        energy loss(z) = dose(z) / cumul_dose * energy
    
    :param cumDose: array of the cumulative dose
    :param energy: Energy in MeV
    
    :returns: Array of the energy loss along the depth 
    
    '''
    energyLoss = np.zeros((cumDose.shape[0]))
    cumDoseTot = cumDose[-1]
    for idx in xrange(cumDose.shape[0]):
        energyLoss[idx] = (cumDose[idx] / cumDoseTot ) * energy 
    return energyLoss
        
def computeEnergyRemaining(energyLoss, energy):
    '''Computes remaining energy along depth:
    
        remainEnergy(z) = energy - energyLoss(z)
    
    :param energyLoss: Energy loss array
    :param energy: energy in MeV
    
    :returns: Array of the remaining energy along the depth
    
    '''
    energyRemain = np.zeros((energyLoss.shape[0]))
    for idx in xrange(energyLoss.shape[0]): 
        energyRemain[idx] = energy - energyLoss[idx]
    return energyRemain
