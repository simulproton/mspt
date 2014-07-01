########################################################################
#
#stoppingPower.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On 06/04/2012
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
##
########################################################################
import numpy as np
import sys,os,fnmatch
from ..mathTools import tools
from scipy.interpolate import interp1d
import bisect
globallistMedia = ('Air' ,'Adipose','Water','Muscle','Bone','Titanium') #Has to be sorted in increasing relative electron density order 
globalmassDensities = (0 ,0.001205,0.920,1.0,1.04,1.85,4.54)# relative mass densities from PSTAR. 
                                                            #Value 0 is added and does not correspond to a medium 
                                                            #but is used as a "minimum" value. It should be in increasing order

globalListFilenamesMassStopWater= ['/RefData/MassStoppingPower_PSTAR/MassStopPwrData-Adipose.txt','/RefData/MassStoppingPower_PSTAR/MassStopPwrData-Air.txt','/RefData/MassStoppingPower_PSTAR/MassStopPwrData-CorticalBone.txt','/RefData/MassStoppingPower_PSTAR/MassStopPwrData-Water.txt','/RefData/MassStoppingPower_PSTAR/MassStopPwrData-MuscleSkeletal.txt','/RefData/MassStoppingPower_PSTAR/MassStopPwrData-Titanium.txt']

patternFilenames = 'MassStopPwrData-*.txt'

pathToCorrecFactor = '/RefData/CorrecFactor_StopPwr/'
patternFilenameCorrecFactor = 'dataCorrecFactorStopPwr_*.txt'

class StoppingPower(object):
    '''
    The stopping power data has been extracted from the `PSTAR database <http://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html>`_ 
    
    :param configData: dictionary (globalVariables) storing the data configuring the simulator. \
    Mandatory key:
        
        * 'typeFloat': type of numpy float to use: 'float32' or 'float64'.
        * 'nameStopPwrCorrecFactor' : name of the correction factor for the mass stopping power to use.
        * 'importNewMassStopPwr' : True to use new stopping power data, false otherwise. Default is False
    
    .. note::
    
        MSPT is provided with Mass Stopping Power tables obtained from the PSTAR database for several media: air, adipose, water, muscle and bone. \
        However, if one wants to add mass stopping power data here is the procedure:
        
            #. Create the table (Energy, Mass stopping power) into a tab delimited '.txt' file (1 file per medium):
            
                #. The 1st line must contain the headers: 1st column: 'Kinetic Energy MeV' , 2nd column: 'Total Stp. Pow. MeV cm2/g'
                #. Then, for each new line:  1st column: the energy in MeV , 2nd column: the corresponding mass stopping power in MeV.cm2/g. \
                Please note that the energies must be stored in increasing order.
                #. Once all the energies and mass stopping powers have been entered, add a new line where: 1st column: medium density in g/cm3 and 2nd column: -1.0.
                #. The name of each file must respect the pattern: 'MassStopPwrData-MediumName.txt'
                
            #. Once all the files have been created (1 per medium!), add these files into the folder *MSPT-MotionSimulatorInProtonTherapy/src/trunk/RefData/MassStoppingPower_NewData/*
            #. Set the configuration variable *importNewMassStopPwr* to True to be able to use this new mass stopping power data.
    
   
    .. note::
    
        To use another correction factor for the stopping power:
        
            #. Create a table with 2 columns into a tab delimited '.txt' file:
                
                #. The fist line must contains the headers: 1st column: 'Density', 2nd column: Correction factor
                #. Then, for each new line:  1st column: the density in g/cm3 , 2nd column: the corresponding correction factor unitless.
                #. The filename must respect this pattern: 'dataCorrecFactorStopPwr_NAMEDATA.txt' where 'NAMEDATA' will be the name \
                you want to give to this data. For example, by default the name for the correction factors used in MSPT is 'MSPT', so the file name is 'dataCorrecFactorStopPwr_MSPT.txt'
            
            #. Place the file in the folder *MSPT-MotionSimulatorInProtonTherapy/src/trunk/RefData/CorrecFactor_StopPwr/*
            #. Set the MSPT configuration variable *nameStopPwrCorrecFactor* to 'NAMEDATA'
            
        And you're done!
            
    '''


    
    def __init__(self, configData):
        '''
        
        '''
        self._globalVar = configData
        self._minEner = None
        self._maxEner = None
        self._importedDensities = dict() 
        self.loadMassStoppingPowerFromPSTARData()
        if self._globalVar.importNewMassStopPwr:
            self.importNewStopPwrData()
        self.precomputeMassStopPowerRatios()
        self.importCorrecStopPwrData()
        


        
    def loadMassStoppingPowerFromPSTARData(self):
        '''Loads Mass Stopping Power data from the PSTAR data
        
        '''
        self._massStopPwrAndRangeData = dict()
        for filename in globalListFilenamesMassStopWater:
            if filename.find('Adipose') > 0 :
                key = 'Adipose'
            elif filename.find('Air') > 0 :
                key = 'Air'
            elif filename.find('Muscle') > 0 :
                key = 'Muscle'
            elif filename.find('Bone') > 0 :
                key = 'Bone'
            elif filename.find('Water') > 0 :
                key = 'Water'
            elif filename.find('Titanium') > 0:
                key = 'Titanium'
            else:
                raise ValueError('loadMassStoppingPowerFromPSTARData unrecognized filename.')
            filename = os.getcwd()+filename
            loadedArray = np.loadtxt(filename,dtype="|S",delimiter='\t',skiprows = 3,usecols=(0,3,5))
            
            self._massStopPwrAndRangeData[key] = loadedArray[84:109].astype(self._globalVar.typeFloat) #select only rows for energies from 20MeV to 275MeV
            if self._minEner is None or self._minEner < self._massStopPwrAndRangeData[key][0,0]:
                self._minEner = self._massStopPwrAndRangeData[key][0,0]
            if self._maxEner is None or self._maxEner > self._massStopPwrAndRangeData[key][-1,0]:
                self._maxEner = self._massStopPwrAndRangeData[key][-1,0]
 
    def loadMassStoppingPowerFromNewData(self,filename,medium):
        '''Loads Mass Stopping Power data from new data
        
        :param filename: File where the data is stored. The file should be using "tab delimited" format. \
        The first row should be the headers : 'Kinetic Energy MeV' and 'Total Stp. Pow. MeV cm2/g'. The \
        last row should be: 1st column: the density value in g/cm3 , 2nd column: -1.
        :param medium: Name of the medium
        
        
        '''

        loadedArray = np.loadtxt(filename,dtype="|S",delimiter='\t',skiprows = 1,usecols=(0,1))
        
        if (float(loadedArray[-1,1]) != -1):
            print "Warning -- Mass Stopping Power for medium %s was not imported:"%medium
            print "\t\t The last row of the file should be: density \t-1"
            print "\t\t Instead it is: %s"%str(floadedArray[-1])
            return
        if (float(loadedArray[-1,0]) < 0):
            print "Warning -- Mass Stopping Power for medium %s was not imported:"%medium
            print "Density : %0.3f < 0!"%float(loadedArray[-1,0])
            return
        for idx,(ener1,ener2) in enumerate( zip(loadedArray[:-2,0],loadedArray[1:-1,0])):
            if float(ener2) < float(ener1):
                print "Warning -- Mass Stopping Power for medium %s was not imported:"%medium
                print "Energies are not in increasing order : lines %i and %i : %0.3f < %0.3f!"%(idx+1,idx,float(ener2),float(ener1))
                return
            
        
        self._importedDensities[medium] = loadedArray[-1,0]
        loadedArray = np.delete(loadedArray,-1,0)
        
        if (float(loadedArray[0,0]) > self._minEner):
            newRow = np.array([self._minEner,loadedArray[0,1]])
            loadedArray = np.vstack( (newRow,loadedArray))
        elif (float(loadedArray[-1,0]) < self._maxEner):
            newRow = np.array([self._maxEner,loadedArray[-1,1]])
            loadedArray = np.vstack( (loadedArray,newRow))
            
        
        self._massStopPwrAndRangeData[medium] = loadedArray.astype(self._globalVar.typeFloat) #select only rows for energies from 20MeV to 275MeV
        
            

 
    def importNewStopPwrData(self):
        '''Function used to import new mass stopping power data. It scans the folder */RefData/MassStoppingPower_NewData/*\
        for mass stopping power files. 
        
        Each file matching the pattern 'MassStopPwrData-MediumName.txt' is imported and is assigned to the medium "MediumName".\
        The imported data is added to the dictionary storing the stopping power data for each medium.
        
        '''
             
             
        dirPath = os.getcwd()+'/RefData/MassStoppingPower_NewData/'
        print "#### Importing new mass stoppping power data ####"
        for r,d,f in os.walk(dirPath):
            for files in f:
                if fnmatch.fnmatch(files, patternFilenames):
                    medium = os.path.splitext(files)[0].split('MassStopPwrData-')[1]
                    filename = os.path.join(r,files)
                    self.loadMassStoppingPowerFromNewData(filename,medium)
                    print "New Mass stopping power for %s imported"%medium
            
            print "#### End Importing new mass stoppping power data ####"
            return # stop after the first level of files has been processed


    def importCorrecStopPwrData(self):
        '''Function used to import stopping power correction factor data. It scans the folder */RefData/CorrecFactor_StopPwr/*\. 
        
        Each file matching the pattern 'dataCorrecFactorStopPwr_NAMEDATA.txt'.\
        NAMEDATA must also be provided by the mspt global variable : nameStopPwrCorrecFactor = 'NAMEDATA'
        
        '''
             
             
        dirPath = os.getcwd()+'/RefData/CorrecFactor_StopPwr/'
        print "#### Importing correction factor stopping power ####"
        self._CorrecFactorTable = None
        for r,d,f in os.walk(dirPath):
            for files in f:
                if fnmatch.fnmatch(files, patternFilenameCorrecFactor):
                    nameData = os.path.splitext(files)[0].split('dataCorrecFactorStopPwr_')[1]
                    if nameData == self._globalVar.nameStopPwrCorrecFactor:
                        filename = os.path.join(r,files)
                        loadedArray = np.loadtxt(filename,dtype="|S",delimiter='\t',skiprows = 1,usecols=(0,1))
                        loadedArray = loadedArray.astype(self._globalVar.typeFloat)
                        loadedArray = loadedArray.transpose()
                        loadedArray.view('%s,%s'%(self._globalVar.typeFloat,self._globalVar.typeFloat)).sort(order=['f1'], axis=0)
                        print "Stopping power correction factor %s imported"%nameData
                        self._CorrecFactorTable = np.array(loadedArray,dtype = self._globalVar.typeFloat, order = 'C')
            if self._CorrecFactorTable is None:
                print "Correction factor data (stopping power) for name : %s , has not been found"%self._globalVar.nameStopPwrCorrecFactor
                raise ValueError('Error - No stopping power correction factor data has been found in /RefData/CorrecFactor_StopPwr/. Please check the data and try again..')
            print "#### End Importing correction factor mass stopping power ####"
            return # stop after the first level of files has been processed


    def getStopPwrCorrecFactor(self):
        '''Get the stopping power correction factor table
        
        '''
        return self._CorrecFactorTable
            
    def precomputeMassStopPowerRatios(self):
        '''Compute relative mass stopping power (relative to water). For each medium divide the mass stopping power \
        by the mass stopping power of water.
        
        '''
    
        self._mstpRatios = dict()
        nbRowsData = self._massStopPwrAndRangeData['Water'].shape[0]
        
        if self._globalVar.importNewMassStopPwr:
            funcInterpWater = interp1d(self._massStopPwrAndRangeData['Water'][:,0], self._massStopPwrAndRangeData['Water'][:,1], kind = 'linear')
        
        for key in self._massStopPwrAndRangeData:
            if not key in self._importedDensities:
                newArray = np.zeros((nbRowsData,2))
                newArray[:,0] = self._massStopPwrAndRangeData[key][:,0] # Energies
                newArray[:,1] = self._massStopPwrAndRangeData[key][:,1] / self._massStopPwrAndRangeData['Water'][:,1]
                self._mstpRatios[key] =  newArray 
            else:
                newArray = np.zeros((self._massStopPwrAndRangeData[key].shape[0],2))
                newArray[:,0] = self._massStopPwrAndRangeData[key][:,0] # Energies
                for idx,ener in enumerate(newArray[:,0]):
                    newArray[idx,1] = self._massStopPwrAndRangeData[key][idx,1]/ funcInterpWater(ener)
                self._mstpRatios[key] =  newArray  
                    
        
    def getConversionTableStpPwrForEnergy(self,energy):
        '''Return a 3 x n array (n being the number of media):
            
            * row 0 : low bound of density range ,
            * row 1: high bound of density range,
            * row 2: rel mass stop pwr associated. Medium: Air,Adipose,Muscle,Bone
        
        :param energy: Energy in MeV
        
        :returns: COnversion table Density to Rel.ative Mass Stopping power: 
        
        .. note::
            
            Mass stopping power [MeV.cm^-1 / (g.cm^-3)] = [MeV.cm^2.g^-1]
        
        '''
        nbRows = len(self._mstpRatios['Water'][:,0])
        if energy < self._minEner or energy > self._maxEner:
            print "Energy = "+str(energy)
            strErr = "Energy out of range (%i-%i) in getConversionTableStpPwrForEnergy"%(int(self._minEner),int(self._maxEner))
            raise ValueError()
        
        if energy == self._minEner:
            idxEnergy = 0
        elif energy == self._maxEner:
            idxEnergy = nbRows - 1
        else:
            flagEnd = 0
            for idx,row in enumerate(self._mstpRatios['Water']):
                if flagEnd == 0:
                    if idx == nbRows - 1:
                        print "idx :%i , nbRows: %i, energy: %f"%(idx,nbRows,energy)
                        raise ValueError("Error indice in getConversionTableStpPwrForEnergy")
                    elif row[0] <= energy < self._mstpRatios['Water'][idx+1,0]:
                        idxEnergy = idx
                        flagEnd = 1

        if self._globalVar.importNewMassStopPwr:
            listMedia = list(globallistMedia[:])
            listDensities = list(globalmassDensities[:])
            
            for key in self._importedDensities:
                idxInsert = bisect.bisect_left(listDensities,self._importedDensities[key])
                listDensities.insert(idxInsert,self._importedDensities[key])
                listMedia.insert(idxInsert-1,key)
        else: 
            listMedia = globallistMedia
            listDensities = globalmassDensities
        
        
        
        table   = None
        idxMaxEnergy = nbRows - 1
        for key,low,high in zip(listMedia , listDensities[:-1],listDensities[1:]):
            
            #Process imported data
            if key in self._importedDensities:
                idxMaxEnergy_newData = len(self._mstpRatios[key][:,0])
                flagEnd = 0
                for idx,row in enumerate(self._mstpRatios[key]):
                    if flagEnd == 0:
                        if idx == nbRows - 1:
                            print "idx :%i , nbRows: %i, energy: %f - medium: %s"%(idx,nbRows,energy,key)
                            raise ValueError("Error indice in getConversionTableStpPwrForEnergy")
                        elif row[0] <= energy < self._mstpRatios[key][idx+1,0]:
                            idxEnergy_newData = idx
                            flagEnd = 1
                curTable = self.buildMediumToTable(key, low, high ,energy, idxEnergy_newData, idxMaxEnergy_newData)
            else: # Process regular data
                curTable = self.buildMediumToTable(key, low, high ,energy, idxEnergy, idxMaxEnergy)
            if table is None:
                table = curTable
            else:
                table =np.hstack((table,curTable))
        #WARNING DATA: ratio of mass stopping powers and not ratio of stopping powers!!
        return np.ascontiguousarray(table.astype(self._globalVar.typeFloat))
      
        
    def buildMediumToTable(self , key, low, high ,energy, idxEnergy, idxMaxEnergy):
        '''Function that builds the table to add to the conversion table for a given medium
        
        :param key: Name of the medium
        :param low: low bound of the density range
        :param high: high bound of the density range
        :param energy: energy in MeV
        :param idxEnergy: index of the energy in the mass stopping power table (PSTAR data).
        :param idxMaxEnergy: index of the maximum energy used in the mass stopping power table (PSTAR data).
        
        :returns: A table (3,1) :
            
            * row 0 : low bound of density range , \
            * row 1: high bound of density range, \
            * row 2: rel mass stop pwr associated to the medium for the given energy. 
        
        '''
                
        table = np.zeros((3,1))
        table[0,0] = low
        table[1,0] = high
        
        idxP1 = idxEnergy+1
        if idxEnergy != idxMaxEnergy:
            
            x1 = self._mstpRatios[key][idxEnergy,0]
            x2 = self._mstpRatios[key][idxP1,0]
            y1 = self._mstpRatios[key][idxEnergy,1]
            y2 = self._mstpRatios[key][idxP1,1]

            table[2,0] = tools.interp(energy, x1,x2,y1,y2)
        else:
            table[2,0] = self._mstpRatios[key][idxEnergy,1]

        return table