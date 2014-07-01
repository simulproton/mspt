########################################################################
#
# tools.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On 07/11/2012
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

import numpy as np
import os, sys
import dicom
from ..dataViewer import array2DViewer as viewer2D 

class Tools(object):
    '''Class Tools aims at providing utilities to the simulator.
    
    :param configData: Dictionary storing the MSPT configuration. Mandatory key is 'typeFloat' : 'float32' or 'float64'.

    '''
    
    
    def __init__(self, configData):
        self._globalVar = configData
        return


        
    def evaluateDifference(self,vol1, vol2): 
        """!Function that computes the difference between two 3D numpy arrays and calculate statistics: absDiff = abs( vol1 - vol2).\
        It also saves the results into a text file (EvalDiffDoseDistrib.txt) in the simulation directory.
        
        :param vol1: First 3D numpy array
        :param vol2: second 3D numpy array
        
        :returns: Dictionary with keys:
        
            * 'vol1': Dictionary with keys:
                
                * "min" : minimum value of vol1
                * "max" : maximum value of vol1
                * "std" : standard deviation value of vol1 
                * "avg" : average value of vol1

            * 'vol2': Dictionary with keys:
                
                * "min" : minimum value of vol2
                * "max" : maximum value of vol2
                * "std" : standard deviation value of vol2
                * "avg" : average value of vol2

            * 'absDiff': Dictionary with keys:
                
                * "min" : minimum value of absDiff
                * "max" : maximum value of absDiff
                * "std" : standard deviation value of absDiff
                * "avg" : average value of absDiff
        
        
        """
        #Comparision planned dose and computed dose
        if len(np.where(vol1 > 0)[0]) == 0 or len(np.where(vol2 > 0)[0]) == 0:
            print "In differences one of the volumes is entirely null."
            return
        diff = np.abs(vol1 - vol2)
        
        minVol1 = np.min(vol1[np.where(vol1>0)])
        maxVol1 = np.max(vol1)
        avgVol1 = np.mean(vol1)
        stdVol1 = np.std(vol1)
        
        minVol2 = np.min(vol2[np.where(vol2>0)])
        maxVol2 = np.max(vol2)
        avgVol2 = np.mean(vol2)
        stdVol2 = np.std(vol2)
        
        minDiff = np.min(diff)
        maxDiff = np.max(diff)
        avgDiff = np.mean(diff)
        stdDiff = np.std(diff)
        
        print "Evaluate differences:"
        print "\t\tmin\t\tmax\t\tmean\t\tstd"
        print "Truth set\t%e\t%e\t%e\t%e"%(minVol1,maxVol1,avgVol1,stdVol1)
        print "Comput. set\t%e\t%e\t%e\t%e"%(np.min(vol2[np.where(vol2>0)]),np.max(vol2),np.mean(vol2),np.std(vol2))
        print "Abs. Diff.\t%e\t%e\t%e\t%e"%(minDiff,maxDiff,avgDiff,stdDiff)
        print "\n"
        

        pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation
        if not os.path.exists(pathToSave):
            os.makedirs(pathToSave)


        with open(pathToSave + 'EvalDiffDoseDistrib.txt', 'w') as myFile:
                myFile.write("Evaluate differences:\n")
                myFile.write("Truth set: Expected dose distribution\n")
                myFile.write("Comput. set: Computed dose distribution\n")
                myFile.write("\t\t\tmin\t\t\t\tmax\t\t\t\tmean\t\t\tstd\n")
                myFile.write("Truth set\t%e\t%e\t%e\t%e\n"%(minVol1,maxVol1,avgVol1,stdVol1))
                myFile.write("Comput. set\t%e\t%e\t%e\t%e\n"%(np.min(vol2[np.where(vol2>0)]),np.max(vol2),np.mean(vol2),np.std(vol2)))
                myFile.write("Abs. Diff.\t%e\t%e\t%e\t%e\n"%(minDiff,maxDiff,avgDiff,stdDiff))
                


        res = dict()
        res['vol1'] = { "min": minVol1, \
                        "max": maxVol1, \
                        "avg": avgVol1, \
                        "std": stdVol1}
                
        res['vol2'] = { "min": minVol2, \
                        "max": maxVol2, \
                        "avg": avgVol2, \
                        "std": stdVol2}
        res['absDiff'] = {  "min": minDiff, \
                            "max": maxDiff, \
                            "avg": avgDiff, \
                            "std": stdDiff}
        return res


    def getMatrixHotAndColdSpot(self,volRef, volTest):
        '''Function that computes the 3D array of hot and cold spots as: diff = volTest - volRef
        
        :param volRef: 3D array of reference
        :param volTest: 3D array tested
        
        :returns: 3D array diff.
        
        '''
        diff = volTest - volRef
        return diff

    def showHotColdSpots(self, volRef, volTest, title = None , pathToSaveImg = None):
        '''Function that display the hot and cold spots on a frame located at the "middle" of the dose distribution.
        
        :param volRef: 3D array of reference
        :param volTest: 3D array tested
        :param title: Title for the figure. Default: None
        :param pathToSaveImg: Path where to save the figure. Default: None
        
        
        '''
        
        diff = self.getMatrixHotAndColdSpot(volRef, volTest)
        slope = (1.0/(np.max(diff.astype(self._globalVar.typeFloat)) - np.min(diff.astype(self._globalVar.typeFloat))))
        intercept = -slope * np.min(diff.astype(self._globalVar.typeFloat))
        diffScaled =  diff.astype(self._globalVar.typeFloat) * slope + intercept
        indImg = np.where( diffScaled != intercept)
        frameImg = np.round((( np.min(indImg[0])+np.max(indImg[0]))/2.0))    
        viewer2D.displayImageFromArray(diffScaled[frameImg], title+' -frame %i'%frameImg)
        if pathToSaveImg is not None:
            path, filename = os.path.split(pathToSaveImg)
            viewer2D.saveFigure( path+'/',filename , ext = 'png', subplotIdx = 111)
        viewer2D.closePlot()
           
        
    def saveCompensatedRTPlan( self,rtPlanPath, newRTPlanPath, dictPlan):
        '''Function that save a dictionary storing a treatment plan to a RP dicom file.
        
        :param rtPlanPath: Path to the input RP dicom file.
        :param newRTPlanPath: Path to the new RP dicom file.
        :param dictPlan: Dictionary storing the new treatment plan : dictionary[beam angle][Energy Rounded][ "ScanSpotPositionMap" / "ScanSpotMetersetWeights"]
        
        '''
        rp = dicom.read_file(rtPlanPath)
        nbBeams = len(rp.IonBeamSequence)
        beamsToDel = []
        for idxBeam in range(nbBeams):
            angle = rp.IonBeamSequence[idxBeam].IonControlPointSequence[0].GantryAngle
            if angle in dictPlan.keys():
                newDataANgle = dictPlan[angle]
                isocenter=rp.IonBeamSequence[idxBeam].IonControlPointSequence[0].IsocenterPosition
                listLayersToRm = []
                for idxLayer,enerLayer in enumerate(rp.IonBeamSequence[idxBeam].IonControlPointSequence[::2]):
                    ener = np.floor(enerLayer.NominalBeamEnergy)
                    if ener in dictPlan[angle].keys():
                        enerLayer.ScanSpotPositionMap = dictPlan[angle][ener]["ScanSpotPositionMap"]
                        enerLayer.ScanSpotMetersetWeights = dictPlan[angle][ener]["ScanSpotMetersetWeights"]
                        nbWeights = len( dictPlan[angle][ener]["ScanSpotMetersetWeights"])
                        rp.IonBeamSequence[idxBeam].IonControlPointSequence[idxLayer+1].ScanSpotPositionMap = dictPlan[angle][ener]["ScanSpotPositionMap"]
                        rp.IonBeamSequence[idxBeam].IonControlPointSequence[idxLayer+1].ScanSpotMetersetWeights = list(np.zeros(nbWeights))
                    else:
                        print "Energy layer %s not in new plan"%str(ener)
                        listLayersToRm.append(idxLayer)
                        listLayersToRm.append(idxLayer+1)
                i = 0     
                listLayersToRm.sort()
                for idxLayer in listLayersToRm:
                    del rp.IonBeamSequence[idxBeam].IonControlPointSequence[idxLayer-i]
                    i += 1
                    
            else:
                print "Beam %s degree not in new plan"%str(angle)
                beamsToDel.append(idxBeam)
    
        i = 0
        for idx in beamsToDel:
            del rp.IonBeamSequence[idx - i]
            i = i+1

        rp.save_as(newRTPlanPath)
        
        
        