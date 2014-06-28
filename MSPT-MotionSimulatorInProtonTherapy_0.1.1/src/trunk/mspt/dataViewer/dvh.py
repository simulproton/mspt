########################################################################
#
# dvh.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# July 2013
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


import numpy as np
import os,sys
import plotDataGraphs as pltGraph
import time
import csv

pathToSave =  None  

class DVH(object):
    '''Class DVH is used to build and draw/save a Dose Volume Histogram (DVH)
    
    :param configData: configuration dictionary (*globalVariables*) for the simulator loaded in *simulator.simulatorMSPT*. It should contain\
    the keys: 
    
        * 'mainPath' : The path where the simulation directory will be saved.
        * 'nameSimulation' : The simulation name.
        * 'dvhRelativeVolume' : True or False: True to plot the DVH using a relative volume (0-100%), False to \
        plot with real volumes.
        
    :param dose: 3D dose distribution (in Gy) stored in a numpy array.
    :param listVOIs: List of the volumes of interest (VOIs). List of tuples:  ( mask 3D Numpy array , ROI name, ROI observation, ROI index, Type).\
    See the function *getMasksForAssociatedDicomImage()* in the module *dicomReader.dicomRSReader.* for more details.
    
    .. note::
        
            Output dose in cGy.
    
    
    '''
    def __init__(self,configData,dose,listVOIs):
        self._dictDVH = None #Dictionnary: 1 key / VOI:  self._dictDVH[key] = [xData , yData]
        self._figure = None # figure used to display the DVHs
        self._subplot = None
        self._dictDv = None # Dictionnary containing the Dv data for each VOI: self._dictDv[voi] = list[ [dose1,vol1],[dose2,vol2]  .....]
        self._dictVd = None # Dictionnary containing the Dv data for each VOI: self._dictVd[voi] = list[ [dose1,vol1],[dose2,vol2]  .....]
        self._globalVar = configData
        global pathToSave
        pathToSave = self._globalVar.mainPath + self._globalVar.nameSimulation + "DVH_Data/"
        self.buildDVHForVOIs(dose , listVOIs)
         
        
        
        
    def getDVHForVOI(self, voi):
        '''Get the DVH for a given VOI name. 
        
        :param voi: Name of the voi. If voi='all' returns a dictonary with all the DVHs.
        
        :returns: The DVH as a dictionary: dictDVH[voi]= [dose bins array , volume data array].

        .. note::
        
            Dose in cGy.
        
        '''
        if voi == 'all':
            return self._dictDVH
        elif voi in self._dictDVH.keys():
            return self._dictDVH[voi]
        else:
            print "Unknown voi %s"%str(voi)
            raise ValueError("VOI unknown in getDVHForVOI")
        
    def getDvDataForVOI(self, voi):
        '''Get the Dv (Dose received by a given volume) data  for a given VOI name. For example D50, would be the dose received \
        by 50% of the VOI volume.
        
        :param voi: Name of the voi. If voi='all' return a dictonary with all the DVHs.

        :returns: The Dv as a dictionary: dictDv[voi]= [ [dose1,vol1], [dose2,vol2], ... [dosen,voln] ].
      
        .. note::
        
            Dose in cGy.
            
        '''
        if voi == 'all':
            return self._dictDv
        elif voi in self._dictDv.keys():
            return self._dictDv[voi]
        else:
            print "Unknown voi %s"%str(voi)
            raise ValueError("VOI unknown in getDvDataForVOI")  
        
        
    def getVdDataForVOI(self, voi):
        ''' Get the Vd (portion of the volume receiving a given dose) Data for a given VOI name. For example V10: volume of VOI receiving \
        10 dose units.

        :param voi: Name of the voi. If voi='all' return a dictonary with all the DVHs.

        :returns: The Vd as a dictionary: dictVd[voi]= [ [dose1,vol1], [dose2,vol2], ... [dosen,voln] ].

        .. note::
        
            Dose in cGy.

        '''
        if voi == 'all':
            return self._dictVd
        elif voi in self._dictVd.keys():
            return self._dictVd[voi]
        else:
            print "Unknown voi %s"%str(voi)
            raise ValueError("VOI unknown in getVdDataForVOI")  
    
   
    def buildDVHForVOIs(self, doseDistr, listVOIs, spacing = None):
        '''Function that receives a list of VOI and build DVH for each of these VOIs.
        
        :param doseDistr: 3D dose distribution (3D array) in Gy.
        :param listVOIs: List of VOIs : each VOI in the list should be stored in the shape provided by *dicomRSReader*: \
        [  ( mask 3D Numpy array , ROI name, ROI observation, ROI index) ,   ( ... )...   ]
        :param spacing: Voxels spacing: used if relativeVolume is False. If spacing is None and relativeVolume is False, \
        then it will be considered that voxel spacing is {'frames':10,'rows':10'cols':10}. Spacing is in mm.
        
        .. note::
            
            A **Dose Volume Histogram** (DVH) is a cumulative histogram for a specified volume of interest (VOI). It \
            represents the cumulative volume irradiated as a function of the dose. See ICRU report: `5 Geometric Terms, and Dose \
            and Dose-Volume Definitions <http://jicru.oxfordjournals.org/content/7/2/83.extract>`_  for more details.\
            To compute it:
        
                #. build the histogram nb voxels vs received dose
                #. build the cumulative histogram from 1): this represents the number of voxels receiving a dose <= tp the given dose
                #. build the DVH from 2): DVH = nb_tot_voxels_in_VOI - cumulative_histogram: this represent the number of voxels receiving  at least the given dose
    
        
        '''
        #Convert dose distribution in Gy into cGy
        doseDistr = 100.0 * doseDistr
        
        self._dictDVH = dict()
        self._dictDVH['relativeVolume'] = self._globalVar.dvhRelativeVolume
        for voi in listVOIs:
            if voi[4] == 'CLOSED_PLANAR' :# i.e. if not 'POINT'     
                doseVOI , indMask = getSubDoseDataForMask( doseDistr , voi[0])
                cumSum , bins = buildCumulHisto(doseVOI, indMask = indMask)
                xData =  bins
                yData = cumSum
                self._dictDVH[voi[1]] = [xData , yData]   
                
        self.normalizeDVH(listVOIs, spacing)
        print "DVHs built"
    
    def normalizeDVH(self,  listVOIs, spacing):
        '''Normalizes the DVH such that the volumes of each VOI is set in the range 0% - 100%.

        :param listVOIs: List of VOIs : each VOI in the list should be stored in the shape provided by\
        DicomRSReader: [  ( mask 3D Numpy array , ROI name, ROI observation, ROI index, Type) ,   ( ... )...   ]
        :param spacing: Voxels spacing: used if relativeVolume is False. If spacing is None and relativeVolume is False, \
        then it will be considered that voxel spacing is {'frames':10,'rows':10'cols':10}. spacing in mm

        '''
        for voi in listVOIs:
            if voi[4] == 'CLOSED_PLANAR' :# i.e. if not 'POINT'       

                xData = self._dictDVH[voi[1]][0]
                yData = self._dictDVH[voi[1]][1]
                maxY = np.max(yData)
                if self._dictDVH['relativeVolume']:
                    if maxY != 0:
                        yData = yData * 100 / maxY
                else:
                    if spacing is None:
                        spacing = {'frames':10,'rows':10,'cols':10}
                    yData = yData * (spacing['frames'] / 10.0 )* (spacing['rows'] / 10.0) * (spacing['cols']/ 10.0)
                self._dictDVH[voi[1]][1] = yData 

    
    
    def computeDvData(self,listVolumes, voiName = None):
        ''' Computes Dv (Dv indicates that a volume V of VOI receives at least a dose equal to Dv.\
        See ICRU report: `5 Geometric Terms, and Dose and Dose-Volume Definitions <http://jicru.oxfordjournals.org/content/7/2/83.extract>`_ \
        for more details.

        :param listVolumes: List of volume values to consider for Dv data.
        :param voiName: Name of the VOI. If voiName is None, it will be computed for all the VOIs.
        
        '''
        
        if isinstance(listVolumes,(float,np.float32,np.float64)) or isinstance(listVolumes,int):
            listVolumes = [listVolumes]
        if isinstance(listVolumes,list):
            self._dictDv = dict()
            allVOIs = True
            if voiName is not None and voiName not in self._dictDVH.keys():
                print "VOI %s not found in : %s.\n Dv data will be computed for all the VOIs."%(str(voiName),str(self._dictDVH.keys()))
            elif voiName in self._dictDVH.keys():
                allVOIs = False
                
            for voi in self._dictDVH:
                if voi != 'relativeVolume':
                    self._dictDv[voi] = list()
                    dvhData = self._dictDVH[voi]
                    for vol in listVolumes:
                        if vol >= 0:
                            dose = getDoseForVolume(  dvhData[0], dvhData[1] ,vol )
                            self._dictDv[voi].append([dose, vol])
                    if not allVOIs and voiName == voi:
                        break
        else:
            print "listVolumes in computeDvData is not a list or a float or an int"
            raise ValueError('Wrong type for listVolumes in computeDvData')
       
       
    def computeVdData(self,listDoses, voiName = None):
        '''Computes Vd (Vd indicates that the volume Vd of VOI that received at least D) values. \
        See ICRU report: `5 Geometric Terms, and Dose and Dose-Volume Definitions <http://jicru.oxfordjournals.org/content/7/2/83.extract>`_  \
        for more details.
        
        :param listDoses: List of dose values to consider for Vd data.
        :param voiName: Name of the VOI. If voiName is None, it will be computed for all the VOIs.
        
        
        '''
        
        if isinstance(listDoses,(int,float,np.float64,np.float32)) :
            listDoses = [listDoses]
        if isinstance(listDoses,list):
            self._dictVd = dict()
            allVOIs = True
            if voiName is not None and voiName not in self._dictDVH.keys():
                print "VOI %s not found in : %s.\n Vd data will be computed for all the VOIs."%(str(voiName),str(self._dictDVH.keys()))
            elif voiName in self._dictDVH.keys():
                allVOIs = False
                
            for voi in self._dictDVH:
                if voi != 'relativeVolume':
                    self._dictVd[voi] = list()
                    dvhData = self._dictDVH[voi]
                    for dose in listDoses:
                        if dose >= 0:
                            vol = getVolumeForDose(  dvhData[0], dvhData[1] ,dose )
                            self._dictVd[voi].append([dose, vol])
                    if not allVOIs and voiName == voi:
                        break
        else:
            print "listDoses in computeVdData is not a list or a float or an int"
            raise ValueError('Wrong type for listDoses in computeVdData')
    
    
    def drawDVHs(self, voiName = None):
        '''Draws the DVH in a pyplot figure.
        
        
        :praram voiName:  None by default in order to display DVH for all VOIs. If not none it should contain the name of the VOI to display.\
        If the name is unknown it will display the DVH for all the VOIs.
        
        
        '''
        if self._dictDVH is None:
            print "DVH: nothing to draw..."
            return
        else:
            self._figure = pltGraph.getCurrentFigure()
            
            allVOIs = True
            if voiName is not None and voiName not in self._dictDVH.keys():
                print "VOI %s not found in : %s.\n All the DVHs will be plotted."%(str(voiName),str(self._dictDVH.keys()))
            elif voiName in self._dictDVH.keys():
                allVOIs = False
    
            supTitle = 'DVH for volumes'
            count = 0
            self._subplot = pltGraph.newSubplot(self._figure )
            xMax = None
            for voi in self._dictDVH:
                if voi != 'relativeVolume':
                    dvhData = self._dictDVH[voi]
                    maxi = np.max(dvhData[0])
                    if xMax is None or maxi > xMax:
                        xMax = maxi
            
            for voi in self._dictDVH:
                if voi != 'relativeVolume':
                    if not allVOIs and voiName == voi:
                        pltGraph.clearSubplot(self._subplot)
                        supTitle = 'DVH '+ str(voi)
                    dvhData = self._dictDVH[voi]
                    if self._dictDVH['relativeVolume']:
                        yMax = 100
                        yLabel = 'Relative Volume %'
                    else:
                        yMax = np.max( dvhData[1]) + 10
                        yLabel = 'Volume (cm3)'
                    pltGraph.drawPlot(self._subplot , dvhData[0], dvhData[1],xmin = 0, xmax = xMax, ymin = 0, ymax =yMax,xlabel = 'Dose (cGy)', ylabel=yLabel, titlePlot = 'DVH '+ str(voi), singlePlot = True, nameSinglePlot = supTitle, grid = True)
                    supTitle = None
                    pltGraph.displayPlot()
                    if not allVOIs and voiName == voi:
                        break
                    
            pltGraph.makePlotLegend(self._figure,self._subplot)
            pltGraph.displayPlot()
    
    def drawDvData(self , voiName = None):
        '''Displays Dv Data in the pyplot figure.
        

        :praram voiName:  None by default in order to display DVH for all VOIs. If not none it should contain the name of the VOI to display.\
        If the name is unknown it will display the DVH for all the VOIs.

        
        '''
        if self._dictDv is None:
            print "Nothing to draw for Dv data"
        else:
            if self._figure is None or self._subplot is None:
                print "No dvh has been drawn yet. No Dv data will be displayed."
            else:
                pltGraph.displayPlot()
            
                allVOIs = True
                if  voiName is not None and voiName not in self._dictDv.keys():
                    print "VOI %s not found in : %s.\n All Dv data will be plotted."%(str(voiName),str(self._dictDv.keys()))
                elif voiName in self._dictDv.keys():
                    allVOIs = False           
                for voi in self._dictDv:
                    for dvData in self._dictDv[voi]: 
                        pltGraph.plotHorizontalLine(self._subplot, dvData[1] , xEnd = dvData[0] ,color = 'b')
                        pltGraph.plotVerticalLine(self._subplot, dvData[0], yEnd =  dvData[1] , color = 'b', label = r'$D_{%i}$'%int(dvData[1]), font = globalFont)
                        pltGraph.displayPlot()
                    if not allVOIs and voiName == voi:
                        break        
            
    def drawVdData(self, voiName = None):
        '''Displays Vd Data in the pyplot figure.
        
        :praram voiName  None by default in order to display DVH for all VOIs. If not none it should contain the name of the VOI to display.\
        If the name is unknown it will display the DVH for all the VOIs.

        '''
        if self._dictVd is None:
            print "Nothing to draw for Vd data"
        else:
            if self._figure is None or self._subplot is None:
                print "No dvh has been drawn yet. No Vd data will be displayed."
            else:
                pltGraph.displayPlot()
            
                allVOIs = True
                if voiName is not None and voiName not in self._dictVd.keys():
                    print "VOI %s not found in : %s.\n All Dv data will be plotted."%(str(voiName),str(self._dictVd.keys()))
                elif voiName in self._dictVd.keys():
                    allVOIs = False           
                for voi in self._dictVd:
                    for vdData in self._dictVd[voi]: 
                        pltGraph.plotHorizontalLine(self._subplot, vdData[1] , xEnd = vdData[0] ,color = 'r',label = r'$V_{%i}$'%int(vdData[0]), font = globalFont)
                        pltGraph.plotVerticalLine(self._subplot, vdData[0], yEnd =  vdData[1] , color = 'r' )
                        pltGraph.displayPlot()
                    if not allVOIs and voiName == voi:
                        break      
        
    def saveDVHAsImage(self, filename):
        '''Saves the current DVH (in a pyplot figure) as an image.
        
        :param filename: Name of the image file to save.
        
        '''
        if self._figure is not None:
            pltGraph.savePlotImage(self._figure, 'DVH_%s'%(str(filename)) , dirPath = pathToSave, ext = 'png', subplot = None)
        else:
            print "Plot as not been drawn yet"
        
    def saveDVHAsCSV(self, filename):
        '''Saves DVH data to CSV files: 1 per VOI.
        
        :param filename: Basis name of the CSV files to save.
        
        '''
        if self._dictDVH is not None:
            exportDVHToCSVFile(self._dictDVH, filename)
        else:
            print "DVHs has not been computed yet, so they cannot be exported to CSV files"
     
     
    def saveEntireDVHAsCSV(self, filename):
        '''Saves an entire DVH data to a CSV file.
        
        :param filename: Filename basis of the CSV file to save.
        
        '''
        if self._dictDVH is not None:
            exportAllDVHToCSV(self._dictDVH, filename)
        else:
            print "DVHs has not been computed yet, so they cannot be exported to CSV files"
     
     
     
    def saveDvAsCSV(self, filename):
        '''
        Saves Dv data to CSV files: 1 per VOI
        
        :param filename: Filename basis of the CSV files to save.
        
        '''
        if self._dictDv is not None:
            volumeType = 'Absolute_cm3'
            if self._dictDVH['relativeVolume']:
                volumeType = 'Relative_%'
            exportDvDataToCSVFile(self._dictDv, volumeType, filename)
        else:
            print "Dv data has not been computed yet, so it cannot be exported to CSV files"   
            
    def saveVdAsCSV(self, filename):
        '''Saves Dv data to CSV files: 1 per VOI.
        
        :param filename: Basis name of the CSV files to save.
        
        '''
        if self._dictVd is not None:
            volumeType = 'Absolute_cm^3'
            if self._dictDVH['relativeVolume']:
                volumeType = 'Relative_%'
            exportVdDataToCSVFile(self._dictVd, volumeType, filename)
        else:
            print "Vd data has not been computed yet, so it cannot be exported to CSV files"   
            
    def closeDVH(self):
        '''Closes pyplot figure.
        
        '''
        pltGraph.closePlot()
    

#######################################################
##                  Utilities: save data DVH
#######################################################

def exportDVHToCSVFile(dictDVH, name):
    '''Saves a DVH dictionary to csv files. One file per VOI.
    
    :param dictDVH: DVH dictionary to save.
    :param name: Name of the DVH. Used as the filename basis.
    
    '''
    volumeType = 'Absolute_cm^3'
    if dictDVH['relativeVolume']:
        volumeType = 'Relative_%'
    for voi in dictDVH:
        if voi != 'relativeVolume':
            tmpListDose = list()
            tmpListDose[:] = dictDVH[voi][0][:]
            tmpListDose.insert(0,'%s_Dose_cGy'%str(voi))
            tmpListVol = list()
            tmpListVol[:] = dictDVH[voi][1][:]
            tmpListVol.insert(0,'%s_Volume_%s'%(str(voi),volumeType))
            dataList = [tmpListDose , tmpListVol]
            save2DListToCSV(dataList, 'DVH_%s_Vol_%s_%s'%((str(voi).replace(" ","_")),volumeType,str(name)))



def exportAllDVHToCSV(dictDVH , name):
    '''Saves a DVH dictionary to a csv files. One file for all the VOIs.
    
    :param dictDVH: DVH dictionary to save.
    :param name: Name of the DVH. Used as the filename.
    
    '''
    if name.endswith(".csv"):
        nane = name[0:-4]
    dataToExport = list()
    title = ("DVH %s"%(name)).replace(" ","_")
    xLabel = "Dose-cGy_"
    yLabel = "Volume-%_"

    idx = 0
    for roi in dictDVH:
        if roi != 'relativeVolume':
            listDataDose = list()
            listDataDose[:] = dictDVH[roi][0][:]
            listDataVol = list()
            listDataVol[:] = dictDVH[roi][1][:]
            listDataDose.insert(0,str(xLabel+roi))
            listDataVol.insert(0,str(yLabel+roi))
            dataToExport.append(listDataDose)
            dataToExport.append(listDataVol)
        
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)             

    with open(pathToSave +title+'.csv',"wb") as f:
        f.write("\n".join(",".join(map(str, x)) for x in dataToExport))
    
    return   pathToSave +title+'.csv'
    



def exportDvDataToCSVFile(dictDv, volumeType, name):
    '''Saves Dv data to csv files. One file per VOI.
    
    :param dictDV: Dictionary containing the (Dv, Volume) tuples for each VOI.
    :param volumeType: 'Relative_%' or 'Absolute_cm^3'
    :param name: Filename basis.
    
    '''
    for voi in dictDv:
        if voi != 'relativeVolume':
            tmpList = list()
            tmpList[:] = dictDv[voi][:]
            tmpList.insert(0,['%s_Dose_cGy'%(str(voi)) , '%s_Volume_%s'%(str(voi),volumeType) ])
            save2DListToCSV(tmpList,'DvData_%s_Vol_%s_%s'%(str(voi),volumeType,str(name)))
            
def exportVdDataToCSVFile(dictVd,volumeType, name):
    '''Saves Vd data to csv files. One file per VOI.
    
    :param dictDV: Dictionary containing the (Dv, Volume) tuples for each VOI.
    :param volumeType: 'Relative_%' or 'Absolute_cm^3'
    :param name: Filename basis.
    
    '''
    for voi in dictVd:
        if voi != 'relativeVolume':
            tmpList = list()
            tmpList[:] = dictVd[voi][:]
            tmpList.insert(0,['%s_Dose_cGy'%(str(voi)) , '%s_Volume_%s'%(str(voi),volumeType) ])
            save2DListToCSV(tmpList,'VdData_%s_Vol_%s_%s'%(str(voi),volumeType,str(name)))


def save2DListToCSV(data,filename):
    '''Saves a 2D list to a CSV file.
    
    :param data: 2D list of values.
    :param filename: CSV filename.
        
    '''
    #Save a 2D python list to csv file
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)             
    if not filename.endswith('.csv'):    
            filename = filename + ".csv"
    myFile = open(pathToSave+filename, 'wb')
    wr = csv.writer(myFile, quoting=csv.QUOTE_ALL)
    for row in data:
        wr.writerow(row)
    myFile.close()
        

    



#######################################################
##                  Utilities: data histogram
#######################################################
def getSubDoseDataForMask( doseData , mask):
    '''Get the dose distribution in the volume defined by the mask.
    
    :param doseData: 3D dose distribution (3D numpy array).
    :param mask: 3D numpy binary array (1 in the VOI, 0 itherwise).
    
    :returns: A tuple with:
    
        #. A 3D volume where doseData is preserved in the VOI (using mask) and is set to 0 outside.
        #. Indices where the mask is 1 (VOI).
 
    '''
    newDoseData = doseData * mask
    return (newDoseData , np.where(mask == 1))
    
    
def buildHistogram(data3D , indMask = None, maxValue = None):
    '''Builds an histogram for data3D. 
    
    :param data3D: 3D array for which we want to build a histogram.
    :param indMask: Indices where the VOI is located. If indMask is not None, the histogram is built only for the voxels \
    specified by indMask, otherwise it builds an histogram for the entire array.
    :param maxValue: If maxValue is None a maxValue is defined as the maximum value of the 3D volume used ( limited to indMask \
    if indMask is not None), otherwise the value provided is used. 
    
    
    The histogram is built with bins ranging from 0 to maxValue with a step of 0.5

    :returns: A tuple with:
    
        #. y values of the histogram (array)
        #. x values (bins) of the histogram (array)

    
    '''
    if indMask == None:
        data = data3D
    else:
        data = data3D[indMask]
        
    # Extend range for values in order to make sure that count will reach 0    
    if maxValue == None:
        maxVal = np.max(data) + 10 
    else:
        maxVal = maxValue + 1
    hist, bins = np.histogram(data,bins = np.arange(start = 0,stop = maxVal, step = 0.5  ))
    return (hist,bins)
    
    
    
def buildCumulHisto( data3D , indMask = None, maxValue = None ):
    '''Builds a cumulative histogram for data3D.
    
    :param data3D: 3D array for which we want to build a histogram.
    :param indMask: Indices where the VOI is located. If indMask is not None, the histogram is built only for the voxels \
    specified by indMask, otherwise it builds an histogram for the entire array.
    :param maxValue: If maxValue is None a maxValue is defined as the maximum value of the 3D volume used ( limited to indMask \
    if indMask is not None), otherwise the value provided is used.

    The histogram is built with bin ranging from 0 to maxValue with a step of 0.5.\
    To build the cumulative histogram the function *buildHistogram()* is called.
    
    :returns: A tuple with:
    
        #. y values of the histogram (array)
        #. x values (bins) of the histogram (array)
        
    '''

    if indMask == None:
        data = data3D
    else:
        data = data3D[indMask]
    
    if maxValue == None:
        max = np.max(data) + 10
    else:
        max = maxValue        
    hist, bins = buildHistogram(data3D , indMask = indMask ,maxValue = max)
    
    doseBins = []
    cumSumHist = []
    
    count = 0
    for dose in bins[:-1]:
        doseBins.append(dose)
        cumSum = np.sum(hist[count:]) 
        cumSumHist.append(cumSum)
        count = count + 1
    
    return (np.array(cumSumHist),np.array(doseBins))
    
    
def getDoseForVolume( doseData, volData ,vol ):
    '''Calculate the Dv (Dose received by a given volume) Data  for a given VOI name. For example D50, would be the dose received by 50% \
    of the VOI volume. 
    
    See ICRU report: `5 Geometric Terms, and Dose and Dose-Volume Definitions <http://jicru.oxfordjournals.org/content/7/2/83.extract>`_  \
    for more details. 
    
    This function is called by computeDvData(). 
    
    :param doseData: 1D array with dose values of the DVH. It should be in inscreasing order.
    :param volData: 1D array with volume values of the DVH. It should be in decreasing order.
    :param vol: Value of the volume considered.
    
    :returns: Dv value calculated.
   
    '''
    
    #volData should be decreasing and dose data inscreasing
    if len(doseData) != len(volData):
        raise ValueError( "volData and doseData don't have the same size")
    
    #Invert data so that volData is in increasing order
    volData = volData[::-1]
    doseData =  doseData[::-1]
    
    
    if volData[0] == volData[-1]:
        return doseData[0] 
    
    count = 0
    for (v , vp1) in zip(volData[0:-1],volData[1:]):
        if vol > 0:
            if v < vol <= vp1:
                if (doseData[count]-doseData[count+1]) != 0:
                    slope = (doseData[count]-doseData[count+1])/(v - vp1)
                    intersect = doseData[count] - slope * v
                    value = (slope * vol + intersect)
                    if value < 1e-5:
                        value = 0
                    return value
                else:
                    return doseData[count]
        else: # vol == 0
            if v == vol and vp1 > vol:
                return doseData[count]
        count = count+1
    return 0
    
    
def getVolumeForDose(doseData, volData ,dose):
    '''Computes Vd (Vd indicates that the volume Vd of VOI that received at least D) values.
     
    See ICRU report: `5 Geometric Terms, and Dose and Dose-Volume Definitions <http://jicru.oxfordjournals.org/content/7/2/83.extract>`_  \
    for more details. 
    
    This function is called by computeVdData(). 
    
    :param doseData: 1D array with dose values of the DVH. It should be in inscreasing order.
    :param volData: 1D array with volume values of the DVH. It should be in decreasing order.
    :param dose: Value of the dose considered.
    
    :returns: Vd value calculated.
    
    '''
    #doseData should be inscreasing and volData decreasing
    if len(doseData) != len(volData):
        raise ValueError( "volData and doseData don't have the same size")

    if doseData[0] == doseData[-1]:
        return volData[0] 
        
    if dose == doseData[-1]:
        return volData[-1]
    
    count = 0
    for (d , dp1) in zip(doseData[0:-1],doseData[1:]):
        if d <= dose < dp1:
            if (volData[count]-volData[count+1]) != 0:
                slope = (volData[count]-volData[count+1])/(d - dp1)
                intersect = volData[count] - slope * d
                value = (slope * dose + intersect)
                if value < 1e-5:
                    value = 0
                return value
            else:
                return volData[count]

        count = count+1
    return 0  

        
#######################################################
##                  Utilities: external call to build easily a DVH 
#######################################################        
        
def processDVHForDoseDistribution(dose, listVOIs , listDv, listVd, nameDVH, dataConfig):
        '''Builds the DVH for a dose distribution and save the data computed: DVH, Dv, Vd in an image and in CSV files for each VOIs. 
        
        :param listVOIs: List of the volumes of interest (VOIs). List of tuples:  ( mask 3D Numpy array , ROI name, ROI observation, ROI index, type).\
        See *dicomReader.dicomRSReader.getMasksForAssociatedDicomImage()* for more details. 
        :param dose: 3D dose distribution (3D numpy array).
        :param listDv: List of volume values for which to compute the Dv value.
        :param listVd: List of dose values for which to compute the Vd value.
        :param nameDVH: Name for the files to save.
        :param dataConfig: configuration dictionary (*globalVariables*) for the simulator loaded in *simulator.simulatorMSPT*. It should contain\
    the keys: 
    
        * 'mainPath' : The path where the simulation directory will be saved.
        * 'nameSimulation' : The simulation name.
        * 'dvhRelativeVolume' : True or False: True to plot the DVH using a relative volume (0-100%), False to \
        plot with real volumes.
    
        :returns: (Dv data, Vd data)
    
            
        '''
        print ">>>>>>>>>>>> Processing DVH %s <<<<<<<<<<<<<<"%str(nameDVH)
        newDVH = DVH(dataConfig,dose , listVOIs)
        newDVH.computeDvData(listDv)
        dvData = newDVH.getDvDataForVOI('all')        
        newDVH.saveEntireDVHAsCSV(nameDVH)
        newDVH.saveDVHAsCSV(nameDVH)
        newDVH.saveDvAsCSV(nameDVH)
        newDVH.saveVdAsCSV(nameDVH)
        newDVH.computeVdData(listVd)
        vdData = newDVH.getVdDataForVOI('all')
        newDVH.drawDVHs()
        #newDVH.drawDvData()
        #newDVH.drawVdData()
        newDVH.saveDVHAsImage(nameDVH)
        newDVH.closeDVH()
        print ">>>>>>>>>>>> End Processing DVH %s <<<<<<<<<<<<<<"%str(nameDVH)
        return (dvData , vdData)
        
        