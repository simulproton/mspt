########################################################################
#
# scanningPath.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# June 2013
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
from scipy import ndimage
import dicom
import os,sys

'''

The module scanning path aims to manage scanning path from RP dicom files for proton-therapy. It is composed of classes:

    * ScanningPathMultipleBeams : Scanning path defining gantry angles: for each gantry angle a ScanningPathSingleBeam is created.
    * ScanningPathSingleBeam : Scanning path for 1 specific gantry angle
    

'''

#Orders for scanning path. More can be added
scanningPath = ['Random', 'Regular']
repaintingMethodSlice = ['IsolayerRepaint','ScaledRepaint', 'FullWeightRepainting']#Note: 'FullWeightRepainting': we deliver the full weight at once instead of o maximum weight or a scaled weight.
repaintingMethodVolume = ['Volumetric','NonVolumetric']

#Margin for compensation
typeAdded = 2
typeOrigin = 1
typeOutside = 0
typeIntRepaint = 'int16' #int8 :-128 to 127 
maxRepFactor = 32767 # max 127 because arrays of int16: -32768 to 32767



class ScanningPathMultipleBeams(object):
    '''Scanning path defining gantry angles: for each gantry angle a ScanningPathSingleBeam is created.
    
    :param rpData: Dicom data from pydicom for the RP dicom file.
    :param dictSettings: Dictionary of settings defining the behavior of the scanning path. It is explained in the class ScanningPathSingleBeam.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'. 
   
    **__iter__( )**: 
    
    Iterate through each beam and build the ScanningPathSingleBeam structure for each beams: *Ion Beam Sequence* tag in dicom documentation (i.e. each gantry angle). 
    
    :returns: A tuple (the current scanning path:ScanningPathSingleBeam, the angle, the isocenter, virtual source axis distances, strAngle). strAngle = 'NewAngle' if the gantry angle stored in the previous *Ion Beam Sequence* is different from the current gantry angle.
    
    **__len__( )**:
    
    :returns: Number of *Ion Beam Sequence*, i.e. number of beams.
    
    
    '''

    def __init__(self, rpData, dictSettings,typeFloat ):
        self._rpData = rpData
        self._dictSettings = dictSettings
        self._typeFloat = typeFloat
        
    def __iter__(self):
        '''
        Iterate through each beam and build the ScanninPath structure for each beams. 
        
        :returns: A tuple (the current scanning path:ScanningPathSingleBeam, the angle, the isocenter, virtual source axis distances, strAngle). strAngle = 'NewAngle' if the gantry angle stored in the previous *Ion Beam Sequence* is different from the current gantry angle.

        
        '''
        prevAngle = None
        for indAngle , beam  in enumerate(self._rpData.IonBeamSequence):
            angle = beam.IonControlPointSequence[0].GantryAngle
            isocenter = beam.IonControlPointSequence[0].IsocenterPosition
            vsad = beam.VirtualSourceAxisDistances
            sPath = ScanningPathSingleBeam(self._dictSettings,self._typeFloat)
            sPath.buildScanningPath(self._rpData,beamIdx = indAngle) 
            if prevAngle is None:
                prevAngle = angle
                strAngle = 'NewAngle'
            elif prevAngle != angle:
                prevAngle = angle
                strAngle = 'NewAngle'
            else:
                strAngle = 'SameAngle'
            
            yield (sPath, angle,isocenter,vsad,strAngle)
    def __len__(self):
        '''
        :returns: Number of *Ion Beam Sequence*, i.e. number of beams.
    
        '''
    
        return len(self._rpData.IonBeamSequence)        


class ScanningPathSingleBeam(object):
    '''Scanning path defining the energies, spots positions and spots weights for each gantry angle. The key idea of this class is to represent the\
    scanning path of a gantry angle as a set of 3D matrices where each frame corresponds to an energy layer. Then each pixel in a frame (or each voxel \
    if we consider the 3D matrix entirely) represent spatially a possible scan spot position. The information contained in the set of 3D matrices will tell if the scan spot positions\
    are part of the treatment plan (i.e., the proton beam should irradiate what is along the pencil beam).
    
    
    :param dictSettings: Dictionary of settings defining the behavior of the scanning path. It is explained below.
    :param typeFloat: The type of numpy float to use. It should  be either 'float32' or 'float64'. 
    
    The set of 3D matrices representing a scanning path is composed of:
        
        * 1 matrix storing the **X positions** of the scan spots. If a voxel (frame0,row0,column0) has a value of 0, it means that this voxel does not \
        belong to the treatment plan. If the value is non-null then it is part of the plan. Spatially, moving of 1 frame means changing energy, and moving of\
        1 row (resp. column) means moving the pencil beam of row_spacing mm (resp. column_spacing mm). The spacing is the smallest y (resp. x) difference between\
        two consecutive scan spots.
        * 1 matrix storing the **Y positions** of the scan spots. Similar to the previous matrix.
        * 1 matrix storing the **planned weight** for each spot position. The weight is set to 0 for scan spots that do not belong to the plan.
        * 1 matrix storing the **planned weight for each spot position at the time of the delivery**. Therefore the weight is set to 0 for scan spots that do not \
        belong to the plan or that have already been delivered. This matrix is updated during the delivery.
        * 1 matrix storing the **type of scan spot position**. If the type is '0', the pencil beam will never stop at this position. If the 'type' is '1', it means\
        that the pencil beam will stop at this position and that this spot is part of the treatment plan. If the 'type' is '2', it means that under certain circumstances\
        (i.e. if 'AddMapMargins' if dictSettings is True - see below) the pencil beam will go by this position, but the weight affected to it is null. This \
        allow to have more spot positions. This is especially used if the delivery is compensated.
        * 1 matrix storing the **number of repainting** of each scan spot: '0' if the pencil beam will not go by a specific scan position, a positive integer otherwise.\
        This matrix is updated during the delivery.
    
   
    **dictSettings**:
    
    dictSettings is a dictionary with mandatory keys and optional keys that depend on the desired behavior of the scanning path. Keys and values are defined here:
    
        * **'TypeScanning'** (mandatory):  values:  'Random' or  'Regular'.  'Random': the spot positions are chosen randomly both for the energy and the position. 'Regular': the spot positions are chosen from the highest energy to the smallest and in the order provided by the RP dicom file.
        * **'Repainting'** (mandatory):  values:  True, if one wants to deliver the treatment plan using the repainting technique, False otherwise.
        * **'SliceRepaint'** (mandatory if 'Repainting' is True): method to use to repaint an energy layer.  values: ['IsolayerRepaint','ScaledRepaint', 'FullWeightRepainting'] :
                * *'IsolayerRepaint'*: a maximum weight will be calculated as the minimum weight of all the weight defined in the treatment plan.   
                * *'ScaledRepaint'* : a fixed number of repaintings is defined.
                * *'FullWeightRepainting'* : used for compensation: we deliver the full weight instead of a maximum weight or a scaled weight. If not used with the compensation it is equivalent to no repainting. It can also bee seen as a scaled repainting of repainting factor 1.
        * **'VolumeRepainting'** (mandatory if 'Repainting' is True): method used to repaint the entire volume. Values: ['Volumetric','NonVolumetric']:
                * *'Volumetric'*: Scan energy layer by energy layer to paint the volume a first time then repaint the volume layer by layer and so on.
                * *'NonVolumetric'*: Repaint each layer the number of time desired before moving to the next energy layer.
        * **'UnlimitedRescan'** (mandatory if 'ScaledRepaint' or  'FullWeightRepainting' and if 'RepaintFactor' is not present): values: \
        True if one wants to allow an unlimited number of repainting, i.e. not controlled by an maximum number of repaintings.\
        However, to avoid infinite loops, a maximum number of repainting is however set to 32767 (max value of an int16). \
        One must note that if 'RepaintFactor' is present, the value of 'UnlimitedRescan' will have a higher priority.
        * **'RepaintFactor'** (mandatory if 'ScaledRepaint' or  'FullWeightRepainting' and if 'UnlimitedRescan' is not present): \
        value: integer value (1-127) defining the number of repaintings desired.\
        One must note that if 'UnlimitedRescan' is present, the value of 'UnlimitedRescan' will have a higher priority.
        * **'AddMapMargins'** (optional): One has to note that we call "map" the set of spot positions within an energy layer (name used in Dicom documentation Scan Spot Position Map).\
        Values: True if one wants to add "margins" (i.e. spot positions of weight 0) around the native maps defined in the energy layers of the RP dicom file. This option is more dedicated to the compensation technique. \
        False (default value): if one do not want to add margins to the map. The margin is created using a maximum thickness of expansion (defined by 'Margin Params') and by performing a dilation. Each spot of the scanning path is viewed as\
        a voxel in a 3D binary matrix, the dilation is applied to the volume represented by these non-null voxels.
        * **'MarginParams'** (mandatory if 'AddMapMargins' is True): value: an iterable of positive real values. It represents the width of the margin in mm along the [x,y,z] axis in the IEC Gantry coordinate system (i.e. coordinate system in which spot positions are defined).
        * **'3DMargins'** (mandatory if 'AddMapMargins' is True): values: True if the the map is also extended along the z axis of the IEC gantry coordinate system. In other words, this adds energy layers in addition to spot positions. \
        False, if one do not want to extend the map in 3D. In this case the 3rd element of 'MarginParams' is discarded.
        * **'UpdateMargins'** (mandatory if 'AddMapMargins' is True): controls whether the map with the margins are updated between repaintings. From one repainting to another the map of an energy layer can change depending on the spots that have been entirely deliverd. \
        Values: True or False.
        * **'FindStartPos'** (optional): True or False. It is used for the compensation. If set to True, it will find the first spot, of an energy layer, that could be delivered in a map for a given motion vector.
        
     
     **__len__( )**: 
     
     :returns: The number of energy layers.
     
     **__iter__( )**:
     
     Iterates through the spots following the order defined in dictSettings.       
     
     :returns: A tuple:
     
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

    def __init__(self, dictSettings,typeFloat):
        
        self._typeFloat = typeFloat
        
        #Process settings in dictSettings:
        if "TypeScanning" not in dictSettings.keys():
            raise ValueError("No TypeScanning key in scanning path settings")
        typeSet =  dictSettings["TypeScanning"]
        if typeSet not in scanningPath:
            strErr = "Unknown type in scanning path settings: %s "%str(typeSet)
            raise ValueError(strErr)
    
        #Initialize scanning path variable for non random delivery
        self._scanPath = typeSet
        if self._scanPath != scanningPath[0]:# : random , regular
            self._currEnergySlice = None #store the index of the current energy slice
            self._sliceScanningOrder = None #store the list of the spots of the current energy slice. The order in the list defines to order that will be
                                            #used when iterating through the list of spots.
            self._nextIndex = None # store the index of the next spot in the current energy slice
            
        
        #Verify that the key "Repainting" is present and that is a boolean.
        if "Repainting" not in dictSettings.keys():
            raise ValueError("No Repainting key in scanning path settings")
        if not isinstance(dictSettings["Repainting"],bool):
            strErr = "Value for Repainting key in scanning path should be a bool not a %s"%str(type(dictSettings["Repainting"]))
            raise ValueError(strErr)
        
        
        # If the user wants to use the repainting
        if dictSettings["Repainting"]:
            self._repaintingIndex = None
            #Slice repainting
            if "SliceRepaint" not in dictSettings.keys():
                raise ValueError("No SliceRepaint key in scanning path settings")
            if dictSettings["SliceRepaint"] not in repaintingMethodSlice:
                strErr = "Value for SliceRepaint is not in %s : value: %s"%(str(repaintingMethodSlice),str(dictSettings["SliceRepaint"]))
                raise ValueError(strErr)
            
            #Method slice repainting
            self._repaintingMethodSlice = dictSettings["SliceRepaint"]
                
                #'IsolayerRepaint'
            if self._repaintingMethodSlice == repaintingMethodSlice[0]: #['IsolayerRepaint','ScaledRepaint', 'FullWeightRepainting']
                self._maxWeight = None
                #Method Volume repainting
                if "VolumeRepainting" not in dictSettings.keys():
                    raise ValueError("VolumeRepainting is not in the scanning path settings")
                if dictSettings["VolumeRepainting"] not in repaintingMethodVolume:
                    strErr = "Value for VolumeRepainting is not in %s : value: %s"%(str(repaintingMethodVolume),str(dictSettings["VolumeRepainting"]))
                    raise ValueError(strErr)
                self._repaintingMethodVolume = dictSettings["VolumeRepainting"]
                    
                    #'ScaledRepaint' or  'FullWeightRepainting'
            elif self._repaintingMethodSlice == repaintingMethodSlice[1] or \
                    self._repaintingMethodSlice == repaintingMethodSlice[2]:
                
                if "UnlimitedRescan" not in dictSettings.keys():
                    if "RepaintFactor" not in dictSettings.keys():
                        raise ValueError("UnlimitedReapint and RepaintFactor are not in the scanning path settings.")
                    else:
                        if not isinstance(dictSettings["RepaintFactor"],int):
                            strErr = "RepaintFactor must be an integer, not a %s"%str(type(dictSettings["RepaintFactor"]))
                            raise ValueError(strErr)
                        self._repaintFactor = dictSettings["RepaintFactor"]
                else:
                    if not isinstance(dictSettings["UnlimitedRescan"],bool):
                        strErr = "UnlimitedRescan in scanning path setting must be a bool, not a %s"%str(type(dictSettings["UnlimitedRescan"]))
                        raise ValueError(strErr)
                    if dictSettings["UnlimitedRescan"]:
                    #if infiniteRepaintFactor:
                        self._repaintFactor = maxRepFactor
                    else:
                        if "RepaintFactor" not in dictSettings.keys():
                            raise ValueError("RepaintFactor is not in the scanning path settings.")
                        else:
                            if not isinstance(dictSettings["RepaintFactor"],int):
                                strErr = "RepaintFactor must be an integer, not a %s"%str(type(dictSettings["RepaintFactor"]))
                                raise ValueError(strErr)
                            self._repaintFactor = dictSettings["RepaintFactor"]
                            
                    #verify volume repainting  for 'ScaledRepaint' or  'FullWeightRepainting'
                if "VolumeRepainting" not in dictSettings.keys():
                    raise ValueError("VolumeRepainting is not in the scanning path settings")
                if dictSettings["VolumeRepainting"] not in repaintingMethodVolume:
                    strErr = "Value for VolumeRepainting is not in %s : value: %s"%(str(repaintingMethodVolume),str(dictSettings["VolumeRepainting"]))
                    raise ValueError(strErr)
                self._repaintingMethodVolume = dictSettings["VolumeRepainting"]
                
            ######## End if  'IsolayerRepaint' else ( 'ScaledRepaint' or  'FullWeightRepainting')
            
        else: # If the user don't want to use the repainting
            self._repaintingMethodSlice = repaintingMethodSlice[1]
            self._repaintingMethodVolume = repaintingMethodVolume[1]
            self._repaintFactor = 1
        
        if "AddMapMargins" in dictSettings.keys():
            #raise ValueError("AddMapMargins is not in the scanning path settings")
            if not isinstance(dictSettings["AddMapMargins"],bool):
                strErr = "AddMapMargins in scanning path setting must be a bool, not a %s"%str(type(dictSettings["AddMapMargins"]))
                raise ValueError(strErr)
            self._addMapMargins = dictSettings["AddMapMargins"]
            if self._addMapMargins:
                self._arraySpotType = None
                if "MarginParams" not in dictSettings.keys(): 
                    raise ValueError("MarginParams is not in the scanning path settings")
                if not isinstance(dictSettings["MarginParams"],list):
                    strErr = "MarginParams in scanning path setting must be a list, not a %s"%str(type(dictSettings["MarginParams"]))
                    raise ValueError(strErr)
                
                checkList = [isinstance(item,(int, long, float)) for item in dictSettings["MarginParams"]]
                if not all(item == True for item in checkList):
                    raise ValueError("All elements in MarginParams should be numbers...")
                
                if "3DMargins" not in dictSettings.keys():
                    raise ValueError("3DMargins is not in the scanning path settings")
                if not isinstance(dictSettings["3DMargins"],bool):
                    strErr = "3DMargins in scanning path setting must be a bool, not a %s"%str(type(dictSettings["3DMargins"]))
                    raise ValueError(strErr)
                self._map3DMargin = dictSettings["3DMargins"]
                if self._map3DMargin:
                    if len(dictSettings["MarginParams"]) != 3:
                        strErr = "MarginParams should contain 3 parameters since 3DMargins is True. It contains %i parameters."%len(dictSettings["MarginParams"])
                        raise ValueError(strErr)
                    self._marginParams = dictSettings["MarginParams"]# extension in mm for [x , y, z] -> this means for ex. add a margin of (x mm) at the right of the map and (x mm) at the left of the map.
                else:
                    if len(dictSettings["MarginParams"]) > 3 or len(dictSettings["MarginParams"]) < 3:
                        strErr = "MarginParams should contain 2 or 3 parameters. It contains %i parameters."%len(dictSettings["MarginParams"])
                        raise ValueError(strErr)
                    self._marginParams = [0,0,0]
                    self._marginParams[0:2] = dictSettings["MarginParams"][0:2]# extension in mm for [x , y, z] -> this means for ex. add a margin of (x mm) at the right of the map and (x mm) at the left of the map.
                if "UpdateMargins" not in dictSettings.keys():
                    raise ValueError("UpdateMargins is not in the scanning path settings")
                if not isinstance(dictSettings["UpdateMargins"],bool):
                    strErr = "UpdateMargins in scanning path setting must be a bool, not a %s"%str(type(dictSettings["3DMargins"]))
                    raise ValueError(strErr)
                self._UpdateDMargin = dictSettings["UpdateMargins"]

        else:
            self._addMapMargins =  False
            
        if "FindStartPos" in dictSettings.keys():
            if not isinstance(dictSettings["FindStartPos"],bool):
                    strErr = "FindStartPos in scanning path setting must be a bool, not a %s"%str(type(dictSettings["FindStartPos"]))
                    raise ValueError(strErr)
            self._findStartPos = dictSettings["FindStartPos"]
        else:
            self._findStartPos = False
                
            #Arrays
        self._arraySpots = None
        self._initialPath = None
        self._arrayX = None
        self._arrayY = None
        self._arrayNumberOfRepainting = None
        
            #Variables
        self._gantryAngle = None
        self._nbEnergyLayers = None
        self._boundingBox = None
        self._spacing = None
        self._listSpotSize = None
        self._listEnergies = None
        self._idxInScanOrder = None
        
        print "--------------------------------------------"
        print "Scanning delivery settings:"
        print "\t\tType Scanning: %s"%str(self._scanPath)
        print "\t\tRepainting: %s"%str(dictSettings["Repainting"])
        if hasattr(self, '_repaintingMethodSlice'):
            print "\t\tSlice repaint: %s"%str( self._repaintingMethodSlice)
        if hasattr(self, '_repaintingMethodVolume'):
            print "\t\tVolume repaint: %s"%str( self._repaintingMethodVolume)
            #self._repaintFactor
        if hasattr(self, '_repaintFactor'):
            if self._repaintFactor == maxRepFactor:
                print "\t\tRepainting factor: %s"%'Unlimited'
            else:
                print "\t\tRepainting factor: %i"%self._repaintFactor
        print "\t\tMargins added: %s"%str(self._addMapMargins)
        if self._addMapMargins:
            print "\t\tMargin parameters (mm): %s"%str(self._marginParams)
        print "\t\tFind a good starting position: %s"%str(self._findStartPos)
        print "--------------------------------------------"


   
            

    def getRepaintingFactor(self):
        '''
        :returns: If a repainting factor exists, it returns it, otherwise it returns None.
    
        '''
        if hasattr(self, '_repaintFactor'):
            return self._repaintFactor
        return None
    
    @property
    def repaintingMethodSlice(self):
        '''
        :returns: The energy layer repainting method : 'IsolayerRepaint','ScaledRepaint' or 'FullWeightRepainting'
    
        '''
        return self._repaintingMethodSlice
        
    @property
    def repaintingMethodVolume(self):
        '''
        :returns: The volume repainting method : 'Volumetric' or 'NonVolumetric'
    
        '''
        return self._repaintingMethodVolume
        
    @property
    def arrayX(self):
        '''
        :returns: 3D numpy array containing X coordinate of scan spot positions. Each frame of the matrix correspond to an energy layer. \
        X value is set to 0 in the matrix if the spot is not a scan spot position.
        '''    
        return self._arrayX

    @property
    def arrayY(self):
        '''
        :returns: 3D numpy array containing Y coordinate of scan spot positions. Each frame of the matrix correspond to an energy layer. \
        Y value is set to 0 in the matrix if the spot is not a scan spot position.
        '''  
        return self._arrayY
        
    @property
    def arraySpots(self):
        '''
        :returns: 3D numpy array containing the planned weights of scan spot positions. Each frame of the matrix correspond to an energy layer. \
        Weight value is set to 0 in the matrix if the spot is not a scan spot position or if it has already been delivered.
        '''  
        return self._arraySpots

    @property
    def listEnergies(self):
        '''
        :returns: A list of energies used. The order is from higher energy to lowest.
        ''' 
        return self._listEnergies

    
    def __len__(self):
        '''
        :returns: Number of energy layers.
        ''' 
        return self._nbEnergyLayers
        

    def __iter__(self):
        ''' Iterates through the spots following the order defined in dictSettings. 
        :returns: A tuple: \
        ( weight spot ,index spot in the 3D matrices : frame, row, column , x positions,y positions,energy,spot size along x and along y , \
        strEnergy : 'NewEnergy' if the energy is changing or 'SameEnergy' otherwise ,typeSpot = 1 (native spot) or 2 (spot from added margin) , repainting index)

        '''
        prevEnergy = None
        if self._scanPath == 'random':
            indexArray = np.arange(self._arraySpots.size)
            np.random.shuffle(indexArray)
            arrayShape = self._arraySpots.shape
            for index in indexArray:
                print ">>>>>>>>> New Spot <<<<<<<<<"
                ind = np.unravel_index(index,arrayShape)
                if prevEnergy is None:
                    prevEnergy = self._listEnergies[ind[0]]
                    strEnergy = 'NewEnergy'
                elif prevEnergy != self._listEnergies[ind[0]]:
                    prevEnergy = self._listEnergies[ind[0]]
                    strEnergy = 'NewEnergy'
                else:
                    strEnergy = 'SameEnergy'
                if self._addMapMargins:
                    typeSpot = self._arraySpotType[ind]
                else:
                    typeSpot = typeOrigin
                if hasattr(self, '_repaintingIndex' ) and self._repaintingIndex is not None:
                    repaintIndex = self._repaintingIndex[ind[0]]
                else:
                    repaintIndex = 0
                dataSpot = (self._arraySpots[ind] ,ind , self._arrayX[ind],self._arrayY[ind],self._listEnergies[ind[0]],self._listSpotSize[ind[0]] ,strEnergy,typeSpot , repaintIndex)
                if dataSpot is None:
                    print "Dataspot is None:"
                    print "self._arraySpots[ind]: %s"%str( self._arraySpots[ind] )
                    print "ind: %s"%str( ind )
                    print "self._arrayX[ind],self._arrayY[ind]: %s , %s"%(str( self._arrayX[ind]  ), str(self._arrayY[ind]))
                    print "self._listEnergies[ind[0]]: %s"%str( self._listEnergies[ind[0]] )
                    print "self._listSpotSize[ind[0]]: %s"%str( self._listSpotSize[ind[0]] )
                    print "strEnergy: %s"%str( strEnergy )
                    print "typeSpot : %s"%str( typeSpot )
                    print "repaintIndex: %s"%str( repaintIndex )
                    raise ValueError ("Data spot is None for random scanning path")
                yield dataSpot
                
                
        
        elif self._scanPath == 'Regular':
            print "Regular delivery"
            arrayShape = self._arraySpots.shape
            it = self.nextSpotIndex()
            
            for idx,ind in enumerate(it):
                print ">>>>>>>>> New Spot <<<<<<<<<"
                print "Spot idx. in Scan Path: %i , (f,r,c):%s"%(idx,str(ind))
                if prevEnergy is None:
                    prevEnergy = self._listEnergies[ind[0]]
                    strEnergy = 'NewEnergy'
                elif prevEnergy != self._listEnergies[ind[0]]:
                    prevEnergy = self._listEnergies[ind[0]]
                    strEnergy = 'NewEnergy'
                else:
                    strEnergy = 'SameEnergy'
                dataSpot = self.buildDataSpotToReturn(ind , strEnergy)
                if dataSpot is None:
                    print "Dataspot is None:"
                    print "self._arraySpots[ind]: %s"%str( self._arraySpots[ind] )
                    print "ind: %s"%str( ind )
                    print "self._arrayX[ind],self._arrayY[ind]: %s , %s"%(str( self._arrayX[ind]  ), str(self._arrayY[ind]))
                    print "self._listEnergies[ind[0]]: %s"%str( self._listEnergies[ind[0]] )
                    print "self._listSpotSize[ind[0]]: %s"%str( self._listSpotSize[ind[0]] )
                    print "strEnergy: %s"%str( strEnergy )
                    print "typeSpot : %s"%str( typeSpot )
                    print "repaintIndex: %s"%str( repaintIndex )
                    raise ValueError ("Data spot is None for regular scanning path")
                yield dataSpot

                        
        else:
            raise ValueError('Unknown spotOrder in scanningPath')






        



    def nextSpotIndex(self):
        '''
        
        :yields: The next spot positions index (frame, rows, column in the 3D matrices representing the scanning path) if 'TypeScanning' (in 'dictSettings') \
        is not 'random'. This method is called from __iter__()
        
        '''
    
        if self._repaintingMethodVolume == repaintingMethodVolume[1]: #'Non Volumetric'
            for ind in self.repaintingNonVolOrder():
                yield ind
        elif self._repaintingMethodVolume == repaintingMethodVolume[0]: #'Volumetric'
            for ind in self.repaintingVolOrder():
                yield ind
        else:
            strErr = "Wrong repaintingMethodVolume: %s"%(self._repaintingMethodVolume)
            raise ValueError(strErr)
    
    
    
    def setNextSpotIndex(self, ind3D):
        '''Set the next spot index : frame, rows, column in the 3D matrices representing the scanning path
        
        :param ind3D: list of 3 elements : [frame,row,column]
        
        '''
        self._nextIndex =  [ind3D[0],ind3D[1],ind3D[2]]
    

    def buildSliceScanningOrder(self):
        '''Build the spot scanning order within the current energy slice. It sets the class attribute 'self._sliceScanningOrder'.
        
        '''
        if not self._addMapMargins:
            ind = self.getSliceOrderFromPlan()
            self._sliceScanningOrder = [ list(ind[0]),list(ind[1]) ]
            
        else:
            ind = self.getSliceOrderFromPlan()
            if len(ind[0]) == 0: #If nothing left to deliver initialize the variable and the iterations will stop
                self._sliceScanningOrder = [ list(ind[0]),list(ind[1]) ]
            else:
                if self._UpdateDMargin:
                    self._arraySpotType = np.zeros(self._arraySpots.shape,dtype='int8',order='C')
                    self.extendMap(self._arraySpots.shape, updateXYAndNbRepatings = False)

                ind = np.where( ((self._arraySpotType[self._currEnergySlice,:,:] == typeOrigin) | \
                                (self._arraySpotType[self._currEnergySlice,:,:] == typeAdded)) & \
                                (self._arrayNumberOfRepainting[self._currEnergySlice,:,:] > 0) )
                self._sliceScanningOrder = [ list(ind[0]),list(ind[1]) ]
        

    ## Repainting
    def repaintingNonVolOrder(self):
        ''' Generator for the non volumetric repainting order.\
        Order: slice by slice : repaint a slice before going to the next slice. \
        Higher energy to lowest energy.
        
        :yields: The next spot index (frame,row,column). This method is called by nextSpotIndex()
        
        '''
        print "Non Volumetric repainting"
        while True:
            if self._currEnergySlice == None:
                self._currEnergySlice = 0
            else:
                self._currEnergySlice = self._currEnergySlice + 1
                if self._currEnergySlice > self._arraySpots.shape[0]-1:
                    self._currEnergySlice = None
                    self._nextIndex = None
                    break
            self.buildSliceScanningOrder()
            print "Nb spots in scanning path: %i"%len(self._sliceScanningOrder[0])
            self._repaintingIndex[self._currEnergySlice] = self._repaintingIndex[self._currEnergySlice] + 1
            self._idxInScanOrder = 0
            while True:
                if self._idxInScanOrder < len(self._sliceScanningOrder[0]) :
                    print "Index in scanning path: %i"%self._idxInScanOrder
                    self.setNextSpotIndex([self._currEnergySlice, self._sliceScanningOrder[0][self._idxInScanOrder],self._sliceScanningOrder[1][self._idxInScanOrder]])
                    self._idxInScanOrder = self._idxInScanOrder + 1
                    yield self._nextIndex
                else:
                    self.buildSliceScanningOrder()
                    self._idxInScanOrder = 0
                    if self._idxInScanOrder < len(self._sliceScanningOrder[0]):
                        print "------------------Starting new repainting---------------"
                        self._repaintingIndex[self._currEnergySlice] = self._repaintingIndex[self._currEnergySlice] + 1
                        self.setNextSpotIndex([self._currEnergySlice, self._sliceScanningOrder[0][self._idxInScanOrder],self._sliceScanningOrder[1][self._idxInScanOrder]])
                        self._idxInScanOrder = self._idxInScanOrder + 1
                        yield self._nextIndex
                    else:
                        break
                        
    def repaintingVolOrder(self):
        '''Generator for the volumetric repainting order.
        
        :yields: The next spot index (frame,row,column). This method is called by nextSpotIndex()
        
        '''
        while True:
            if self._currEnergySlice == None:
                self._currEnergySlice = 0
                
            else:
                self._currEnergySlice = self._currEnergySlice + 1
                if self._currEnergySlice > self._arraySpots.shape[0]-1:
                    ind = np.where(self._arrayNumberOfRepainting > 0)
                    if len(ind[0]) != 0:
                        self._currEnergySlice = ind[0][0]
                        self._nextIndex = None
                    else:
                        self._currEnergySlice = None
                        self._nextIndex = None
                        break
            self.buildSliceScanningOrder()
            self._repaintingIndex[self._currEnergySlice] = self._repaintingIndex[self._currEnergySlice] + 1
            self._idxInScanOrder = 0
            while True:
                if self._idxInScanOrder < len(self._sliceScanningOrder[0]) :
                    print "------------------Starting new repainting---------------"
                    self.setNextSpotIndex([self._currEnergySlice, self._sliceScanningOrder[0][self._idxInScanOrder],self._sliceScanningOrder[1][self._idxInScanOrder]])
                    self._idxInScanOrder = self._idxInScanOrder + 1
                    yield self._nextIndex
                else:
                    break
                        
    
    def getWeightForSpotInd(self,ind3D):
        '''
        :param ind3D: 3D indices: (frame, row, column) defined in the 3D matrices representing the scanning path.
        
        :returns: The weight to be delivered of the spot stored at (frame, row, column) of the weight matrix. It takes into account the repainting index.
        
        '''
        weight = self._arraySpots[ind3D[0],ind3D[1],ind3D[2]]
        repaintFactor = self._arrayNumberOfRepainting[ind3D[0],ind3D[1],ind3D[2]]
        if self._repaintingMethodSlice == repaintingMethodSlice[0]:
            if repaintFactor >  1 :
                return self._maxWeight
            else:
                return weight
        elif self._repaintingMethodSlice == repaintingMethodSlice[1]:
            if repaintFactor >=  1 :
                return weight/np.float32(repaintFactor)
            else:
                return 0
        elif self._repaintingMethodSlice == repaintingMethodSlice[2]:
            return weight
                
                
    def getPlannedWeightForSpotInd(self,ind3D):
        '''
        :param ind3D: 3D indices: (frame, row, column) defined in the 3D matrices representing the scanning path.
        
        :returns: The planned weight of the spot stored at (frame, row, column) of the weight matrix.
        
        '''
        weight = self._initialPath[ind3D[0],ind3D[1],ind3D[2]]
        return weight

        
        
        
    def updateSpot(self, indices, weight, decrRepaint = True):
        '''Update the weight of a spot located in the 3D array at 'indices':
        remove 'weight' from the value contained in the matrix storing the updated spot weights.
        
        :param indices: 3D indices: (frame, row, column) defined in the 3D matrices representing the scanning path.
        :param weight: weight to be removed from the weight matrix
        :param decrRepaint: True / False : True if one wants to decrement the repainting factor of the current spot, False otherwise.    
        
        '''
        
        self._arraySpots[indices[0],indices[1],indices[2]] = self._arraySpots[indices[0],indices[1],indices[2]] - weight
        
        if (self._arraySpots[indices[0],indices[1],indices[2]] < 0):
            print  "Warning Hot Spot: dose deposited higher than planned"
        
        if self._arraySpots[indices[0],indices[1],indices[2]]  <= sys.float_info.epsilon:
            self._arraySpots[indices[0],indices[1],indices[2]]  = 0
           
        print "Spot updated: %s / new weight: %s"%( str(tuple((self._arrayX[indices[0],indices[1],indices[2]],self._arrayY[indices[0],indices[1],indices[2]])  )), str(float(self._arraySpots[indices[0],indices[1],indices[2]])))
        if decrRepaint:
            self.decrementRepaint(indices)
        
        
    def decrementRepaint(self , indices):
        '''Decrement the repainting index of the spot stored at 'indices' in the matrix storing the repainting factor of each spot.
        
        param indices: 3D indices: (frame, row, column) defined in the 3D matrices representing the scanning path.
        
        '''
        (f,r,c) = indices
        self._arrayNumberOfRepainting[f,r,c] = self._arrayNumberOfRepainting[f,r,c] - 1
    

    
    def getEnergySliceWeights(self, energy):
        '''Get the frame in which is stored all the spots positions for the energy layer 'energy'. 
        
        :param energy: energy 
        
        :returns: 2D numpy array containing the weights of all the spots positions. If the weight is 0, either the spot has been delivered or it was not \
        a spot position.
        
        '''   
        try:
            f = self._listEnergies.index(energy)
        except ValueError:
            strErr =  "Remove spot - Energy: %f , not in list: %s"%(energy,str(self._listEnergies))
            raise ValueError(strErr)      
        return self._arraySpots[f,:,:]
    
    ###############################################################
    #
    #               Methods used for compensation
    #
    #
    ###############################################################
    
    
    def isPointNearOtherSpot(self,energy,xShift, yShift):
        '''Find wether a point (x,y coordinates in the IEC Gantry coordinate system) has a scan spot position in its vicinity ( +/- 1 voxel in x and y directions) in the 3D matrices representing the scanning path.
        
        :param energy: beam energy
        :param xShift: x spot position in mm
        :param yShift: y spot position in mm
            
        :returns: Wether there is a scan spot in the neighborhood of +/- 1 voxel in an energy layer. If there is no scan spot, it returns tuple([False]), otherwise it returns a tuple:
            
            * index 0: True
            * index 1: 3D indices of the scan spot found
            * index 2: (x,y) coordinates of the scan spot found
            * index 3: the euclidean distance between the 2 spots
           
        '''
        shape = self._arraySpots.shape
        (r ,c ) = coordinatesToIndex([xShift] , [yShift],   shape[1], shape[2],self._spacing )
        
        if (r >= shape[1]) or ( c >= shape[2]) or ( r < 0) or (c < 0):
            print "Spot (%f,%f) -> (%i,%i) is outside the plan matrix (%i,%i)"%(xShift , yShift, r,c,shape[1],shape[2])
            return tuple([False])
        
        try:
            f = self._listEnergies.index(energy)
        except ValueError:
            strErr =  "Remove spot - Energy: %s , not in list: %s"%(str(energy),str(self._listEnergies))
            raise ValueError(strErr)
        spotsInLayer = np.where((self._arrayNumberOfRepainting[f,:,:] > 0) & \
                            (self._arraySpots[f,:,:] > 0) )
        
        rows = spotsInLayer[0]
        cols = spotsInLayer[1]
        listSpots = []
        tmpSLiceEnergy = np.zeros((shape[1],shape[2]))
        tmpSLiceEnergy[r,c] = 0.5
        print "Spacing used in scanning path (mm): %s"%str(self._spacing)
        minDist = np.Inf
        xyClosest = None
        closestInd = None
        for ro,co in zip(rows,cols):
            x = self._arrayX[f, ro ,co]
            y = self._arrayY[f, ro ,co]
            dist = euclideanDistance( (xShift, yShift) , (x,y))            
            if dist < minDist and np.abs(x-xShift) < (self._spacing[0]/2.0) and np.abs(y-yShift) < (self._spacing[1]/2.0):
#             if dist < minDist and np.abs(x-xShift) < (self._spacing[0]) and np.abs(y-yShift) < (self._spacing[1]):
                minDist = dist
                xyClosest = (x,y)
                closestInd = (ro , co)
        
        if xyClosest is not None and closestInd is not None:
            print "Closest spot found at : (%f,%f)(%s) w:%f - and original spot was (%f,%f)(%s). Actual distance is:%f "%(xyClosest[0],xyClosest[1],str(closestInd),self._arraySpots[f,closestInd[0] ,closestInd[1]],xShift,yShift,str((int(r),int(c))),minDist)
            tmpSLiceEnergy[closestInd[0] ,closestInd[1]] = 1
            return (True , (f, closestInd[0] ,closestInd[1]),xyClosest,minDist)
            
        else:
            return tuple([False])
 





    def findStartingPointForMotion( self, dispVec_IECG ,ddCurvesData ,strEnergy ):
        '''Find a starting point for the delivery of the current energy slice, i.e. find a point that could be compensated as soon as the delivery\
        of the current energy slice start. 
        
        .. note::
        
            **Compensation principle** :
            
                We assume the current pencil beam position being (x_cur,y_cur). This spot can be compensated\
                if, for a given displacement vector ( dx,dy) , there exists a scan spot position whose weight is >0 in the neighborhood (+/- 1 voxel in\
                the 2D matrices representing a scanning path of the current energy layer) of the position \
                ( x_cur - dx , y_cur - dy) of the current energy layer.
        
        :param dispVec_IECG: displacement vector in the IEC gantry coordinate system in mm.
        :param ddCurvesData: depth dose curves data (dose profiles as a function of depth) for all the energies available. \
        (To be used for further implementation of 3D compensation)
        :param strEnergy: 'NewEnergy' or 'SameEnergy' depending if this is a repainting of the first painting of the energy layer. This is used to return\
        all the needed information to the calling function.
        
        
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
        if self._findStartPos:
            if self._idxInScanOrder != 1:
                strErr = "Error - in findStartingPointForMotion, _idxInScanOrder should be 1 but it is %i!"%self._idxInScanOrder
                raise ValueError(strErr)
            (dx,dy,dz) = dispVec_IECG # in mm
            threshold = np.sqrt( (self._spacing[0])*(self._spacing[0]) + (self._spacing[1])*(self._spacing[1]))
            startPos = None
            for (r0,c0) in zip(self._sliceScanningOrder[0] , self._sliceScanningOrder[1]): 
                energy = self._listEnergies[self._currEnergySlice]
                xShift = self._arrayX[self._currEnergySlice,r0,c0] - dx
                yShift = self._arrayY[self._currEnergySlice,r0,c0] - dy
                energyForCompens = energy
                print "Ener for compens: %s"%str(energyForCompens)
                newSpot = self.isPointNearOtherSpot(energyForCompens,xShift, yShift)
                if newSpot[0] and newSpot[3] < threshold:
                    startPos = list((self._currEnergySlice,r0,c0))
                    print "Found a good starting position: %s can be compensated by %s"%( str(startPos),str(newSpot[1]))
                    break
            if startPos is None or (startPos is not None and startPos == self._nextIndex): 
                print "No good starting spot found"
                return None
            flagIdxFound = False
            for idx , (r,c) in enumerate(zip(self._sliceScanningOrder[0],self._sliceScanningOrder[1])):
                if r == startPos[1] and c == startPos[2]:
                    flagIdxFound = True
                    self._idxInScanOrder = idx+1
                
            if not flagIdxFound:
                strErr = 'Error - no index found in _sliceScanningOrder for (%i,%i)'%(startPos[1],startPos[2])
                raise ValueError(strErr)
            dataSpot = self.buildDataSpotToReturn(startPos , strEnergy)
            self.setNextSpotIndex(startPos)
            return dataSpot
        else:
            return None 
 

    def findCorrespondingEnergy(self, newEnergy):
        '''Find in the scanning path the closest energy to newEnergy.
        
        :param newEnergy: input energy
        
        :returns: None if no energy has been found, return the found energy otherwise.
        
        '''
        
        cpyEnergies0 = self._listEnergies[0:-1]
        cpyEnergies1 = self._listEnergies[1:]
        energyFound = None
        for e0 , e1 in zip(cpyEnergies0,cpyEnergies1):
            if e0 >= newEnergy >= e1:
                if np.abs(e0 - newEnergy) < np.abs(e1 - newEnergy):
                    energyFound = e0
                else:
                    energyFound = e1
                break
        if energyFound is None:
            if np.floor(self._listEnergies[0]) <= newEnergy <= np.ceil(self._listEnergies[0]):
                energyFound = self._listEnergies[0]
            elif np.floor(self._listEnergies[-1]) <= newEnergy <= np.ceil(self._listEnergies[-1]):
                energyFound = self._listEnergies[-1]
            else:
                print "No energy found in scanning path corresponding to %s MeV"%str(newEnergy)
        
        if energyFound is not None:
            print "Energy found in scanning path: %f m for energy: %f"%(energyFound,newEnergy)
        return energyFound
    

    
    
    

  
    ###############################################################
    #
    #               Methods used for scanning path initialization
    #
    #
    ###############################################################
  
    
    def buildScanningPath(self , rpData, beamIdx = 0):
        '''Method that build the scanning path and creates initialize all the attributes (including the set of matrices) of this class.
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)
        
        
        '''
        self._beamIdx = beamIdx
        self._rpData = rpData
        self.findBoundingBoxScanningPath(rpData, beamIdx)
        self.fillSpotSizeList( rpData, beamIdx)
        self.fillListEnergies(rpData, beamIdx)
        self.findSpacing(rpData, beamIdx)
        if self._addMapMargins:
            self.extendBoundingBox()
        
        arrayShape = self.findShape3DArray()
        
        self._arraySpots = np.zeros(arrayShape,dtype=self._typeFloat,order='C')
        self._initialPath = np.zeros(arrayShape,dtype=self._typeFloat,order='C')
        self._arrayX = np.zeros(arrayShape,dtype=self._typeFloat,order='C')
        self._arrayY = np.zeros(arrayShape,dtype=self._typeFloat,order='C')
        
        self._arrayNumberOfRepainting = np.zeros(arrayShape,dtype=typeIntRepaint,order='C')
        self._repaintingIndex = -1*np.ones(arrayShape[0],dtype=typeIntRepaint,order='C')
        self.fillArraysScanningPath(rpData, beamIdx)
        
        self._initialPath[:,:,:] = self._arraySpots[:,:,:]
    
        self.fillRepaintingArray()
        if self._addMapMargins:
            self._arraySpotType = np.zeros(arrayShape,dtype='int8',order='C')
            self.extendMap(arrayShape, updateXYAndNbRepatings = True)
        
            
        
        print "Array spots (min MU,max MU,nb spots): \n%s"%str((np.min(self._arraySpots),np.max(self._arraySpots),len(np.where(self._arraySpots > 0)[0]) ))
        print "Nb spots: %i"%len(np.where(self._arraySpots > 0)[0])
        print "Scanning Path ... ok"
        
        
    def fillArraysScanningPath(self,rpData, beamIdx):
        '''Method that fills the matrices storing X and Y coordinates and the weights from the treatment plan.
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)
        
        '''
        
        beam = rpData.IonBeamSequence[beamIdx]
        shape = self._arraySpots.shape
        count = 0
        (i,j) = coordinatesToIndex([self._boundingBox[0][0],self._boundingBox[1][0]] , [self._boundingBox[0][1],self._boundingBox[1][1]],   shape[1], shape[2],self._spacing )
        for energyLayer in beam.IonControlPointSequence:
                if not isEnergyLayerEmpty(energyLayer):
                    listX = energyLayer.ScanSpotPositionMap[0::2]
                    listY = energyLayer.ScanSpotPositionMap[1::2]
                    (j ,i ) = coordinatesToIndex(listX , listY,   shape[1], shape[2],self._spacing )
#                     print (j ,i )
                    self._arraySpots[count,j,i] =  energyLayer.ScanSpotMetersetWeights 
                    self._arrayX[count,j,i] = listX
                    self._arrayY[count,j,i] = listY
                    count = count + 1
                    
                
    def fillRepaintingArray(self):
        '''Method that fills the matrice storing the repainting factor for each scan spot.
                
        '''        
            
        if self._repaintingMethodSlice == repaintingMethodSlice[0]:
#             self._maxWeight = np.mean(self._arraySpots[np.where(self._arraySpots > 0) ])
            self._maxWeight = np.min(self._arraySpots[np.where(self._arraySpots > 0) ])
            print "Isolayer repainting : max weight = %f"%self._maxWeight
            self._arrayNumberOfRepainting[:,:,:] = np.round(self._arraySpots[:,:,:]/self._maxWeight)
            print "max repaint factor : %i"%np.max(self._arrayNumberOfRepainting)
        elif self._repaintingMethodSlice == repaintingMethodSlice[1] or \
             self._repaintingMethodSlice == repaintingMethodSlice[2]:
            print "Scaled repainting - repainting factor = %i"%self._repaintFactor
            indSpots = np.where(self._arraySpots > 0)
            self._arrayNumberOfRepainting[indSpots] = self._repaintFactor
        else:
            strErr = "in fillRepaintingArray() unknown type %s "%str(self._repaintingMethodSlice)
            raise ValueError (strErr)
        
        
    def findBoundingBoxScanningPath(self, rpData, beamIdx):
        '''The goal of this function to find the number of energy layers, 
        and find xmin, xmax ymin ymax of the scan spots, i.e. a bounding box
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)

        
        '''
        beam = rpData.IonBeamSequence[beamIdx]
        self._gantryAngle = beam.IonControlPointSequence[0].GantryAngle
        minX = min(beam.IonControlPointSequence[0].ScanSpotPositionMap[0::2])
        maxX = max(beam.IonControlPointSequence[0].ScanSpotPositionMap[0::2])
        minY = min(beam.IonControlPointSequence[0].ScanSpotPositionMap[1::2])
        maxY = max(beam.IonControlPointSequence[0].ScanSpotPositionMap[1::2])
        count = 1

        for energyLayer in beam.IonControlPointSequence[1:]:
            if not isEnergyLayerEmpty(energyLayer):
                count = count + 1
                listX = energyLayer.ScanSpotPositionMap[0::2]
                listY = energyLayer.ScanSpotPositionMap[1::2]
                currMinX = min(listX)
                currMaxX = max(listX)
                currMinY = min(listY)
                currMaxY = max(listY)
                if currMinX < minX:
                    minX = currMinX
                if currMaxX > maxX:
                    maxX = currMaxX
                if currMinY < minY:
                    minY = currMinY
                if currMaxY > maxY:
                    maxY = currMaxY
        self._nbEnergyLayers = count
        #print "Bounding box: xmin:%f , xmax:%f , ymin: %f , ymax: %f"%(minX,maxX,minY,maxY)                    
        print "Nb energy layers: %i"%self._nbEnergyLayers
        self._boundingBox = [( minX-0.01,minY-0.01) , ( maxX+0.01,maxY+0.01)] 

    def extendBoundingBox(self):
        '''This method extends the bounding box of the initial scanning path if 'AddMapMargins' (in dictSettings) is True. \
        If 'MarginParams' = [dx,dy,dz], it is extended of dx at the left and the right and dy at\
        the top and bottom of the initial bounding box
        
        '''
        [( minX,minY) , ( maxX,maxY)] = self._boundingBox  
        [dx,dy,dz] = self._marginParams
        self._boundingBox = [( minX-dx,minY-dy) , ( maxX+dx,maxY+dy)]


    def fillSpotSizeList(self , rpData, beamIdx):
        '''This function fills the list of spot sizes (attribute of the class).
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)
      
        '''
        beam = rpData.IonBeamSequence[beamIdx]
        self._listSpotSize = list()
        for energyLayer in beam.IonControlPointSequence:
            if not isEnergyLayerEmpty(energyLayer):
                #sig = FWHM / (2 * sqrt( 2 * ln 2)) ~ FWHM / 2.3548200 
                #self._listSpotSize.append([value/2.3548200 for value in energyLayer.ScanningSpotSize])
                self._listSpotSize.append([value for value in energyLayer.ScanningSpotSize])


    def fillListEnergies(self , rpData, beamIdx):
        '''This function fills the list of energies (attribute of the class).
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)

        '''
        beam = rpData.IonBeamSequence[beamIdx]
        self._listEnergies = list()
        for energyLayer in beam.IonControlPointSequence:
            if not isEnergyLayerEmpty(energyLayer):
                self._listEnergies.append(energyLayer.NominalBeamEnergy)
        
    
    def getListEnergies(self):
        '''Get the list of energies used in the scanning path.
  
        '''        
        
        return self._listEnergies
    
    def findSpacing(self, rpData, beamIdx):
        '''This function finds the spacing to use to represent the scanning path spatially in a 3D matrix. The spacing is the smallest y (resp. x) difference between\
        two consecutive scan spots.
        
        :param rpData: RP dicom data from pydicom
        :param beamIdx: Beam index (i.e. gantry angle index)

        '''    
    
    
        minX = np.finfo(np.float32).max
        minY = np.finfo(np.float32).max
        beam = rpData.IonBeamSequence[beamIdx]
        
        for energyLayer in beam.IonControlPointSequence:
                if not isEnergyLayerEmpty(energyLayer):
                    if len(energyLayer.ScanSpotPositionMap) == 2: #Only 1 spot in the energy layer
                        continue
                    
                    argMinX = None
                    argMinY = None
                    listX = energyLayer.ScanSpotPositionMap[0::2] #X Position in mm in the IEC Gantry coordinate system
                    x1 = np.array(listX[0:-1])
                    x2 = np.array(listX[1:])
                    diffX = np.abs( (x1-x2) )
                    currMinx = np.min(diffX )
                    
                    if currMinx == 0:
                        argMinX = np.argmin(diffX)# Used to know if diffX = 0 somewhere
                        ind = np.where(diffX > 0)
                        if len(ind[0]) > 0: 
                            currMinx = np.min(diffX[ind])
                        else:
                            currMinx = 0  # really small value greater than 0 -> avoid to divide by zero when spacing = 0
                    if currMinx < minX and currMinx > 0:
                        minX = currMinx
                    
                    
                    
                    listY = energyLayer.ScanSpotPositionMap[1::2]#Y Position in mm in the IEC Gantry coordinate system
                    y1 = np.array(listY[0:-1])
                    y2 = np.array(listY[1:])
                    diffY = np.abs( (y1-y2) )
                    currMiny = np.min( diffY ) 
        
                    if currMiny == 0:
                        argMinY = np.argmin(diffY)# Used to know if diffY = 0 somewhere
                        ind = np.where(diffY > 0)
                        if len(ind[0]) > 0: 
                            currMiny = np.min(diffY[ind])
                        else:
                            currMiny = 0  # really small value greater than 0 -> avoid to divide by zero when spacing = 0
                        
                    if currMiny < minY and currMiny >  0:
                        minY = currMiny
                    
                        
                    if argMinX is not None:
                        if diffY[argMinX] == 0:
                            raise ValueError("2 Spots have same X and Y in the plan")
                    
                    if argMinY is not None:
                        if diffX[argMinY] == 0:
                            raise ValueError("2 Spots have same X and Y in the plan")
                        
                            
                    
        if minX ==   np.finfo(np.float32).max:
            minX = 1
        if minY ==   np.finfo(np.float32).max:
            minY = 1        
        self._spacing = [minX , minY]
        print "Spacing scanning path (mm): %s"%str(self._spacing)
        
    
    def findShape3DArray(self):
        '''Function that finds the shape of the 3D matrices that are created to represent the scanning path.\
        It is based on the number of energy layers, the spacing and the bounding box.
        
        '''
        if self._nbEnergyLayers is not None and self._boundingBox is not None and self._spacing is not None:
            nbFrames = self._nbEnergyLayers
            maxAbsY = max([ abs(self._boundingBox[1][1]) , abs(self._boundingBox[0][1]) ])
            maxAbsX = max([ abs(self._boundingBox[1][0]) , abs(self._boundingBox[0][0]) ])
            nbRows = np.ceil(2 * ( maxAbsY + 1)/self._spacing[1])
            nbCols = np.ceil(2 * ( maxAbsX + 1)/self._spacing[0])

            if self._addMapMargins and self._map3DMargin:
                #Using PSTAR ranges in water for energies 30MeV to 250MeV, the average change of range is 1e-1 g.cm-2 per MeV which correspond to 1mm/MeV in water
                deltaRange = 1 # mm/MeV
                #self._listEnergies should be ordered in decreasing order
                deltaEnergy = np.mean( np.array(self._listEnergies[:-1]) - np.array(self._listEnergies[1:]) )
                dz = self._marginParams[2]
                dFrames = int(np.ceil( dz / (deltaRange * deltaEnergy) ))
                nbFrames = nbFrames + 2 * dFrames
                firstEner = self._listEnergies[0]
                lastEner = self._listEnergies[-1]
                for idx in range(dFrames):
                    self._listEnergies.insert(0, firstEner - (idx+1)*deltaEnergy )
                    self._listEnergies.append(lastEner + (idx+1)*deltaEnergy )
                if len(self._listEnergies) == nbFrames:
                    self._nbEnergyLayers = nbFrames
                else:
                    strErr = "Wrong number of frame when trying 3D margin in plan. nb layers = %i , nbFrames = %i , originally %i layers"%(len(self._listEnergies) , nbFrames,self._nbEnergyLayers)
                    raise ValueError(strErr)
                
            #self._spacing = [minSpacingX , minSpacingY]
            shape = (nbFrames , np.round(nbRows), np.round(nbCols))
            print "Shape array scanning path : %s"%str(shape)
            return shape
        else:
            print "in find 3D shape: \nself._nbEnergyLayers : %s \nself._boundingBox : %s \nself._listSpotSize: %s"%(str(self._nbEnergyLayers),str(self._boundingBox),str(self._listSpotSize))
            raise ValueError ("in find 3D shape: self._nbEnergyLayers is None or self._boundingBox is None or self._listSpotSize is None")


    def extendMap(self,shape, updateXYAndNbRepatings = True):
        '''Extend scan spot position maps in the case of margins added.
        
        :param shape: tuple (number of frames, number of rows, number of columns). Shape of the 3D matrices.
        :param updateXYAndNbRepatings: True / False to update the 3D matrices stroring X and Y coordinates and the number of repaintings.
        
        '''
        self.fillSpotTypeArray(shape)
        if updateXYAndNbRepatings:
            self.updateXYArrayForMargin(shape)
            self.updateNbRepaintingsForMargins()
     

    def fillSpotTypeArray(self, shape):
        '''Fill the 3D array of spot types:
            
            * typeOutside = 0
            * typeOrigin = 1
            * typeAdded = 2
            
        :param shape: tuple (number of frames, number of rows, number of columns). Shape of the 3D matrix.
        
        '''
        maskInitialMap = np.zeros(shape,dtype='int8',order='C')
        indMap = np.where(self._arraySpots > 0)
        maskInitialMap[indMap] = 1
        [dx,dy,dz]=self._marginParams
        nbIterations = int(np.max([ np.ceil(dx / self._spacing[0] ), np.ceil(dy / self._spacing[1] )]))
        dilatedMap = ndimage.binary_dilation(maskInitialMap,iterations = nbIterations).astype('int8')
        indAddedSpot = np.where( (maskInitialMap == 0) & ( dilatedMap == 1))
        self._arraySpotType[indMap] = typeOrigin
        self._arraySpotType[indAddedSpot] = typeAdded
    
    def updateXYArrayForMargin(self,shape):
        '''Update X, Y 3D matrices for margins added.
        
        :param shape: tuple (number of frames, number of rows, number of columns). Shape of the 3D matrix.

        
        '''
        indAddedSpots = np.where(self._arraySpotType == typeAdded)
        (listX,listY) = indexToCoordinates( indAddedSpots[1] , indAddedSpots[2] , shape[1], shape[2], self._spacing )
        self._arrayX[indAddedSpots] = listX
        self._arrayY[indAddedSpots] = listY
        
    def updateNbRepaintingsForMargins(self):
        '''Update 3D matrix storing the number of repaintings. To be used if margins are added.
        
        '''    
        indAddedSpots = np.where(self._arraySpotType == typeAdded)
        maxRepaint = np.max(self._arrayNumberOfRepainting)
        self._arrayNumberOfRepainting[indAddedSpots] = maxRepaint
        


    def buildDataSpotToReturn(self, ind , strEnergy):
        '''Build the spot data to return for "regular" delivery.
        It should return a tuple with:
        
        
        Inputs:
        ind: spot indices
        strEnergy: string saying if the is a new energy or not
        
        '''
        weight = self._arraySpots[ind[0],ind[1],ind[2]]
        x = self._arrayX[ind[0],ind[1],ind[2]]
        y = self._arrayY[ind[0],ind[1],ind[2]]
        ener =self._listEnergies[ind[0]]
        size = self._listSpotSize[ind[0]]
        if self._addMapMargins:
            typeSpot = self._arraySpotType[ind[0],ind[1],ind[2]]
        else:
            typeSpot = typeOrigin
        repaintFactor = self._arrayNumberOfRepainting[ind[0],ind[1],ind[2]]
        if hasattr(self, '_repaintingIndex' ) and self._repaintingIndex is not None:
            repaintIndex = self._repaintingIndex[ind[0]]
        else:
            repaintIndex = 0
        if self._repaintingMethodSlice == repaintingMethodSlice[0]:
            if repaintFactor >  1 :
                return( self._maxWeight , ind, x,y,ener,size ,strEnergy,typeSpot,repaintIndex)
            else:
               return( weight , ind, x,y,ener,size ,strEnergy,typeSpot,repaintIndex)
        elif self._repaintingMethodSlice == repaintingMethodSlice[1]:
            if repaintFactor >=  1 :
                return( weight/np.float32(repaintFactor) , ind, x,y,ener,size ,strEnergy,typeSpot,repaintIndex)
            #else:
                #yield( 0 , ind, x,y,ener,size ,strEnergy)
        elif self._repaintingMethodSlice == repaintingMethodSlice[2]:
                return( weight , ind, x,y,ener,size ,strEnergy,typeSpot,repaintIndex)
        else:
            strErr = "Repainting method is set to %s. It should be Isolayer or Scaled or full weight"%str(self._repaintingMethodSlice)
            raise ValueError(strErr)
        
        
    def getSliceOrderFromPlan(self):
        '''Get the regular scanning order (order provided in the dicom file) for the current energy slice.
        
        :returns: tuple of 2 numpy arrays containing respectively the row index and the column index of each spot.
        
        '''
        beam = self._rpData.IonBeamSequence[self._beamIdx]
        ind_rows_cols = [[],[]]
        for energyLayer in beam.IonControlPointSequence:
            if not isEnergyLayerEmpty(energyLayer):
                if energyLayer.NominalBeamEnergy == self._listEnergies[self._currEnergySlice]: 
                    listX = energyLayer.ScanSpotPositionMap[0::2]
                    listY = energyLayer.ScanSpotPositionMap[1::2]
                    for (x,y) in zip(listX,listY):
                        ind = np.where( (self._arrayX[self._currEnergySlice,:,:] == x) & (self._arrayY[self._currEnergySlice,:,:] == y) & \
                                        (self._arrayNumberOfRepainting[self._currEnergySlice,:,:] > 0) & \
                                        (self._arraySpots[self._currEnergySlice,:,:] > 0))
                        if len(ind[0]) > 0:
#                             print ind
                            ind_rows_cols[0].append(ind[0][0])
                            ind_rows_cols[1].append(ind[1][0])
                    return (np.array(ind_rows_cols[0]),np.array(ind_rows_cols[1]))
        return (np.array([]),np.array([]))

    
    def getNumberOfRepaintingsPerSlice(self):
        '''
        Returns the number of time each energy slice hase been repaint
        '''
        return self._repaintingIndex + 1
        
    

    ###############################################################
    #
    #               Methods used for visualization
    #
    #
    ###############################################################

    def locateSpots(self,listCoord):
        '''Create a 2D array where spots defined in listCoord are set to a non-null value. \
        This method is used to visually locate spots by displaying the array afterwards.

        :param listInd: list of coordinates of spots to display
        :param listCoord: list of [x,y] coordinates of the spots of interest.

        :returns: A 2D numpy array where value = 0 if this is not a spot of interest and  1>= value  >= 0.5 otherwise.

        '''
        shape = self._arraySpots.shape

        tmpSLiceEnergy = np.zeros((shape[1],shape[2]))
        nbSpots = len(listCoord[0])
        if nbSpots == 0:
            return None
        for i,(x,y) in enumerate(zip(listCoord[0],listCoord[1])):
            
            (r ,c ) = coordinatesToIndex([x] , [y],   shape[1], shape[2],self._spacing )
            (r,c) = (r[0] ,c[0] )
            print "Locate spots: %f %f , %i %i"%(x,y,r,c)
            tmpSLiceEnergy[r,c] = 0.5 + i * (0.5/nbSpots )
        return tmpSLiceEnergy
    
    
    def getMapWithEnhancedSpot(self,ind):
        '''Create a 2D array representing all the spot position of a energy layer where a given spot, defined by its 3D indices \
        in  the matrices representing the scanning path, is "highlighted". \
        This method is used to visually locate a spot by displaying the array afterwards.

        :param ind: indices (frame,row,column) of the spot of interest.
        

        :returns: A 2D numpy array where value = 0 if this is not a spot positions, value = 1 if this is a spot position in the energy layer and value = 2\
        for the spot of interest.

        '''        
        if self._addMapMargins:
            indSpots = np.where(self._arraySpotType[ind[0],:,:] > 0)
            array = np.zeros(self._arraySpotType[ind[0],:,:].shape, dtype = 'int8',order = 'C')
            array[indSpots] = 1
            array[ind[1],ind[2]] = 2
            return array
        else:
            indSpots = np.where(self._initialPath[ind[0],:,:] > 0)
            array = np.zeros(self._initialPath[ind[0],:,:].shape, dtype = 'int8',order = 'C')
            array[indSpots] = 1
            array[ind[1],ind[2]] = 2
            return array     
            
    def getMapWithMargin(self):
        '''Create a 2D array representing the positions added to the original scan spot position (margins added) and the original positions.\
        This method is used to visually locate the different regions.
        
        :returns: A 2D numpy array with 3 different areas depending on the spot types.
        
        '''
        if self._addMapMargins:
            if self._scanPath != scanningPath[0]:
                return self._arraySpotType[self._currEnergySlice,:,:]
        return None


###############################################################
#
#              Utilities
#
#
###############################################################


def euclideanDistance( xy1 , xy2):
    '''Computes the Euclidean distance between 2 points: 
    :param xy1: (x1, y1) 
    :param xy2: (x2,y2)
    
    :returns: Euclidean distance
    
    '''
    dist = np.sqrt((xy1[0]-xy2[0])* (xy1[0]-xy2[0]) + (xy1[1]-xy2[1])* (xy1[1]-xy2[1]))
    return dist

def isEnergyLayerEmpty(energyLayer):
    '''Check wether an energy layer in the dicom file is empty.
    
    :returns: True if empty, False otherwise.
    
    '''
    if isinstance(energyLayer.ScanSpotMetersetWeights,float) or isinstance(energyLayer.ScanSpotMetersetWeights,int):
        if energyLayer.ScanSpotMetersetWeights == 0:
            return True
    else:
        if sum(energyLayer.ScanSpotMetersetWeights) == 0:
            return True
    return False



def coordinatesToIndex( listX , listY , nRows, nCols, spacing ):
    '''Find indices in 2D.\
    For a given list of X coordinates, return the column index for each of them and for a given list of Y coordinates, return the row index for each of them.
    
    :param listX: list of x coordinates
    :param listX: list of y coordinates
    :param nRows: number of rows in a 2D matrix
    :param nCols: number of columns in a 2D matrix
    :param spacing: list :[xSpacing, y spacing]
    
    Computes lhe list of column indices and the list or rows indices. (x = 0 , y = 0) as being (f = (NumberOfFrames - 1)/2, r = (NumberOfRows - 1)/2, c = (NumberOfCols - 1)/2)
    
    :returns: (round_j  , round_i ), where i are columns and j rows. X is associated to columns and y to rows.
    
    '''
   
    x = np.array(listX)
    y = np.array(listY)

    i = ( x / spacing[0] ) + (nCols-1)/2.0 
    j = ( y / spacing[1] ) + (nRows-1)/2.0
    
    round_i = np.asarray(np.floor(i),dtype = 'int16')
    round_j = np.asarray(np.floor(j),dtype = 'int16')
    
        
    return ( round_j,  round_i)

def indexToCoordinates( listRow , listCols , nRows, nCols, spacing ):
    '''Find coordinates in 2D.\
    For a given list of Row indices, return the x coord for each of them and \
    for a given list of Column indices, return the y coord for each of them.
    
    :param listRow: list of rows indices
    :param listCols: list of columns indices
    :param nRows: number of rows in a 2D matrix
    :param nCols: number of columns in a 2D matrix
    :param spacing: list :[xSpacing, y spacing]

    
    :returns: The list of x coords and the list of y coords. 
    
    '''    
    r = np.array(listRow)
    c = np.array(listCols)

    x = ( c - ((nCols-1)/2.0)) * spacing[0]
    y = ( r - ((nRows-1)/2.0)) * spacing[1]     

    
    return ( x,  y)



