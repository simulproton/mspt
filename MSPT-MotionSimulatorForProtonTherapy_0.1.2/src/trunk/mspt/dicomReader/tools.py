########################################################################
#
# tools.py
# Dicom Reader Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# On Oct, 16 2012
# 
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
from PIL import Image , ImageDraw
import math


'''
Set of useful tools that can be used when using dicom readers.
'''

def tripletwise(iterable):
    '''Generate an triplets from iterable.
    
    :param iterable: An iterable
    
    :return: An iterable able to provide triplets of the input iterable.
    
    '''
    newIterator = iter(iterable)
    return zip(newIterator,newIterator,newIterator)



def coordinatesToIndexFromImagePosition(point, imagePosition, imageOrientation, pixelSpacing):
    '''Computes the indices (frame, row, column) , i.e. voxel location, for a given point (x,y,z) provided in an image coordinate system (mm).

    :param point: (x,y,z) coordinates of a point
    :param imagePosition: refers to coordinates (x0,y0,z0) of the top left corner in the image 
    :param imageOrientation: unit vectors in the 3 directions with respect to the patient coordinate system
    :param pixelSpacing: pixel spacing in the 3 directions
    
    Formula: look at page 410 from Dicom Documentation: "PS 3.3-2011 Digital Imaging and Communications in Medicine (DICOM) Part 3: Information Object Definitions"\
    More information on dicom standards can be found at Dicom `NEMA <http://medical.nema.org/standard.html>`_ and `here <http://medical.nema.org/Dicom/2011/11_03pu.pdf>`_ 
    
    :returns: (frame, row, column) indices

    '''
    X = [float(imageOrientation[i]) for i in range(3)]
    Y = [float(imageOrientation[i+3]) for i in range(3)]
    Z = [float(imageOrientation[i+6]) for i in range(3)]
    di = float(pixelSpacing['cols'])
    dj = float(pixelSpacing['rows'])
    dk = float(pixelSpacing['frames'])
    S = [float(imagePosition[i]) for i in range(len(imagePosition))]
    M = np.mat([[X[0]*di,Y[0]*dj,Z[0]*dk,S[0] ],\
                [X[1]*di,Y[1]*dj,Z[1]*dk,S[1] ],\
                [X[2]*di,Y[2]*dj,Z[2]*dk,S[2] ],\
                [0,0,0,1 ]])

    res = M.I * np.mat([[float(point[0])],[float(point[1])],[float(point[2])],[1.0]])
    return (int(np.round(res[2,0])),int(np.round(res[1,0])),int(np.round(res[0,0])))
    
def indexToCoordinatesFromImagePosition(index, imagePosition, imageOrientation, pixelSpacing):
    '''Computes the coordinates of a point (x, y, z)in an image coordinate system (mm) for a given point indices (pixel location) (frame, row,col).

    :param index: (frame, row,col) indices of  a voxel.
    :param imagePosition: refers to coordinates (x0,y0,z0) of the top left corner in the image 
    :param imageOrientation: unit vectors in the 3 directions with respect to the patient coordinate system
    :param pixelSpacing: pixel spacing in the 3 directions

    
    Formula: look at page 410 from Dicom Documentation: "PS 3.3-2011 Digital Imaging and Communications in Medicine (DICOM) Part 3: Information Object Definitions"\
    More information on dicom standards can be found at Dicom `NEMA <http://medical.nema.org/standard.html>`_
    
    :returns: (x, y, z) coordinates   
    
    '''
     
    index_float = [float(index[i]) for i in range(len(index))]
    X = [float(imageOrientation[i]) for i in range(3)]
    Y = [float(imageOrientation[i+3]) for i in range(3)]
    Z = [float(imageOrientation[i+6]) for i in range(3)]
    di = float(pixelSpacing['cols'])
    dj = float(pixelSpacing['rows'])
    dk = float(pixelSpacing['frames'])
    S = [float(imagePosition[i]) for i in range(len(imagePosition))]
    M = np.mat([[X[0]*di,Y[0]*dj,Z[0]*dk,S[0] ],\
                [X[1]*di,Y[1]*dj,Z[1]*dk,S[1] ],\
                [X[2]*di,Y[2]*dj,Z[2]*dk,S[2] ],\
                [0,0,0,1 ]])
    [[x],[y],[z],[u]] = np.asarray(M * np.mat([[index_float[2]],[index_float[1]],[index_float[0]],[1.0]]))
    return (x,y,z)
    
    
def truncateVal( value):
    '''Truncates decimal value.
    
    :param value: Value to truncate
    
    '''
    if abs(value - math.floor(value)) < 0.5:
        value = math.floor(value)
    else:
        value = math.ceil(value)
    return int(value)
    
    

def getMaskForContour( contour , size2D , outline = 1 , fill = 1 , background = 0 ):
    '''Creates a  2D mask for a given contour:
    
    :param contour: should be a list of coordinates: either [ (x1,y1) , (x2,y2) ...] or [x1, y1, x2 , y2 ..] Note: xi, yj should be indices not real values
    :param size2D: 2D-tuple representing the image size
    :param outline:  int value to use to represent the contour's outline (1 by default)
    :param fill: int value to use to fill the contour's interior(1 by default)
    :param background: value to represent the background (0 by default)
        
    :returns: A 2D mask (numpy array) of size size2D is returned with "fill" value inside the polygon, "outline" on the edges, "background" in the background. 
    
    '''
    if not isinstance(contour,list) :
        raise ValueError("Wrong format for contour")
    if contour == []:
        print "Empty contour"
        return
    size2D = [size2D[1],size2D[0]]
    img = Image.new('L', size2D,background)
    if isinstance(contour[0],tuple):
        if len(contour) == 1:
            contour.append(contour[0])
    elif  isinstance(contour[0],(float,int)):
        contour.append(contour[0])
        contour.append(contour[1])
    else:
        strErr = "Contour: data is neither tuple nor float or int... type: %s"%type(contour[0])
        raise ValueError(strErr)
    try:
        ImageDraw.Draw(img).polygon(contour, outline = outline, fill = fill)
    except:
        print "Contour: %s "%contour
        raise ValueError("Contour Error")
    
    newArray = np.array(img, dtype = 'int8',order = "C")
    return newArray


    
    
def crossProduct( vectU , vectV):
    '''Computes the cross product U X V. U and V should be arrays, lists or tuples.
    
    :param vectU: first vector
    :param vectV: second vector
    
    [w1,w2,w3] = [u1,u2,u3] x [v1,v2,v3] = [u2v3 - u3v2 , u3v1 - u1v3, u1v2 - u2v1]
    
    :returns: The cross product
    
    '''
    vectW = ( vectU[1]*vectV[2] - vectU[2]*vectV[1] ,vectU[2]*vectV[0] - vectU[0]*vectV[2], vectU[0]*vectV[1] - vectU[1]*vectV[0] )
    return vectW
    


