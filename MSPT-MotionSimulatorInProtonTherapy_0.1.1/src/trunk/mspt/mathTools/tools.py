########################################################################
#
# File: tools.py
#  
# Proton Therapy Simulator Project
# Written by by Paul Morel, 
#               LIGM, Universite Paris-Est Marne La Vallee (France)
#               paul.morel@univ-mlv.fr
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
#
########################################################################


########################################################################
#
# Description:
#  
# Tools for mathematical operations.
# 
#
########################################################################


import numpy as np
import scipy as sp


    
def interp(xValue, x1,x2,y1,y2):
    '''Interpolates linearly between (x1,y1) and (x2,y2) and return the value for xValue
    
    :param xValue: value to interpolate
    :param x1: x coordinate of the first point
    :param x2: x coordinate of the second point
    :param y1: y coordinate of the first point
    :param y2: y coordinate of the second point
    
    :returns: Interpolated value.
    
    '''
    slope = (y2-y1) / (x2-x1)
    intercept = y1 - slope * x1
    return (xValue * slope + intercept)
    
def rotateCoordinates(rotMatrix,point, typeFloat):
    '''Computes the rotation of a point.
    
    :param rotMatrix: 3x3 rotation matrix (numpy array)
    :param point: coordinates (x,y,z) of the point to rotate
    :param typeFloat: type of numpy float to use; 'float32' or 'float64'
    
    :returns: The coordinates of the rotated point (x',y',z')
    
    '''
      
      
    row0 = rotMatrix[0] 
    row1 = rotMatrix[1]
    row2 = rotMatrix[2]                    

    rotatedPoint = np.array([row0[0]*point[0] + row0[1]*point[1] + row0[2]*point[2] , \
                            row1[0]*point[0] + row1[1]*point[1] + row1[2]*point[2] ,\
                            row2[0]*point[0] + row2[1]*point[1] + row2[2]*point[2] ],dtype= typeFloat,order='C')

    return rotatedPoint    
    
def distEuclideanP1P2( p1,p2):
    '''Computes the Euclidean distance
    
    :param p1: Coordinates point1 (x,y,z)
    :param p2: Coordinates point2 (x,y,z)

    :returns: The euclidean distance
    
    '''
   
    return np.sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) +\
                            (p1[1]-p2[1])*(p1[1]-p2[1]) +\
                            (p1[2]-p2[2])*(p1[2]-p2[2]))
                            
def rot3D(X,Y,Z,n,nr,floatType):
    '''Computes the rotation of the X,Y,Z coordinates.
    
    :param X: 3D numpy array with x coordinates of each voxel
    :param Y: 3D numpy array with y coordinates of each voxel
    :param Z: 3D numpy array with z coordinates of each voxel
    :param n: unit vectors (1 by column) in the initial coordinate system
    :param nr: unit vectors (1 by column) in the new coordinate system
    :param floatType: Type of numpy float : 'float32' or 'float64'
    
    :returns: (Xr,Yr,Zr) a tuple of 3D arrays corresponding to the rotated coordinates
    
    '''
    #extract the unit vectors for each coordinate system
    i = n[:,0];
    j = n[:,1];
    k = n[:,2];

    ir = nr[:,0];
    jr = nr[:,1];
    kr = nr[:,2];
    
    # calculate the transformation matrix
    T = np.array([ [np.dot(np.transpose(i),ir),  np.dot(np.transpose(j),ir) , np.dot(np.transpose(k),ir)], \
                    [np.dot(np.transpose(i),jr),  np.dot(np.transpose(j),jr),  np.dot(np.transpose(k),jr)], \
                    [np.dot(np.transpose(i),kr),  np.dot(np.transpose(j),kr),  np.dot(np.transpose(k),kr)] ], dtype=floatType)
    #Warning T created is nr^-1 !!! This is due to the rotation matrix that rotates counter clockwise, i.e opposite direction than gantry


    
    #perform the transformation
    Xr = np.dot(T[0,0],X) + np.dot(T[0,1],Y) + np.dot(T[0,2],Z)
    Yr = np.dot(T[1,0],X) + np.dot(T[1,1],Y) + np.dot(T[1,2],Z)
    Zr = np.dot(T[2,0],X) + np.dot(T[2,1],Y) + np.dot(T[2,2],Z)

    return (Xr,Yr,Zr)
    
def fromCTToIECFixedCoord(X,Y,Z,isocenter):
    '''Converts coordinates from Dicom Patient coordinate system to IEC fixed coordinate system:
    
    :param X: matrix of x coordinates for voxels
    :param Y: matrix of y coordinates for voxels
    :param Z: matrix of z coordinates for voxels
    :param isocenter: [x,y,z] coordinates of the isocenter defined in the RP dicom (plan) file 
    
    :returns: A tuple:
    
        #. X_iec_f = X - x_isocenter 
        #. Y_iec_f = Y - y_isocenter
        #. Z_iec_f = Z - z_isocenter
        
    '''
    #Translation to express the origin of the dicom patient coordinate system in the IEC fixed coordinate system
    X_translated = X - isocenter[0] 
    Y_translated = Y - isocenter[1]
    Z_translated = Z - isocenter[2]
    
    #Take into account the fact that in IEC fixed coord system, Z is positive when Y is negative in the dicom patient coordinate system. X doesn't change, and Y in IEC fixed is the same as Z in dicom patient coord systm.
    X_iec_f = X_translated 
    Y_iec_f = Z_translated
    Z_iec_f = -1*Y_translated
    
    return (X_iec_f,Y_iec_f,Z_iec_f)
    
def fromIECFixedToCTCoord(XiecF,YiecF,ZiecF,isocenter):
    '''Converts coordinates from IEC fixed coordinate system to Dicom Patient coordinate system.
    
    :param XiecF: matrix of x coordinates for voxels
    :param YiecF: matrix of y coordinates for voxels
    :param ZiecF: matrix of z coordinates for voxels
    :param isocenter: [x,y,z] coordinates of the isocenter defined in the RP dicom (plan) file 
    
    :returns: A tuple:
    
        #. X_iec_f = X - x_isocenter 
        #. Y_iec_f = Y - y_isocenter
        #. Z_iec_f = Z - z_isocenter
        
    '''
    X = XiecF + isocenter[0]
    Y = -1*ZiecF + isocenter[1]
    Z = YiecF + isocenter[2]
    
    return (X,Y,Z)



def rotationMatrixForAngle( angle): # angle in degrees
                                         #Rotation Around y axis clockwise
    '''Computes and returns a rotation matrix for a given angle in degrees.The rotation is arounf the y axis of the 
    IEC fixed coordinate system.
    
    :param angle: Angle in degrees
    
    :returns: A 3x3 array.
    
    '''
    C = sp.pi/180.0
    rotMatrix = np.array(([np.cos(angle * C) , 0 , np.sin(angle * C)],[0 , 1 , 0],[-np.sin(angle * C) , 0 , np.cos(angle * C)]) )
    assert(rotMatrix.shape == (3,3) )
    return rotMatrix