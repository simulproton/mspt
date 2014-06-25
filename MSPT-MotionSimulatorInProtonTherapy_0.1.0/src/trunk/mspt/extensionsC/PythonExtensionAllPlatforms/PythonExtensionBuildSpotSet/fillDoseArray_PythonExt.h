/*
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
*/

#include "Python.h" // MUST be the first library included!!
#include "arrayobject.h" // from python numpy library 
#include <math.h>


/* ==== Prototypes =================================== */

// .... Python callable functions ..................

static PyObject *fillDoseGrid(PyObject *self, PyObject *args);


/* .... utility functions ..................*/
float trilinearInterpolation(float * valNeigh ,  float  pointNeigh[8][3] , float * pointRD);
float computeBarycenter(float valN1, float valN2, float coordP, float coordN1, float coordN2);
void createTransformationMatrix(float * imagePosition, float * imageOrientation, float * pixelSpacing,float * transformMatrix);
void createTransformationMatrixIndToCoord(float * imagePosition, float * imageOrientation, float * pixelSpacing,float * transformMatrix);
void indexToCoordFromImagePosition( float * point, float * matrix,int *indices);
void  coordinatesToIndexFromImagePosition( float * point, float * invMatrix,int *indices);
bool invertMatrix(float m[16], float invOut[16]);
float absolute(float val);

