/*
########################################################################
#
#  beamletSimulation_PythonExt.cpp 
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# September 2012
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
*/
  
#include "Python.h" // MUST be the first library included!!
#include "arrayobject.h" // from python numpy library 
#include <math.h>
#include <float.h>


/* ==== Prototypes =================================== */
// .... Python callable functions ..................

static PyObject *beamletSimul(PyObject *self, PyObject *args);


/* .... utility functions ..................*/
float radialEmmittance( float deff, float dmax, float sig0);
float interp(float xval, float * xList, float * yList, int lenList);
float centralAxisTerm(float dd, float ssd0, float deff, float z);
float offAxisTerm(float variance, float distPQ);
float offAxisTerm2Sigmas(float varX,float varY, float xrshift,float yrshift, float xBev, float yBev);
float dotProduct(float * p1,float * p2);
float distance( float * p1, float *p2);
void getRotatedPoint( float * newP , float * P, float * rotMat);
void findPointIndicesInGrid( int * ind, float * pointP1, float * spacing, float * x0y0z0,float * rotMat);