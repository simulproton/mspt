/*
########################################################################
#
#  beamletSimulation_PythonExtDouble.h
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# September 2012
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
*/
#include "Python.h" // MUST be the first library included!!
#include "arrayobject.h" // from python numpy library 
#include <math.h>
#include <float.h>


/* ==== Prototypes =================================== */
// .... Python callable functions ..................

static PyObject *beamletSimulDouble(PyObject *self, PyObject *args);


/* .... utility functions ..................*/
double radialEmmittance( double deff, double dmax, double sig0);
double interp(double xval, double * xList, double * yList, int lenList);
double centralAxisTerm(double dd, double ssd0, double deff, double z);
double offAxisTerm(double variance, double distPQ);
double offAxisTerm2Sigmas(double varX,double varY, double xrshift,double yrshift, double xBev, double yBev);
double dotProduct(double * p1,double * p2);
double distance( double * p1, double *p2);
void getRotatedPoint( double * newP , double * P, double * rotMat);
void findPointIndicesInGrid( int * ind, double * pointP1, double * spacing, double * x0y0z0,double * rotMat);