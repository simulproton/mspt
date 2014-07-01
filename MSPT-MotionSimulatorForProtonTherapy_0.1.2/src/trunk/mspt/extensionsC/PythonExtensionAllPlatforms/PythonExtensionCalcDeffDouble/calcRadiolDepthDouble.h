/*
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
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#define SUCCESS 1
#define FAILURE 0


// 3D point
typedef struct
{
    double x;
    double y;
    double z;
} POINTXYZ;

// 3D grid
typedef struct
{
    POINTXYZ start;
    POINTXYZ inc;
    int nbCols;
    int nbRows;
    int nbFrames;
    double *matrix;
} VOLUME;


// Prototypes
int calcRadiolDepthDouble(VOLUME *stpPwr_grid,VOLUME *radiological_depth_grid,
             POINTXYZ point1,POINTXYZ point2, int * tabCount);
double coordPlane(int idx, double coordPlane1, double spacing);
void sortAlphaMinMax(double * alphasMinMax);
int computeAlpha(int idx, double coord1, double coord2,double coordPlane1, double spacing , double * ptAlpha);
void findAlphaMinMax(double * alphasX, double * alphasY, double * alphasZ, double * alphaMinMax);
void findIndMinMax(double coordPlane1,double coordPlaneN,double spacing,double coord1, double coord2, int nbPlanes, double * alphaMinMax , int * indMinMax);
int calculateParametricSets(int * indMinMax , double coord1, double coord2, double coordPlane1, double spacing, double * alphaSet, double *alphaMinMax);
void mergeSets(double * alphaMinMax, double * alphaX, int lenAlphaX, double * alphaY, int lenAlphaY , double * alphaZ, int lenAlphaZ, double * alpha, double d12, double spacing);
void merge2Arrays(double * dest, double * array1, int lenArray1, double * array2, int lenArray2);
void computeVoxelsLengths(double * alpha, int lenAlpha, double * lengthsVoxels, double d12);
double euclideanDistanceDouble( POINTXYZ p1, POINTXYZ p2);
void calculateVoxelIndices(POINTXYZ p1, POINTXYZ p2, POINTXYZ plane1, POINTXYZ spacing, double * alpha, int lenAlpha, int ** indList, int ** indMinMax);
int roundFloat( double value);
double calculateRadiologicalPath(VOLUME * densGrid, VOLUME * radGrid, int ** listInd, double * listAlpha,int  lenAlpha,double d12,int * countPtr);
double absolute(double val);