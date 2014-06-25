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

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <float.h>
#include <assert.h>

#define SUCCESS 1
#define FAILURE 0


// 3D point
typedef struct
{
    float x;
    float y;
    float z;
} POINTXYZ;

// 3D grid
typedef struct
{
    POINTXYZ start;
    POINTXYZ inc;
    int nbCols;
    int nbRows;
    int nbFrames;
    float *matrix;
} VOLUME;


// Prototypes
int calcRadiolDepth(VOLUME *stpPwr_grid,VOLUME *radiological_depth_grid,
             POINTXYZ point1,POINTXYZ point2, int * tabCount);
float coordPlane(int idx, float coordPlane1, float spacing);
void sortAlphaMinMax(float * alphasMinMax);
int computeAlpha(int idx, float coord1, float coord2,float coordPlane1, float spacing , float * ptAlpha);
void findAlphaMinMax(float * alphasX, float * alphasY, float * alphasZ, float * alphaMinMax);
void findIndMinMax(float coordPlane1,float coordPlaneN,float spacing,float coord1, float coord2, int nbPlanes, float * alphaMinMax , int * indMinMax);
int calculateParametricSets(int * indMinMax , float coord1, float coord2, float coordPlane1, float spacing, float * alphaSet, float *alphaMinMax);
void mergeSets(float * alphaMinMax, float * alphaX, int lenAlphaX, float * alphaY, int lenAlphaY , float * alphaZ, int lenAlphaZ, float * alpha, float d12, float spacing);
void merge2Arrays(float * dest, float * array1, int lenArray1, float * array2, int lenArray2);
void computeVoxelsLengths(float * alpha, int lenAlpha, float * lengthsVoxels, float d12);
float euclideanDistance( POINTXYZ p1, POINTXYZ p2);
void calculateVoxelIndices(POINTXYZ p1, POINTXYZ p2, POINTXYZ plane1, POINTXYZ spacing, float * alpha, int lenAlpha, int ** indList, int ** indMinMax);
int roundFloat( float value);
float calculateRadiologicalPath(VOLUME * densGrid, VOLUME * radGrid, int ** listInd, float * listAlpha,int  lenAlpha,float d12, int * countPtr);
float absolute(float val);