/*
########################################################################
#
#  fillDoseArray_PythonExtDouble.cpp
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# September 2012
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
##
########################################################################
*/

#include "fillDoseArray_PythonExtDouble.h"

const int interpol = 0; //0: not interpolate, 1 interpolate initial dose grid

char errstr[200];  // error string that all routines have access to
/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"fillDoseGrid", fillDoseGrid, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _fillDoseArray_PythonExtDouble in compile and linked
extern "C" { 
    void init_fillDoseArray_PythonExtDouble()  {
        (void) Py_InitModule("_fillDoseArray_PythonExtDouble", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * fillDoseGrid(PyObject *self, PyObject *args){
/*fillDoseGrid aims to fill a 3D array with the planned dose at each voxel. The originalDoseGrid and newDoseGrid don't always have the same
    size. That is why we need to associate to each coordinate in newDoseGrid with the coordinates in originalDoseGrid and set the voxel value value in newDoseGrid with the voxel value in originalDoseGrid
        
        - dimNewGrid:Dimensions of the 3D array to fill with the dose contained in originalDoseGrid. Its size is m x n x q
        - originalDoseGrid: 3D array containing the planned dose extracted from RD dicom file. Its size can be different than m x n x q, it can be smaller or greater.  
        - imagePosition: coordinates (in mm) of the top left corner on the first slice of originalDoseGrid
        - imageEndPosition: coordinates (in mm) of the bottom right corner on the last slice of originalDoseGrid
        - imageOrientation: axis information in orginialDoseGrid : 1D array with 9 elements
        - spacingInfoCT : spacing information (mm) in newDoseGrid: x:spacingInfoCT[0] ,y:spacingInfoCT[1] ,z:spacingInfoCT[2] ,
        - spacingInfoRD : spacing information (mm) in originalDoseGrid (from RD File) : x:spacingInfoRD[0] ,y:spacingInfoRD[1] ,z:spacingInfoRD[2] 
        - imagePositionRD : x,y,z coordinates of the first top left voxel of the dose data.
        The goal of the function is to find in originalDoseGrid the position of each spot in newDoseGrid in order to set the voxel of newDoseGrid to the corresponding dose value.
*/

    //Initialize pointers that will receive input data
     
    PyArrayObject * dimNewGrid = NULL;
    PyArrayObject * originalDoseGrid = NULL; 
    PyArrayObject * imagePosition = NULL; 
    PyArrayObject * imageEndPosition = NULL;
    PyArrayObject * imageOrientation = NULL;
    PyArrayObject * spacingInfoCT = NULL;
    PyArrayObject * spacingInfoRD = NULL;
    PyArrayObject * imagePositionRD = NULL; 
    
    
    
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!", 
        &PyArray_Type,&dimNewGrid, 
        &PyArray_Type,&originalDoseGrid, 
        &PyArray_Type,&imagePosition, 
        &PyArray_Type,&imageEndPosition,
        &PyArray_Type,&imageOrientation,
        &PyArray_Type,&spacingInfoCT,
        &PyArray_Type,&spacingInfoRD,
        &PyArray_Type,&imagePositionRD))  return NULL;

    if(!dimNewGrid) return NULL;
    if(!originalDoseGrid) return NULL;
    if(!imagePosition) return NULL;
    if(!imageEndPosition) return NULL;
    if(!imageOrientation) return NULL;
    if(!spacingInfoCT) return NULL;
    if(!spacingInfoRD) return NULL;
    if(!imagePositionRD) return NULL;
    
    //Cast values to double values: numpy creates array with double values
    if (dimNewGrid->descr->type_num != PyArray_INT) dimNewGrid = (PyArrayObject*)PyArray_Cast(dimNewGrid,PyArray_INT);
    if (originalDoseGrid->descr->type_num != NPY_DOUBLE) originalDoseGrid = (PyArrayObject*)PyArray_Cast(originalDoseGrid,PyArray_DOUBLE);
    if (imagePosition->descr->type_num != NPY_DOUBLE) imagePosition = (PyArrayObject*)PyArray_Cast(imagePosition,PyArray_DOUBLE);
    if (imageEndPosition->descr->type_num != NPY_DOUBLE) imageEndPosition = (PyArrayObject*)PyArray_Cast(imageEndPosition,PyArray_DOUBLE);
    if (imageOrientation->descr->type_num != NPY_DOUBLE) imageOrientation = (PyArrayObject*)PyArray_Cast(imageOrientation,PyArray_DOUBLE);
    if (spacingInfoCT->descr->type_num != NPY_DOUBLE) spacingInfoCT = (PyArrayObject*)PyArray_Cast(spacingInfoCT,PyArray_DOUBLE);
    if (spacingInfoRD->descr->type_num != NPY_DOUBLE) spacingInfoRD = (PyArrayObject*)PyArray_Cast(spacingInfoRD,PyArray_DOUBLE);
    if (imagePositionRD->descr->type_num != NPY_DOUBLE) imagePositionRD = (PyArrayObject*)PyArray_Cast(imagePositionRD,PyArray_DOUBLE);
    
    double  * dataOriginalDoseGrid, * dataImagePosition, * dataImageEndPosition,* dataImagePositionRD;
    double * dataImageOrientation, * dataSpacingInfoCT, * dataSpacingInfoRD;
    int * dataDimNewGrid;
    
    dataDimNewGrid = (int*) dimNewGrid -> data;
    dataOriginalDoseGrid = (double*) originalDoseGrid -> data;  
    dataImagePosition = (double*) imagePosition -> data; 
    dataImageEndPosition = (double*) imageEndPosition -> data;
    dataImageOrientation = (double*) imageOrientation -> data;  
    dataSpacingInfoCT = (double*) spacingInfoCT -> data;  
    dataSpacingInfoRD = (double*) spacingInfoRD -> data;
    dataImagePositionRD = (double*) imagePositionRD -> data; 
    
    PyArrayObject * newDoseGrid = NULL;
    if(  ( newDoseGrid = (PyArrayObject *) PyArray_FromDims(originalDoseGrid ->nd, dataDimNewGrid, PyArray_DOUBLE) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in fillDoseGrid - doseGride allocation failed\n");
    }
    double * dataNewDoseGrid = (double*) newDoseGrid -> data;
    
    
    
    int nbFramesCT = dataDimNewGrid[0];
    int nbRowsCT = dataDimNewGrid[1];
    int nbColsCT = dataDimNewGrid[2];
    
    int nbFramesRD = originalDoseGrid->dimensions[0];
    int nbRowsRD = originalDoseGrid->dimensions[1];
    int nbColsRD = originalDoseGrid-> dimensions[2];
    
    int indicesRD[3]; // voxel indices in RD pixel array
    double point[3];
    //Neighbors spots in the RD image: 
    // If we consider a cube centered on the current RD spot we want to interpolate
    // We consider the image orientation, and indices expressed as (frame, row, col):
    // Neighbor1 is on the front face, top left: (-1,-1,-1)
    // Neighbor2 is on the back face, top left: (+1,-1,-1)
    // Neighbor3 is on the back face, top right:(+1,-1,+1)
    // Neighbor4 is on the back face, top right:(-1,-1,+1)
    // Neighbor5 is on the front face, top left: (-1,+1,-1)
    // Neighbor6 is on the back face, top left: (+1,+1,-1)
    // Neighbor7 is on the back face, top right:(+1,+1,+1)
    // Neighbor8 is on the back face, top right:(-1,+1,+1)
    double pointNeigh[8][3];
    int tmpInd[3];
    double valNeigh[8];
    double tmpPoint[3];


    
    
    
    double transformMatrix[16]= { 0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0,
                                 0,0,0,0};
                                 
                                 
    int performInterpol = 0;
    if ( interpol == 1){
        for ( int spIdx = 0; spIdx < 3; spIdx++)
            if(dataSpacingInfoCT[spIdx] != dataSpacingInfoRD[spIdx]) performInterpol = 1;
    }                                 
    createTransformationMatrix(dataImagePositionRD, dataImageOrientation, dataSpacingInfoRD, transformMatrix);
//    createTransformationMatrixIndToCoord(dataImagePosition, dataImageOrientation, dataSpacingInfoRD, transformMatrixIndToCoord);
    for(int f = 0 ; f < nbFramesCT ; f++){ // End - iterate through frames
        for(int r = 0 ; r < nbRowsCT ; r++){ // End - iterate through rows
            for(int c = 0 ; c < nbColsCT ; c++){ // End - iterate through columns
                point[0] = dataImagePosition[0] + c * dataSpacingInfoCT[2];// col
                point[1] = dataImagePosition[1] + r * dataSpacingInfoCT[1];// row
                point[2] = dataImagePosition[2] + f * dataSpacingInfoCT[0]; // frame
                //Find "point" indices in RD dose grid
                coordinatesToIndexFromImagePosition( point,transformMatrix, indicesRD);
                if ( interpol == 1 && performInterpol == 1){
                
                     ///////// Neighbor1 is on the front face, top left: (-1,-1,-1) 
                    tmpPoint[0] = point[0] - dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] - dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] - dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                     valNeigh[0] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];
                    }
                    else{
                    valNeigh[0] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++){
                        pointNeigh[0][i]=tmpPoint[i];
                    }
                      
                    ////// Neighbor2 is on the back face, top left: (+1,-1,-1)
                    tmpPoint[0] = point[0] - dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] - dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] + dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                     valNeigh[1] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                     }
                    else{
                    valNeigh[1] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[1][i]=tmpPoint[i];
                        
                        
                    ///// Neighbor3 is on the back face, top right:(+1,-1,+1) 
                    tmpPoint[0] = point[0] + dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] - dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] + dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD) valNeigh[2] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];
                    else{
                     valNeigh[2] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];

                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[2][i]=tmpPoint[i];
                       
                    //// Neighbor4 is on the back face, top right:(-1,-1,+1)
                    tmpPoint[0] = point[0] + dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] - dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] - dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){ 
                     valNeigh[3] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                    }
                    else{
                    valNeigh[3] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[3][i]=tmpPoint[i];
                        
                    //// Neighbor5 is on the front face, top left: (-1,+1,-1) 
                    tmpPoint[0] = point[0] - dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] + dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] - dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                        valNeigh[4] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                    }
                    else{
                    valNeigh[4] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[4][i]=tmpPoint[i];
                        
                    //// Neighbor6 is on the back face, top left: (+1,+1,-1)
                    tmpPoint[0] = point[0] - dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] + dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] + dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                      valNeigh[5] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                     }
                    else{
                    valNeigh[5] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[5][i]=tmpPoint[i];                   
 
 
                    //// Neighbor7 is on the back face, top right:(+1,+1,+1)
                    tmpPoint[0] = point[0] + dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] + dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] + dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                      valNeigh[6] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                     }
                    else{
                    valNeigh[6] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[6][i]=tmpPoint[i];
                 
                        
                    //// Neighbor8 is on the back face, top right:(-1,+1,+1) 
                    tmpPoint[0] = point[0] + dataSpacingInfoCT[2];// col
                    tmpPoint[1] = point[1] + dataSpacingInfoCT[1];// row
                    tmpPoint[2] = point[2] - dataSpacingInfoCT[0];// frame
                    coordinatesToIndexFromImagePosition( tmpPoint,transformMatrix, tmpInd);
                    if(tmpInd[0] < 0 || tmpInd[0] >= nbFramesRD
                     || tmpInd[1] < 0 || tmpInd[1] >= nbRowsRD
                     || tmpInd[2] < 0 || tmpInd[2] >= nbColsRD){
                    valNeigh[7] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];

                     }
                    else{
                    valNeigh[7] = dataOriginalDoseGrid[ tmpInd[0]* nbColsRD * nbRowsRD + tmpInd[1] * nbColsRD + tmpInd[2] ];
                    }
                    for( int i = 0; i < 3; i++)
                        pointNeigh[7][i]=tmpPoint[i];
                                         
                    double doseValue = trilinearInterpolation(valNeigh ,  pointNeigh , point);
                    if (doseValue < 0) exit(1);
                    dataNewDoseGrid[ f * nbColsCT * nbRowsCT + r * nbColsCT + c ] = doseValue;                       
                    
                }
                else{
                    if (interpol == 0){
                    //Set the dose value in the new dose grid
                    dataNewDoseGrid[ f * nbColsCT * nbRowsCT + r * nbColsCT + c ] = dataOriginalDoseGrid[ indicesRD[0]* nbColsRD * nbRowsRD + indicesRD[1] * nbColsRD + indicesRD[2] ];
                    
                    }
                    else{
                        PyErr_SetString(PyExc_ValueError,"Wrong interpol value in fillDoseArray\n");
                    }
                
                }

        
            } // End - iterate through columns
        } // End - iterate through rows
    } // End - iterate through frames
    

    return PyArray_Return(newDoseGrid);
}

double trilinearInterpolation(double * valNeigh ,  double  pointNeigh[8][3] , double * pointRD){

    
    // p1: between neighbors 1 and 4
    // p2: between neighbors 2 and 3
    // p3: between neighbors 5 and 8
    // p4: between neighbors 6 and 7
    double valP1,valP2,valP3,valP4;
    valP1 = computeBarycenter(valNeigh[0], valNeigh[3], pointRD[0], pointNeigh[0][0], pointNeigh[3][0]);
    valP2 = computeBarycenter(valNeigh[1], valNeigh[2], pointRD[0], pointNeigh[1][0], pointNeigh[2][0]);
    valP3 = computeBarycenter(valNeigh[4], valNeigh[7], pointRD[0], pointNeigh[4][0], pointNeigh[7][0]);
    valP4 = computeBarycenter(valNeigh[5], valNeigh[6], pointRD[0], pointNeigh[5][0], pointNeigh[6][0]);

    // q1: between p1 and p2
    // q2: between p4 and p4
    double valQ1, valQ2;
    valQ1 = computeBarycenter(valP1, valP2, pointRD[2], pointNeigh[0][2], pointNeigh[1][2]);
    valQ2 = computeBarycenter(valP3, valP4, pointRD[2], pointNeigh[4][2], pointNeigh[5][2]);
    
    double ret = computeBarycenter(valQ1, valQ2, pointRD[1], pointNeigh[0][1], pointNeigh[4][1]);

    return ret;

}

double computeBarycenter(double valN1, double valN2, double coordP, double coordN1, double coordN2){
    
    if (valN1 == -1){
        if (valN2 == -1){
            return -1.0;
        }
        else return valN2;
    }
    else{
        if (valN2 == -1){
            return valN1;
        }
        if (coordN1 == coordN2) return valN1;
        
        double alpha1 = absolute(coordP - coordN1)/ absolute(coordN2 - coordN1);
        
        double alpha2 = absolute(coordP - coordN2)/ absolute(coordN2 - coordN1);
        return (alpha1 * valN1 + alpha2 * valN2);
    }
}


void createTransformationMatrix(double * imagePosition, double * imageOrientation, double * pixelSpacing,double * transformMatrix){
        /* Returns the transformation matrix used to go from one coordinates in mm to indices in a Dicom Array
        Image Position: refers to coordinates (x0,y0,z0) of the top left corner in the image 
        Image Orientation: unit vectors in the 3 directions with respect to the patient coordinate system
        Pixel Spacing: pixel spacing in the 3 directions
        transformMatrix: tha matrix to fill
     
        Formula: look at page 410 from Dicom Documentation */
        
        double * X = (double *)malloc(sizeof(double) * 3);
        double * Y = (double *)malloc(sizeof(double) * 3);
        double * Z = (double *)malloc(sizeof(double) * 3);
        double di = pixelSpacing[2] ; // cols spacing
        double dj = pixelSpacing[1] ; // rows spacing
        double dk = pixelSpacing[0] ; // frames spacing
        double * S = imagePosition;
        
        X[0] = imageOrientation[0];
        X[1] = imageOrientation[1];
        X[2] = imageOrientation[2];
        Y[0] = imageOrientation[3];
        Y[1] = imageOrientation[4];
        Y[2] = imageOrientation[5];
        Z[0] = imageOrientation[6];
        Z[1] = imageOrientation[7];
        Z[2] = imageOrientation[8];
    
        double  initialMatrix[16] = { X[0]*di , Y[0]*dj , Z[0]*dk , S[0],
                                       X[1]*di , Y[1]*dj , Z[1]*dk , S[1],
                                       X[2]*di , Y[2]*dj , Z[2]*dk , S[2],
                                       0       ,    0    ,     0   ,  1};
                                    
        
        if (! invertMatrix(initialMatrix, transformMatrix) ){
            PyErr_SetString(PyExc_ValueError,"Error - matrix inversion - det = 0");

        }    
        free(X);
        free(Y);
        free(Z);

    
}



void createTransformationMatrixIndToCoord(double * imagePosition, double * imageOrientation, double * pixelSpacing,double * transformMatrix){
        /* Returns the transformation matrix used to go from indices to coordinates in mm in a Dicom Array
        Image Position: refers to coordinates (x0,y0,z0) of the top left corner in the image 
        Image Orientation: unit vectors in the 3 directions with respect to the patient coordinate system
        Pixel Spacing: pixel spacing in the 3 directions
        transformMatrix: tha matrix to fill
     
        Formula: look at page 410 from Dicom Documentation */
        
        double * X = (double *)malloc(sizeof(double) * 3);
        double * Y = (double *)malloc(sizeof(double) * 3);
        double * Z = (double *)malloc(sizeof(double) * 3);
        double di = pixelSpacing[2] ; // cols spacing
        double dj = pixelSpacing[1] ; // rows spacing
        double dk = pixelSpacing[0] ; // frames spacing
        double * S = imagePosition;
        
        X[0] = imageOrientation[0];
        X[1] = imageOrientation[1];
        X[2] = imageOrientation[2];
        Y[0] = imageOrientation[3];
        Y[1] = imageOrientation[4];
        Y[2] = imageOrientation[5];
        Z[0] = imageOrientation[6];
        Z[1] = imageOrientation[7];
        Z[2] = imageOrientation[8];
    
        double  initialMatrix[16] = { X[0]*di , Y[0]*dj , Z[0]*dk , S[0],
                                       X[1]*di , Y[1]*dj , Z[1]*dk , S[1],
                                       X[2]*di , Y[2]*dj , Z[2]*dk , S[2],
                                       0       ,    0    ,     0   ,  1};
        for(int i = 0 ; i < 16; i ++){
            transformMatrix[i] = initialMatrix[i];
        }

        free(X);
        free(Y);
        free(Z);

    
}

void  indexToCoordFromImagePosition( double * point, double * matrix,int *indices){
        /* Returns (frame, row, column) indices (pixel location))for a given point (x,y,z) provided in an image 
        coordinate system (mm).
        invMatrix : the transformation matrix
     
        Formula: look at page 410 from Dicom Documentation */
        

        
        int tmp[3];
        //Column
        tmp[0] = round(matrix[0] * indices[0] +  matrix[1] * indices[1] + matrix[2] * indices[2] + matrix[3]);
        //Row
        tmp[1] = round(matrix[4] * indices[0] +  matrix[5] * indices[1] + matrix[6] * indices[2] + matrix[7]);
        //Frame
        tmp[2] = round(matrix[8] * indices[0] +  matrix[9] * indices[1] + matrix[10] * indices[2] + matrix[11]);
      
        point[0] = tmp[2];
        point[1] = tmp[1];
        point[2] = tmp[0];
}



void  coordinatesToIndexFromImagePosition( double * point, double * invMatrix,int *indices){
        /* Returns (frame, row, column) indices (pixel location))for a given point (x,y,z) provided in an image 
        coordinate system (mm).
        invMatrix : the transformation matrix
     
        Formula: look at page 410 from Dicom Documentation */
        

        
        int tmp[3];
        //Column
        tmp[0] = round(invMatrix[0] * point[0] +  invMatrix[1] * point[1] + invMatrix[2] * point[2] + invMatrix[3]);
        //Row
        tmp[1] = round(invMatrix[4] * point[0] +  invMatrix[5] * point[1] + invMatrix[6] * point[2] + invMatrix[7]);
        //Frame
        tmp[2] = round(invMatrix[8] * point[0] +  invMatrix[9] * point[1] + invMatrix[10] * point[2] + invMatrix[11]);
      
        indices[0] = tmp[2];
        indices[1] = tmp[1];
        indices[2] = tmp[0];
}

bool invertMatrix(double m[16], double invOut[16])
{   // Functions that invert a 4x4 matrix

    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}


double absolute(double val){
    if (val < 0) return (-val);
    return val;

}
