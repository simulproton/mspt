/*
########################################################################
#
#  fillCTArray_PythonExtDouble.cpp
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
##
########################################################################
*/
#include "fillCTArray_PythonExtDouble.h"

char errstr[200];  // error string that all routines have access to
/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"fillCTGrid", fillCTGrid, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _fillCTArray_PythonExtDouble in compile and linked
extern "C" { 
    void init_fillCTArray_PythonExtDouble()  {
        (void) Py_InitModule("_fillCTArray_PythonExtDouble", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * fillCTGrid(PyObject *self, PyObject *args){
/*
 Function to fill the CT grid based on a new grid shape
 Input:
    ctGrid : 3D numpy array representing the CT
    newShape : grid new shape
    indStart : index in the input CT grid where to start taking values.
 
 Output: 
    3D numpy array representing a new CT grid.
 
 */ 
 
    //Initialize pointers that will receive input data
 
    PyArrayObject * ctGrid = NULL; // We assume that ct Value have been previously correctly scaled.
    PyArrayObject * newShape = NULL; // shape [frames, rows,columns] of the new CT grid
    PyArrayObject * indStart = NULL; // indices in ctGrid of the first voxel of the new CTGrid
   
 
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!", 
      
        &PyArray_Type,&ctGrid, 
        &PyArray_Type,&newShape,
        &PyArray_Type,&indStart))  return NULL;

   
    if(!ctGrid) return NULL;
    if(!newShape) return NULL;
    if(!indStart) return NULL;
  
    
    //Cast values to double values: numpy creates array with double values
    if (ctGrid->descr->type_num != NPY_DOUBLE) ctGrid = (PyArrayObject*)PyArray_Cast(ctGrid,PyArray_DOUBLE);
    if (newShape->descr->type_num != NPY_INT) newShape = (PyArrayObject*)PyArray_Cast(newShape,PyArray_INT);
    if (indStart->descr->type_num != NPY_INT) indStart = (PyArrayObject*)PyArray_Cast(indStart,PyArray_INT);
    
    
    double  * dataNewCTGrid, * dataCTGrid;
    int * datanewShape, *dataIndStart;
    
    dataCTGrid = (double*) ctGrid -> data;  
    datanewShape = (int*) newShape -> data; 
    dataIndStart = (int*) indStart -> data; 

    int nbframesOrigin = ctGrid -> dimensions[0];
    int nbrowsOrigin = ctGrid -> dimensions[1];
    int nbcolsOrigin = ctGrid -> dimensions[2];
  
    int nbFramesCT = datanewShape[0];
    int nbRowsCT = datanewShape[1];
    int nbColsCT = datanewShape[2];
    int dimGrid[] ={nbFramesCT,nbRowsCT,nbColsCT} ;
    
    PyArrayObject * newCTGrid = NULL;
       // create output grid
    if(  ( newCTGrid = (PyArrayObject *) PyArray_FromDims(ctGrid->nd, dimGrid, PyArray_DOUBLE) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in fillCTGrid - newCTGrid allocation failed\n");
    }
    dataNewCTGrid = (double*) newCTGrid -> data;
    
    double minCT= DBL_MAX;
    for(int i= 0 ; i < nbFramesCT *nbRowsCT*nbColsCT ; i++ ){
        if(dataCTGrid[ i ] < minCT ) minCT= dataCTGrid[ i ];
    }    
    
    for(int i= 0 ; i < nbFramesCT *nbRowsCT*nbColsCT ; i++ ){
        dataNewCTGrid[i] = minCT;
    }
    
    
    
    double ctValue;
    int fStart,rStart,cStart;
    if( dataIndStart[0] < 0) fStart = absolute(dataIndStart[0]);
    else fStart = 0;
    if( dataIndStart[1] < 0) rStart = absolute(dataIndStart[1]);
    else rStart = 0;
    if( dataIndStart[2] < 0) cStart = absolute(dataIndStart[2]);
    else cStart = 0;
    
    int fEnd,rEnd,cEnd;
    if(  dataIndStart[0] + (nbFramesCT-1) >= nbframesOrigin ) fEnd = nbframesOrigin - dataIndStart[0];
    else fEnd = nbFramesCT;
    if(  dataIndStart[1] + (nbRowsCT-1) >= nbrowsOrigin ) rEnd = nbrowsOrigin - dataIndStart[1];
    else rEnd = nbRowsCT;    
     if(  dataIndStart[2] + (nbColsCT-1) >= nbcolsOrigin ) cEnd = nbcolsOrigin - dataIndStart[2];
    else cEnd = nbColsCT;    
   
    
    for(int f = fStart ; f < fEnd ; f++){ // End - iterate through frames
        for(int r = rStart ; r < rEnd ; r++){ // End - iterate through rows
            for(int c = cStart ; c < cEnd ; c++){ // End - iterate through columns
                int idx =  (f+dataIndStart[0]) * nbcolsOrigin * nbrowsOrigin + (r+dataIndStart[1]) * nbcolsOrigin + (c+dataIndStart[2]);
                ctValue = dataCTGrid[idx];
                dataNewCTGrid[ f * nbColsCT * nbRowsCT + r * nbColsCT + c ] = ctValue;
            } // End - iterate through columns
        } // End - iterate through rows
    } // End - iterate through frames

    return PyArray_Return(newCTGrid);
}

int absolute(int val){
    if (val < 0) return (-val);
    return val;

}

