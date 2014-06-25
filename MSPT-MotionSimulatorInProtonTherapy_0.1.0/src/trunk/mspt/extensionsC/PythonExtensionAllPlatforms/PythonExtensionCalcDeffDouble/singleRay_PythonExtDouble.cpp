/*
########################################################################
#
#  singleRay_PythonExtDouble.cpp 
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# 21 Feb. 2014
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
#   
########################################################################

Python extension coded in C:

The extension receives computes and return the voxels encountered by a ray in a ray.

*/

#include "singleRay_PythonExtDouble.h"

char errstr[200];  // error string that all routines have access to
/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"singleRayTrace", singleRayTrace, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _singleRay_PythonExtDouble in compile and linked
extern "C" { 
    void init_singleRay_PythonExtDouble()  {
        (void) Py_InitModule("_singleRay_PythonExtDouble", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * singleRayTrace(PyObject *self, PyObject *args){
    /*
    This function is the python extension. It calculates the voxels of a volume encountered by a ray:
    Inputs:
        grid : 3D numpy array representing the volume 
        startVec : numpy array (1,3) : Coordinates of the first voxel (0,0,0)
        incVec : numpy array (1,3): Increment vector (spacing) when incrementing (+1,+1,+1) the voxels indices.
        sourceVec : numpy array (1,3): Coordinates of the source of the ray
        targetVec : numpy array (1,3): Coordinates of the target of the ray
    
    Output:
        Numpy array (nbVoxels,4) in which is stored a list of voxel: for each voxel (frame,row,column) indices are stored along 
        with the length of the ray going through the voxel.
    */

    //Initialize pointers that will receive input data
    PyArrayObject * grid = NULL;
    PyArrayObject * startVec = NULL;
    PyArrayObject * incVec = NULL;
    PyArrayObject * sourceVec = NULL;
    PyArrayObject * targetVec = NULL;
    
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!", 
        &PyArray_Type,&grid, 
        &PyArray_Type,&startVec, 
        &PyArray_Type,&incVec, 
        &PyArray_Type,&sourceVec,
        &PyArray_Type,&targetVec))  return NULL;

    if(!grid) return NULL;
    if(!startVec) return NULL;
    if(!incVec) return NULL;
    if(!sourceVec) return NULL;
    if(!targetVec) return NULL;
    
    //Cast values to double values: numpy creates array with double values
    if (grid->descr->type_num != NPY_DOUBLE) grid = (PyArrayObject*)PyArray_Cast(grid,PyArray_DOUBLE);
    if (startVec->descr->type_num != NPY_DOUBLE) startVec = (PyArrayObject*)PyArray_Cast(startVec,PyArray_DOUBLE);
    if (incVec->descr->type_num != NPY_DOUBLE) incVec = (PyArrayObject*)PyArray_Cast(incVec,PyArray_DOUBLE);
    if (sourceVec->descr->type_num != NPY_DOUBLE) sourceVec = (PyArrayObject*)PyArray_Cast(sourceVec,PyArray_DOUBLE);
    if (targetVec->descr->type_num != NPY_DOUBLE) targetVec = (PyArrayObject*)PyArray_Cast(targetVec,PyArray_DOUBLE);


    // scalars
    int i,j,k,M,N,Q;

    // pointers to vectors that do not need to be freed
    int *gridDim, *startDim, *incDim ,*sourceVecDim, *targetVecDim;
    double *gridPtr,*startPtr,*incPtr, *sourceVecPtr, *targetVecPtr;
    // int *indicesPtr;
    double *indicesPtr;
    int dimList[2];

    // pointers to vectors that need to be freed
    double *x,*y,*z;

    // other
    int numDimGrid;
    
    // point structures that will be used for raytrace operation
    POINTXYZ p1;
    POINTXYZ p2;

    // grid structures for raytracing operation
    VOLUME gridValues;
    
    INDICES indices;

    
    //Output data 
    PyArrayObject * listIndices = NULL;

    // get dimensions of the input arrays   
    gridDim = fillDimensionTable(grid);
    startDim = fillDimensionTable(startVec);
    incDim = fillDimensionTable(incVec);
    sourceVecDim = fillDimensionTable(sourceVec);
    targetVecDim = fillDimensionTable(targetVec);
    
    //Check validity of input arguments
    if( grid -> nd != 3) PyErr_SetString(PyExc_ValueError,"Error in singleRayTrace - size density grid is not 3.\n");
    numDimGrid = grid -> nd;
    if( startVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in singleRayTrace - size start vector is not 1.\n");
    if( startDim[0] != 3){ 
        sprintf(errstr,"Error in singleRayTrace - start vector must have 3 elements (%d given).\n",startDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( incVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in singleRayTrace - size increment vector is not 1.\n");
    if( incDim[0] != 3){
        sprintf(errstr,"Error in singleRayTrace - increment vector must have 3 elements (%d given).\n",incDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( sourceVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in singleRayTrace - size source vector is not 1.");
    if( sourceVecDim[0] != 3){ 
        sprintf(errstr, "Error in singleRayTrace - source vector must have 3 elements(%d given).\n",sourceVecDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( targetVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in singleRayTrace - size target vector is not 1.");
    if( targetVecDim[0] != 3){ 
        sprintf(errstr, "Error in singleRayTrace - target vector must have 3 elements(%d given).\n",sourceVecDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }// End validity check


    // assign data pointers
    gridPtr = (double *) grid -> data;
    startPtr = (double *)startVec -> data;
    incPtr = (double *) incVec -> data;
    sourceVecPtr = (double *)sourceVec -> data;
    targetVecPtr = (double *)targetVec -> data;

    
    M = gridDim[2];// Number of columns
    N = gridDim[1];// Number of rows
    Q = gridDim[0];// Number of frames
    

    // create spatial grid location vectors 
    x = (double *)malloc(sizeof(double)*M);
    if (!x)  PyErr_SetString(PyExc_ValueError,"Error - malloc in singleRayTrace - x allocation failed\n");
    y = (double *)malloc(sizeof(double)*N);
    if (!y)  PyErr_SetString(PyExc_ValueError,"Error - malloc in singleRayTrace - y allocation failed\n");
    z = (double *)malloc(sizeof(double)*Q);
    if (!z)  PyErr_SetString(PyExc_ValueError,"Error - malloc in singleRayTrace - z allocation failed\n");

    // calculate spatial grid location vectors
    for (i=0;i<M;i++) x[i] = startPtr[0] + ((double)i)*incPtr[0];
    for (j=0;j<N;j++) y[j] = startPtr[1] + ((double)j)*incPtr[1];
    for (k=0;k<Q;k++) z[k] = startPtr[2] + ((double)k)*incPtr[2];
        
    
    gridValues.matrix = gridPtr;
    gridValues.start.x = x[0];
    gridValues.start.y = y[0];
    gridValues.start.z = z[0];
    gridValues.inc.x = incPtr[0];
    gridValues.inc.y = incPtr[1];
    gridValues.inc.z = incPtr[2];
    gridValues.nbCols = M;
    gridValues.nbRows = N;
    gridValues.nbFrames = Q;

    
    p1.x = sourceVecPtr[0];
    p1.y = sourceVecPtr[1];
    p1.z = sourceVecPtr[2];
    
    p2.x = targetVecPtr[0];
    p2.y = targetVecPtr[1];
    p2.z = targetVecPtr[2];
    
    // Calculate voxel indices on the ray
    if (singleRaytraceDouble(&gridValues,p1,p2,&indices) == FAILURE)
        PyErr_SetString(PyExc_ValueError,errstr);
    
    
    // create output grid
    dimList[0] = indices.nbVoxels;
    dimList[1] = 4; 
    if(  ( listIndices = (PyArrayObject *) PyArray_FromDims(2, dimList, PyArray_DOUBLE) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in singleRayTrace - listIndices allocation failed\n");
    }
    indicesPtr = (double *)listIndices -> data;
    
    
    // Fill the PyArray structure to return
    for (k = 0 ; k < dimList[0]; k++){
         //printf("(");
        for( i = 0 ; i < dimList[1]-1; i++){
            double val = indices.indList[k][i];
            indicesPtr[i + k * dimList[1]] = val;
        }
        indicesPtr[dimList[1]-1 + k * dimList[1]] = indices.lenVoxels[k];
    }
    
    
    // Free allocated memory
    for (k = 0 ; k < dimList[0]; k++){
        
            free(indices.indList[k]);
    }
    free(indices.indList);
    free(indices.lenVoxels);

    free(gridDim);
    free(startDim);
    free(incDim);
    free(targetVecDim);
    free(sourceVecDim);
    free(x);
    free(y);
    free(z);


    return PyArray_Return(listIndices);
}

int * fillDimensionTable( PyArrayObject * array){
// fillDimensionTable is used to return a table of dimensions for a given array. For instance, let's assume array shape is (4,7,6),
// array -> nd will return 3, therefore tabDim[0] = 4, tabDim[1] = 7,tabDim[2] = 6. tabDim will be returned.
    int i;
    if ( !array) return NULL;
    if( array -> nd <= 0) return NULL;// check that the number of dimension composants is at least 1.
    int * tabDim = (int *) malloc(sizeof(int) * (array -> nd));
    if(!tabDim) PyErr_SetString(PyExc_ValueError,"Error - malloc in fillDimanesionTable (singleRayTrace) - tabDim allocation failed\n");
    for (i = 0 ; i < array -> nd; i++){
        tabDim[i] =     array -> dimensions[i];
        if (tabDim[i] <= 0) PyErr_SetString(PyExc_ValueError,"Error in fillDimensionTable - a dimension is < or = 0.\n"); 
    }

    return tabDim;
}


