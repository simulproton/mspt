/*
 ########################################################################
#
#  fewRays_RadDepth_PythonExt.cpp
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# 22 May 2014
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
 */
#include "fewRays_RadDepth_PythonExt.h"
char errstr[200];  // error string that all routines have access to

/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"fewRaysRadDepth", fewRaysRadDepth, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _fewRays_RadDepth_PythonExt in compile and linked
extern "C" { 
    void init_fewRays_RadDepth_PythonExt()  {
        (void) Py_InitModule("_fewRays_RadDepth_PythonExt", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * fewRaysRadDepth(PyObject *self, PyObject *args){
    
    /*
    This function is the python extension. It finds the voxels of a volume encountered by a ray and their radiological depth 
    for few target spots:
    Inputs:
        grid : 3D numpy array representing the volume: density or relative stopping power (depending if this is photons or protons) 
        startVec : numpy array (1,3) : Coordinates of the first voxel (0,0,0)
        incVec : numpy array (1,3): Increment vector (spacing) when incrementing (+1,+1,+1) the voxels indices.
        sourceVec : numpy array (1,3): Coordinates of the source of the ray
        targetVec : numpy array (1,3): Coordinates of the target of the ray // This should be the "central" target
        auxilaryTargetVecs : numpy array (nbTargets,3): Coordinates of the targets around the "central" target
    
    Output:
        Numpy array (3D) in which voxels encountered by the ray have a value of radiological depth, and -1 otherwise
    */
     
    
    //Initialize pointers that will receive input data
    PyArrayObject * grid = NULL;
    PyArrayObject * startVec = NULL;
    PyArrayObject * incVec = NULL;
    PyArrayObject * sourceVec = NULL;
    PyArrayObject * targetVec = NULL;
    PyArrayObject * auxilaryTargetVecs = NULL;
    
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!", 
        &PyArray_Type,&grid, 
        &PyArray_Type,&startVec, 
        &PyArray_Type,&incVec, 
        &PyArray_Type,&sourceVec,
        &PyArray_Type,&targetVec,
        &PyArray_Type,&auxilaryTargetVecs))  return NULL;

    if(!grid) return NULL;
    if(!startVec) return NULL;
    if(!incVec) return NULL;
    if(!sourceVec) return NULL;
    if(!targetVec) return NULL;
    if(!auxilaryTargetVecs) return NULL;
    
    //Cast values to float values: numpy creates array with float values
    if (grid->descr->type_num != NPY_FLOAT) grid = (PyArrayObject*)PyArray_Cast(grid,PyArray_FLOAT);
    if (startVec->descr->type_num != NPY_FLOAT) startVec = (PyArrayObject*)PyArray_Cast(startVec,PyArray_FLOAT);
    if (incVec->descr->type_num != NPY_FLOAT) incVec = (PyArrayObject*)PyArray_Cast(incVec,PyArray_FLOAT);
    if (sourceVec->descr->type_num != NPY_FLOAT) sourceVec = (PyArrayObject*)PyArray_Cast(sourceVec,PyArray_FLOAT);
    if (targetVec->descr->type_num != NPY_FLOAT) targetVec = (PyArrayObject*)PyArray_Cast(targetVec,PyArray_FLOAT);
    if (auxilaryTargetVecs->descr->type_num != NPY_FLOAT) auxilaryTargetVecs = (PyArrayObject*)PyArray_Cast(auxilaryTargetVecs,PyArray_FLOAT);


    // scalars
    int i,j,k,M,N,Q;

    // pointers to vectors that do not need to be freed
    int *gridDim, *startDim, *incDim ,*sourceVecDim, *targetVecDim,*auxTargetVecDim;
    float *gridPtr,*startPtr,*incPtr, *sourceVecPtr, *targetVecPtr,*auxTargetVecPtr,*deffPtr;
    int * countPtr;

    // pointers to vectors that need to be freed
    float *x,*y,*z;
    
    // other
    int numDimGrid = -1;//Number of dimensions of the density grid - 3 is expected
    
    // point structures that will be used for raytrace operation
    POINTXYZ p1;
    POINTXYZ p2;

    // grid structures for raytracing operation
    VOLUME deff;
    VOLUME gridValues;

    //Output data 
    PyArrayObject * deffGrid = NULL;

    // get dimensions of the input arrays   
    gridDim = fillDimensionTable(grid);
    startDim = fillDimensionTable(startVec);
    incDim = fillDimensionTable(incVec);
    sourceVecDim = fillDimensionTable(sourceVec);
    targetVecDim = fillDimensionTable(targetVec);
    auxTargetVecDim = fillDimensionTable(auxilaryTargetVecs);
    
    //Check validity of input arguments
    if( grid -> nd != 3) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size density grid is not 3.\n");
    numDimGrid = grid -> nd;
    if( startVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size start vector is not 1.\n");
    if( startDim[0] != 3){ 
        sprintf(errstr,"Error in fewRaysRadDepth - start vector must have 3 elements (%d given).\n",startDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( incVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size increment vector is not 1.\n");
    if( incDim[0] != 3){
        sprintf(errstr,"Error in fewRaysRadDepth - increment vector must have 3 elements (%d given).\n",incDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( sourceVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size source vector is not 1.");
    if( sourceVecDim[0] != 3){ 
        sprintf(errstr, "Error in fewRaysRadDepth - source vector must have 3 elements(%d given).\n",sourceVecDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( targetVec -> nd != 1) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size target vector is not 1.");
    if( targetVecDim[0] != 3){ 
        sprintf(errstr, "Error in fewRaysRadDepth - target vector must have 3 elements(%d given).\n",targetVecDim[0]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }
    if( auxilaryTargetVecs -> nd != 2) PyErr_SetString(PyExc_ValueError,"Error in fewRaysRadDepth - size auxiliary target vectors is not 2.");
    if( auxTargetVecDim[1] != 3){ 
        sprintf(errstr, "Error in fewRaysRadDepth - auxiliary target vectors must have 3 elements per target(%d given).\n",auxTargetVecDim[1]);
        PyErr_SetString(PyExc_ValueError,errstr);
    }// End validity check

    // assign data pointers
    gridPtr = (float *) grid -> data;
    startPtr = (float *)startVec -> data;
    incPtr = (float *) incVec -> data;
    sourceVecPtr = (float *)sourceVec -> data;
    targetVecPtr = (float *)targetVec -> data;
    auxTargetVecPtr = (float *)auxilaryTargetVecs -> data;

    M = gridDim[2];// Number of columns
    N = gridDim[1];// Number of rows
    Q = gridDim[0];// Number of frames
    
    countPtr = (int *) malloc(sizeof(int)*M*N*Q);
    for ( int i = 0; i < M*N*Q; i++) countPtr[i] = 0;


    // create spatial grid location vectors 
    x = (float *)malloc(sizeof(float)*M);
    if (!x)  PyErr_SetString(PyExc_ValueError,"Error - malloc in fewRaysRadDepth - x allocation failed\n");
    y = (float *)malloc(sizeof(float)*N);
    if (!y)  PyErr_SetString(PyExc_ValueError,"Error - malloc in fewRaysRadDepth - y allocation failed\n");
    z = (float *)malloc(sizeof(float)*Q);
    if (!z)  PyErr_SetString(PyExc_ValueError,"Error - malloc in fewRaysRadDepth - z allocation failed\n");

    // calculate spatial grid location vectors
    for (i=0;i<M;i++) x[i] = startPtr[0] + ((float)i)*incPtr[0];
    for (j=0;j<N;j++) y[j] = startPtr[1] + ((float)j)*incPtr[1];
    for (k=0;k<Q;k++) z[k] = startPtr[2] + ((float)k)*incPtr[2];


    // create output grid
    if(  ( deffGrid = (PyArrayObject *) PyArray_FromDims(grid -> nd, gridDim, PyArray_FLOAT) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in fewRaysRadDepth - Radiol. depth Grid allocation failed\n");
    }
    deffPtr = (float *)deffGrid -> data;
    
    //initialize grids
    for (i=0;i<M*N*Q;i++) deffPtr[i] = -1.0;
        
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


    deff.matrix = deffPtr;
    deff.start.x = x[0];
    deff.start.y = y[0];
    deff.start.z = z[0];
    deff.inc.x = incPtr[0];
    deff.inc.y = incPtr[1];
    deff.inc.z = incPtr[2];
    deff.nbCols = M;
    deff.nbRows = N;
    deff.nbFrames = Q;




    
    p1.x = sourceVecPtr[0];
    p1.y = sourceVecPtr[1];
    p1.z = sourceVecPtr[2];
    
    p2.x = targetVecPtr[0];
    p2.y = targetVecPtr[1];
    p2.z = targetVecPtr[2];

    p2.x = p1.x + (p2.x-p1.x)*(float)100.0;
    p2.y = p1.y + (p2.y-p1.y)*(float)100.0;
    p2.z = p1.z + (p2.z-p1.z)*(float)100.0;
    //Fill radiological depth for first central target
    if ( calcRadiolDepth(&gridValues,&deff,p1,p2,countPtr) == FAILURE)
                 PyErr_SetString(PyExc_ValueError,errstr);
    
    //Fill radiological depth for auxiliary targets
    for ( int k = 0; k <  auxTargetVecDim[0]; k++){
        p2.x = auxTargetVecPtr[k*3 +0];
        p2.y = auxTargetVecPtr[k*3 +1];
        p2.z = auxTargetVecPtr[k*3 +2];

        p2.x = p1.x + (p2.x-p1.x)*(float)100.0;
        p2.y = p1.y + (p2.y-p1.y)*(float)100.0;
        p2.z = p1.z + (p2.z-p1.z)*(float)100.0;

        if ( calcRadiolDepth(&gridValues,&deff,p1,p2,countPtr) == FAILURE)
                 PyErr_SetString(PyExc_ValueError,errstr);    
    
    }

    for (i=0;i<M*N*Q;i++){
        double count = countPtr[i];
        if (count > 1) deffPtr[i] = deffPtr[i] / count;
    }


    free(gridDim);
    free(startDim);
    free(incDim);
    free(sourceVecDim);
    free(targetVecDim );
    free(auxTargetVecDim );

    free(x);
    free(y);
    free(z);

    return PyArray_Return(deffGrid);
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


