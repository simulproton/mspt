/*
########################################################################
#
#  fillStopPwrArray_PythonExt.cpp
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
########################################################################
*/ 
#include "fillStopPwrArray_PythonExt.h"

#define BONE_HIGH 20 // Density max from CT numbers 

char errstr[200];  // error string that all routines have access to

/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"fillStopPwrGrid", fillStopPwrGrid, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _fillStopPwrArray_PythonExt in compile and linked
extern "C" { 
    void init_fillStopPwrArray_PythonExt()  {
        (void) Py_InitModule("_fillStopPwrArray_PythonExt", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * fillStopPwrGrid(PyObject *self, PyObject *args){
/*
 Function that converts the density grid into relative stopping power grid, using a conversion table
 Inputs:
    densityGrid : 3D numpy array in which each voxel is represented by its density
    conversionTable : table used for the conversion. Array of 3 rows and n columns. 
                        Row 0 : lowest bound of the density values for each media, 
                        Row 1 : upper bound of the density values for each madia , 
                        Row 3: corresponding relative mass Stopping Power value. 
                        Each column corresponds to 1 type of medium: 'air','adipose,','muscle','bone' ...
                        sorted by increasing density
    tableCorrectionFactor : Table used to correct the mass stopping power.
                        Row 0: Densities
                        Row 1: Correction factor
 
 Output: 
    3D numpy array in which the relative stopping power is stored.
 
 */ 

    //Initialize pointers that will receive input data
    PyArrayObject * densityGrid = NULL; 
    PyArrayObject * conversionTable = NULL;  
    PyArrayObject * correcTable = NULL;  
 
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!", 
        &PyArray_Type,&densityGrid, 
        &PyArray_Type,&conversionTable,
        &PyArray_Type,&correcTable))  return NULL;
        
    if(!densityGrid) return NULL;
    if(!conversionTable) return NULL;
    if(!correcTable) return NULL;

    
    //Cast values to float values: numpy creates array with float values
    if (densityGrid->descr->type_num != NPY_FLOAT) densityGrid = (PyArrayObject*)PyArray_Cast(densityGrid,PyArray_FLOAT);
    if (conversionTable->descr->type_num != NPY_FLOAT) conversionTable = (PyArrayObject*)PyArray_Cast(conversionTable,PyArray_FLOAT);
    if (correcTable->descr->type_num != NPY_FLOAT) correcTable = (PyArrayObject*)PyArray_Cast(correcTable,PyArray_FLOAT);


    
    float  * dataStpPwrGrid, * dataDensGrid, * dataConversion, * dataCorrecTable;
    
    dataDensGrid = (float*) densityGrid -> data;  
    dataConversion = (float*) conversionTable -> data; 
    dataCorrecTable = (float*) correcTable -> data; 
   
    int nbColsConvTab = conversionTable->dimensions[1];
  
    int nbFramesDens = densityGrid->dimensions[0];
    int nbRowsDens = densityGrid->dimensions[1];
    int nbColsDens = densityGrid-> dimensions[2];
    int dimGrid[] ={nbFramesDens,nbRowsDens,nbColsDens} ;
    
    int shapeCorrecTable[2];
    shapeCorrecTable[0] = correcTable->dimensions[0];
    shapeCorrecTable[1] = correcTable->dimensions[1];
    
      
    // create output grid
    PyArrayObject * stpPwrGrid = NULL;
    if(  ( stpPwrGrid = (PyArrayObject *) PyArray_FromDims(densityGrid->nd, dimGrid, PyArray_FLOAT) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in fillDensityGrid - stpPwrGrid allocation failed\n");
    }
    dataStpPwrGrid = (float*) stpPwrGrid -> data;
    
    float * slope, * intercept;
    slope = (float *) malloc(sizeof(float) * (nbColsConvTab) );
    if(!slope) {
        PyErr_SetString(PyExc_ValueError,"Error allocation slope tab in fillStpPwr\n");
        
    }
    intercept = (float *) malloc(sizeof(float) * (nbColsConvTab) );
    if(!intercept) {
        PyErr_SetString(PyExc_ValueError,"Error allocation intercept tab in fillStpPwr\n");
        
    }
    for(int i = 0; i < nbColsConvTab; i++){
        slope[i] = FLT_MAX;
        intercept[i] = FLT_MAX;
    }
    
    float densValue, densValueMass;
    float convertedValue;
    int idx;
    for(int f = 0 ; f < nbFramesDens ; f++){ // End - iterate through frames
        for(int r = 0 ; r < nbRowsDens ; r++){ // End - iterate through rows
            for(int c = 0 ; c < nbColsDens ; c++){ // End - iterate through columns
                idx = f * nbColsDens * nbRowsDens + r * nbColsDens + c;
                densValueMass = dataDensGrid[idx];
                densValue =densValueMass;
                convertedValue = 0;
                if (densValue < dataConversion[0]){
                    printf("Error- fill Stop Pwr Ratio array - density value < airLow ( 0 ): densValue = %f\n",densValue);
                    PyErr_SetString(PyExc_ValueError,"Error density value in fill Stop Pwr Ratio array...");
                                
                }
                if (densValue > BONE_HIGH){
                     printf("Error- fill Stop Pwr Ratio array - density value > boneHigh ( %i ): densValue = %f\n",BONE_HIGH,densValue);
                    PyErr_SetString(PyExc_ValueError,"Error density value in fill Stop Pwr Ratio array...");
                                                   
                }
                for( int i = 0 ; i< nbColsConvTab ; i++){
                    if(slope[i] == FLT_MAX){
                        if(i==0){
                            slope[i] =  0;
                        }
                        else{
                            slope[i] = (dataConversion[2 * nbColsConvTab + i-1 ] - dataConversion[2 * nbColsConvTab + i ])/ (dataConversion[i] - dataConversion[1 * nbColsConvTab + i ]);
                        }
                    }                        
                    if(intercept[i] == FLT_MAX){
                        if(i==0){
                            intercept[i] = dataConversion[2 * nbColsConvTab + i ]; 
                        }
                        else{
                            intercept[i] = dataConversion[2 * nbColsConvTab + i ] - dataConversion[1 * nbColsConvTab + i ] * slope[i];
                        }                            
                    }
                    if(i < nbColsConvTab - 1)
                    {
                        if (densValue >= dataConversion[i] && densValue < dataConversion[1 * nbColsConvTab + i ]){
                            convertedValue = densValue * slope[i] + intercept[i];
                            break;
                        }
                    }
                    else{
                        if (densValue >= dataConversion[i] && densValue <= dataConversion[1 * nbColsConvTab + i ]){
                            convertedValue = densValue * slope[i] + intercept[i];
                            break;
                        }
                        if(densValue > dataConversion[1 * nbColsConvTab + i ] && densValue < BONE_HIGH)
                        {
                            convertedValue = densValue * slope[i] + intercept[i] ;
                            break;
                        }
                        else{
                            PyErr_SetString(PyExc_ValueError,"Dens value >= 20! \n");
                           
                        }
                    }
                }
            dataStpPwrGrid[ idx ] = convertedValue * densValueMass * scalingFactor(densValueMass,dataCorrecTable,shapeCorrecTable); // if converted value is mass stopping power ratio
            if (convertedValue < 0 ){
                    printf("Before exit: idx = %i, converted Value= %f , Dens Value = %f , minDensity = %f , max Density = %f \n",idx,convertedValue,densValue,dataConversion[0],dataConversion[1 * nbColsConvTab+nbColsConvTab-1]);
                    PyErr_SetString(PyExc_ValueError,"Wrong Smw value < 0 ");
                    
                }
            } // End - iterate through columns
        } // End - iterate through rows
    } // End - iterate through frames
    free(slope);
    free(intercept);
    return PyArray_Return(stpPwrGrid);
}


float scalingFactor(float density , float * factCorrecTable, int * tableShape){
    /*
        density : density in g.cm-3
        factCorrecTable : table containing the correction factors
        tableShape: shape (rows , columns) of factCorrecTable. It should be (2,n)

    */
    

    if(density <= 0 || density > BONE_HIGH){
        PyErr_SetString(PyExc_ValueError,"Density < 0 density > BONE_HIGH or in fillStopPwrArray\n");
        
    }
//     float  lDens[] = {0.01,	0.05,	0.1,	0.15,	0.2,	0.25,	0.3,	0.35,	0.4,	0.45,	0.5,	0.825,	0.85,	0.857,	0.9,	0.925,	0.95,	0.975,  0.976,	1.0,	1.025,	1.05,	1.075,	1.1,	1.14,	1.16,	1.2,	1.25,	1.3,	1.35,	1.4,	1.45,	1.5};
//     float  lFact[] = {15.0,	5.48,	1.9,	1.7,	1.4,	1.28,	1.22,	1.18,	1.16,	1.17,	1.11,	1.12,	1.04,	1.02,	1.04,	1.019,	1.03,	1.05,	1.0,    1.0,	1.0,	1.06,	1.03,	1.03,	1.02,	1.02,	1.02,	1.03,	1.04,	1.04,	1.05,	1.05,	1.04};
//     int lenList = 33;

    if (factCorrecTable && tableShape){ 
        int lenList = tableShape[1]; // Number of columns
        float * lDens = factCorrecTable;
        float * lFact = factCorrecTable + lenList;
        return interpList(density,lDens, lFact , lenList);
    }
    return 1.0;
    
}






float interp( float x1, float y1, float x2, float y2, float val){

    float slope = (y1 - y2)/(x1-x2);
    float intercept = y2 - x2 * slope;
    return slope * val + intercept;


}


float interpList(float xval, float * xList, float * yList, int lenList){
    /* Piecewise linear interpolation function:
     Function that interpolates a function defined by (xList,yList) and returns 'yval' corresponding to 'xval' for the interpolated function.
     */

    float x0 = xList[0];
    int idx=-1;
    //If xval <= x0 return y0 , if xval >= x_end return y_end
    if( xval <= x0) return yList[0];
    if (xval >= xList[lenList -1]) return yList[lenList -1];

    //Find place of xval in xList : find index
    //Search by dichotomy
    int left = 0, right = lenList -1;
    int middle;
    while( left <= right){
        middle = (left + right)/2;
        if (xList[middle] <= xval && xList[middle+1] > xval){
          idx = middle;
          break;
        } 
        if( xList[middle] >  xval){
            right = middle;
        }
        else left = middle+1;
    }
    if (idx == -1){
        PyErr_SetString(PyExc_ValueError,"Error - index not found in xList for linear search\n");
       
    }

    float slope = (yList[idx+1] - yList[idx]) /  (xList[idx+1] - xList[idx]);
    
    return slope * (xval - xList[idx]) + yList[idx]; 
}
