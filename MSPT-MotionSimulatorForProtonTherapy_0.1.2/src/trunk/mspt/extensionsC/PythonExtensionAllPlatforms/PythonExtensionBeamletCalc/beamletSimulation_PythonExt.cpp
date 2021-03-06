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


#include "beamletSimulation_PythonExt.h"


/*Constant values used to compute the radial emittance
This is from the paper by Hong (PMB 1996)
*/
const double Afit= 0.69;
const double Bfit= 0.33;
const double Cfit =0.12085e-4;
const double Dfit =0.02275;
const double PI = 3.14159265359;

#define Use2Sigmas 0 // set to 1 to use 2 sigmas and set to 0 to use only 1 sigma.

/* ==== Set up the methods table ====================== */
static PyMethodDef _CPP_extensionMethods[] = {
    {"beamletSimul", beamletSimul, METH_VARARGS},

    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};


/* ==== Initialize the CPP functions ====================== */
// Module name must be _beamletSimulation_PythonExt in compile and linked
extern "C" { 
    void init_beamletSimulation_PythonExt()  {
        (void) Py_InitModule("_beamletSimulation_PythonExt", _CPP_extensionMethods);
        import_array();  // Must be present for NumPy.  Called first after above line.
    }
}


static PyObject * beamletSimul(PyObject *self, PyObject *args){
    /*    BeamletSimul simulates the effect of a single pencil beam on a 3D volume.
    The simulation is based on Hong et al. paper, "Pencil beam algorithm for proton dose calculations",  Phys. Med. and Biol. 41,1996, and on Szymanowski et al. paper, "Two-dimensional pencil beam scaling: an improved proton dose algorithm for heterogeneous media", Phys. Med. and Biol. 47,2002.
    Inputs:
        - rhomw : density ratio to water
        - Smw :stopping power ratio to water
        - Xr : 3D matrix representing x coordinates in the IEC Gantry coordinate system (cm).
        - Yr : 3D matrix representing y coordinates in the IEC Gantry coordinate system (cm).
        - Zr : 3D matrix representing z coordinates in the IEC Gantry coordinate system (cm).
        - dmax : maximal depth in water reached for the beam energy (cm).
        - z_dd : "depth" coordinates for the depth dose curve (cm).
        - dd : "dose" coordinates for the depth dose curve. (Gy.cm ^2 )
        - sig0 : beam size at patient entrance (cm).
        - xrshift, yrshift : x and y beam shifts (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
        - sad : source-axis distance: distance between the source and the isocenter.(cm)
        - Deff : 3D matrix representing the radiological depth.(cm)
        - ssd0 : source surface distance
        - indicesRay : indices of the voxels along the ray
        - x0y0z0IECF: beam target point in the IEC fixed coordinate system.
        - spacingIECF : voxel spacing in the IEC fixed coord system
        - invRotMat : Matrix to go from the IECG to IECF coordinate system
      
    Output: 
        3D numpy array irradiated.
    Coordinates used by Hong et al are defined in the Beam Eye View 
    */
     
    //Initialize pointers that will receive input data   
    PyArrayObject * Deff = NULL;
    PyArrayObject * rhomw = NULL; 
    PyArrayObject * Smw = NULL;
    PyArrayObject * Xr = NULL;
    PyArrayObject * Yr = NULL;
    PyArrayObject * Zr = NULL;
    PyArrayObject * dmax = NULL;
    PyArrayObject * z_dd = NULL;
    PyArrayObject * dd = NULL;
    PyArrayObject * sig0 = NULL;
    PyArrayObject * xrshift = NULL;
    PyArrayObject * yrshift = NULL;
    PyArrayObject * sad = NULL;
    PyArrayObject * ssd0 = NULL;
    PyArrayObject * indicesRay = NULL;
    PyArrayObject * x0y0z0IECF = NULL;
    PyArrayObject * spacingIECF = NULL;
    PyArrayObject * invRotMat = NULL; //Matrix to go from the IECG to IECF coordinate system
    
 
    //Parse input args:
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", 
        &PyArray_Type,&rhomw, 
        &PyArray_Type,&Smw,
        &PyArray_Type,&Xr,
        &PyArray_Type,&Yr,
        &PyArray_Type,&Zr,
        &PyArray_Type,&dmax,
        &PyArray_Type,&z_dd,
        &PyArray_Type,&dd,
        &PyArray_Type,&sig0,
        &PyArray_Type,&xrshift,
        &PyArray_Type,&yrshift,
        &PyArray_Type,&sad,
        &PyArray_Type,&ssd0,
        &PyArray_Type,&Deff,
        &PyArray_Type,&indicesRay,
        &PyArray_Type,&x0y0z0IECF,//Coord Voxel(0,0,0) in the IEC fixed coord system
        &PyArray_Type,&spacingIECF,
        &PyArray_Type,&invRotMat
        ))  return NULL;

    
    if(!rhomw) return NULL;// 3D Array
    if(!Smw) return NULL;// 3D Array
    if(!Xr) return NULL;// 3D Array
    if(!Yr) return NULL;// 3D Array
    if(!Zr) return NULL;// 3D Array
    if(!dmax) return NULL;// Array with a single value
    if(!z_dd) return NULL; // 1D Array
    if(!dd) return NULL; // 1D Array
    if(!sig0) return NULL; // Array with a single value
    if(!xrshift) return NULL; // Array with a single value
    if(!yrshift) return NULL; // Array with a single value
    if(!sad) return NULL;// Array with a single value
    if(!ssd0) return NULL;// Array with a single value
    if(!Deff) return NULL; // 3D Array
    if(!indicesRay) return NULL; // Indices of the ray
    if(!x0y0z0IECF) return NULL; // Array with 3 values
    if(!spacingIECF) return NULL; //Array with 3 values
    if(!invRotMat) return NULL; // 3X3 array
  
    
    if (rhomw->descr->type_num != NPY_FLOAT) rhomw = (PyArrayObject*)PyArray_Cast(rhomw,PyArray_FLOAT);
    if (Smw->descr->type_num != NPY_FLOAT) Smw = (PyArrayObject*)PyArray_Cast(Smw,PyArray_FLOAT);
    if (Xr->descr->type_num != NPY_FLOAT) Xr = (PyArrayObject*)PyArray_Cast(Xr,PyArray_FLOAT);
    if (Yr->descr->type_num != NPY_FLOAT) Yr = (PyArrayObject*)PyArray_Cast(Yr,PyArray_FLOAT);
    if (Zr->descr->type_num != NPY_FLOAT) Zr = (PyArrayObject*)PyArray_Cast(Zr,PyArray_FLOAT);
    if (dmax->descr->type_num != NPY_FLOAT) dmax = (PyArrayObject*)PyArray_Cast(dmax,PyArray_FLOAT);
    if (z_dd->descr->type_num != NPY_FLOAT) z_dd = (PyArrayObject*)PyArray_Cast(z_dd,PyArray_FLOAT);
    if (dd->descr->type_num != NPY_FLOAT) dd = (PyArrayObject*)PyArray_Cast(dd,PyArray_FLOAT);
    if (sig0->descr->type_num != NPY_FLOAT) sig0 = (PyArrayObject*)PyArray_Cast(sig0,PyArray_FLOAT);
    if (xrshift->descr->type_num != NPY_FLOAT) xrshift = (PyArrayObject*)PyArray_Cast(xrshift,PyArray_FLOAT);
    if (yrshift->descr->type_num != NPY_FLOAT) yrshift = (PyArrayObject*)PyArray_Cast(yrshift,PyArray_FLOAT);
    if (sad->descr->type_num != NPY_FLOAT) sad = (PyArrayObject*)PyArray_Cast(sad,PyArray_FLOAT);
    if (ssd0->descr->type_num != NPY_FLOAT) ssd0 = (PyArrayObject*)PyArray_Cast(sad,PyArray_FLOAT);
    if (Deff->descr->type_num != NPY_FLOAT) Deff = (PyArrayObject*)PyArray_Cast(Deff,PyArray_FLOAT);
    if (indicesRay->descr->type_num != NPY_INT) indicesRay = (PyArrayObject*)PyArray_Cast(indicesRay,PyArray_INT);
    if (x0y0z0IECF->descr->type_num != NPY_FLOAT) x0y0z0IECF = (PyArrayObject*)PyArray_Cast(x0y0z0IECF,PyArray_FLOAT);
    if (spacingIECF->descr->type_num != NPY_FLOAT) spacingIECF = (PyArrayObject*)PyArray_Cast(spacingIECF,PyArray_FLOAT);
    if (invRotMat->descr->type_num != NPY_FLOAT) invRotMat = (PyArrayObject*)PyArray_Cast(invRotMat,PyArray_FLOAT);
    
    
    float * dataDoseGrid, * dataRhomwGrid , * dataSmwGrid, * dataXrGrid,* dataYrGrid;
    float * dataZrGrid, * dataDmax, * dataZDDArray, *dataDDArray, * dataSig0 , * dataXrShift;
    float * dataYrShift, * dataSad ,* dataSSD0, * dataDeffGrid, * datax0y0z0IECF, * dataSpacingIECF, * dataInvRotMat;
    int * dataIndicesRay;
    
    dataRhomwGrid = (float*) rhomw -> data; 
    dataSmwGrid = (float*) Smw -> data;
    dataXrGrid = (float*) Xr -> data;
    dataYrGrid = (float*) Yr -> data;
    dataZrGrid = (float*) Zr -> data;
    dataDmax = (float*) dmax -> data;
    dataZDDArray = (float*) z_dd -> data;
    dataDDArray = (float*) dd -> data;
    dataSig0 = (float*) sig0 -> data;
    dataXrShift = (float*) xrshift -> data;
    dataYrShift = (float*) yrshift -> data;
    dataSad = (float*) sad -> data;
    dataSSD0  = (float*) ssd0 -> data;
    dataDeffGrid = (float*) Deff -> data;
    
    dataIndicesRay = (int *) indicesRay-> data;

    datax0y0z0IECF = (float *) x0y0z0IECF -> data;
    dataSpacingIECF = (float*) spacingIECF -> data;
    dataInvRotMat = (float*) invRotMat -> data;
    
    
    
    int nbFrames = rhomw -> dimensions[0];
    int nbRows = rhomw -> dimensions[1];
    int nbCols = rhomw -> dimensions[2];
    

    int dimGrid[] ={nbFrames,nbRows,nbCols} ;
    
    PyArrayObject * doseGrid = NULL;
       // create resulting grid
    if(  ( doseGrid = (PyArrayObject *) PyArray_FromDims(rhomw->nd, dimGrid, PyArray_FLOAT) ) == NULL ){
        PyErr_SetString(PyExc_ValueError,"Error - malloc in beammlet simulation - doseGrid allocation failed\n");
    }
    dataDoseGrid = (float*) doseGrid -> data; // dataDoseGrid initialized to 0
    

       
    float source[3]; // Coordinates in the IEC gantry coord system
    source[0] = 0;
    source[1] = 0;
    source[2] = dataSad[0];
    float target[3];// Coordinates in the IEC gantry coord system of the target point 
    target[0] = dataXrShift[0];
    target[1] = dataYrShift[0];
    target[2] = 0;
    float vecST[3]; // vector source-target
    vecST[0] = target[0] - source[0];
    vecST[1] = target[1] - source[1];
    vecST[2] = target[2] - source[2];
    float distST = distance( source , target);
    
    
    float orthoP[3];  // Point orthogonal projection of P on the Ray S->T
    int indOrthoP[3]; // indices of orthoP in the matrix
    float distSOrthoP;
    float pointP[3]; //Coordinates of point P
    float vecSP[3];

    float maxSig = dataSig0[0] > dataSig0[1] ? dataSig0[0] : dataSig0[1];

    float variance,variance2, ddose_z, cz,xBev,yBev,zBev,oxyz, radDepth, dotProd;
    if (Use2Sigmas == 0)
        variance2 = 1; // This allows to pass the test if (variance > 0 and variance2 > 0) in the for loops
    int idx, idxVoxelBeam;
    int lenList = z_dd -> dimensions[0];
    
    
    for(int f = 0 ; f < nbFrames ; f++){ // End - iterate through frames
        for(int r = 0 ; r < nbRows ; r++){ // End - iterate through rows
            for(int c = 0 ; c < nbCols ; c++){ // End - iterate through columns
                idx = f * nbRows * nbCols + r  * nbCols + c;
                pointP[0] = dataXrGrid[idx];
                pointP[1] = dataYrGrid[idx];
                pointP[2] = dataZrGrid[idx];
                
                
                vecSP[0] = pointP[0] - source[0];
                vecSP[1] = pointP[1] - source[1];
                vecSP[2] = pointP[2] - source[2];
                        
                dotProd = dotProduct(vecST,vecSP);
                distSOrthoP =dotProd / distST;
                if (distST > 0){
                    orthoP[0] = source[0] + (distSOrthoP/distST)*vecST[0];
                    orthoP[1] = source[1] + (distSOrthoP/distST)*vecST[1];
                    orthoP[2] = source[2] + (distSOrthoP/distST)*vecST[2];
                }
                else{
                    PyErr_SetString(PyExc_ValueError,"Source and target are the same\n");
                   
                }
                
                if (distance( orthoP , pointP) > 9 * maxSig ){ // 3* 3sigmas
                    dataDoseGrid[idx]=0;
                    continue;
                }

                findPointIndicesInGrid( indOrthoP, orthoP, dataSpacingIECF, datax0y0z0IECF,dataInvRotMat);
                if ( (indOrthoP[0] < 0 || indOrthoP[0] >= nbFrames) ||   
                        (indOrthoP[1] < 0 || indOrthoP[1] >= nbRows) ||
                        (indOrthoP[2] < 0 || indOrthoP[2] >= nbCols)){
                        dataDoseGrid[idx]=0;
                        continue;
                        
                } 
                idxVoxelBeam = (indOrthoP[0]) * nbRows * nbCols + ( indOrthoP[1])  * nbCols + (indOrthoP[2]);
                
                radDepth = dataDeffGrid[idxVoxelBeam];
                
                //Find a neighbor with radiol depth computed if the voxel found on the ray did not have a radiological depth already computed
                 if(radDepth == -1){
                    float newRadDepth = -1;
                    for ( int f = -2; f <= 2; f++){
                        if ( indOrthoP[0] + f >= 0 && indOrthoP[0] + f < nbFrames){
                            for ( int r = -2; r <= 2; r++){
                                if ( indOrthoP[1] + r >= 0 && indOrthoP[1] + r < nbRows){
                                    for ( int c = -2; c <= 2; c++){
                                        if ( indOrthoP[2] + c >= 0 && indOrthoP[2] + c < nbCols){    
                                            int newIdx = (indOrthoP[0] + f) * nbRows * nbCols + ( indOrthoP[1] + r)  * nbCols + (indOrthoP[2] + c);
                                            if ( dataDeffGrid[newIdx] != -1){
                                                newRadDepth = dataDeffGrid[newIdx];
                                                break;
                                            }
                        
                                        }
                                    }
                                    if (newRadDepth != -1) break;
                                }                       
                            }
                            if (newRadDepth != -1) break;
                        }
                    }
                    if (newRadDepth == -1){
                         PyErr_SetString(PyExc_ValueError,"No radiol depth computed\n");
                         
                    } 
                    radDepth = newRadDepth;
                 }

                
                
                //Compute radial emmittance 
                variance = radialEmmittance(radDepth , dataDmax[0], dataSig0[0]);
                if (Use2Sigmas == 1)
                    variance2 = radialEmmittance(radDepth, dataDmax[0], dataSig0[1]);
                if (variance > 0 && variance2 > 0){ //Else let dose value to 0  
                    //Compute depth-dose term
                    ddose_z = interp(radDepth, dataZDDArray, dataDDArray,lenList);
                    if (ddose_z > 0){
                        //Get distance Source to the orthogonal projection of point P on the ray  
                        zBev = distSOrthoP;
                        cz = centralAxisTerm(ddose_z, dataSSD0[0], radDepth, zBev) ;
                        if ( cz > 0){  //Else let dose value to 0
                            xBev = dataXrGrid[idx] - orthoP[0];
                            yBev = dataYrGrid[idx] - orthoP[1];
                            //Compute off-axis term
                            if (Use2Sigmas == 1)
                                oxyz = offAxisTerm2Sigmas(variance,variance2,dataXrShift[0],dataYrShift[0], xBev,yBev);
                            else
                                oxyz = offAxisTerm(variance,distance( orthoP , pointP));
                            if (oxyz > 0){ //Else let dose value to 0
                                    //Compute deposited dose according to equation (9) in Hong 1996.
                                    dataDoseGrid[idx] = cz * oxyz; 

                            } // end if oxy > 0
                            else  dataDoseGrid[idx]=0;
                        } //Enf if Cz > 0
                        else  dataDoseGrid[idx]=0;
                    }// end ddose_z > 0
                    else dataDoseGrid[idx]=0;
                }//End if variance > 0
                else  dataDoseGrid[idx]=0;
            } // End - iterate through columns
        } // End - iterate through rows
    } // End - iterate through frames


    return PyArray_Return(doseGrid);
}



float radialEmmittance( float deff, float dmax, float sig0){
    
    /*Computes the radial emittance: 
     - first part depth dependent: (beam spread with depth) due to multiple scattering inside patient body.
    The beam spread parameter is calculated using these values as:
    y_t = y_R*(Afit*(t/R)^2 + Bfit*(t/R))
    where y_R = Cfit*R^2 + Dfit*R, R = dmax (Bragg peak depth)
    This is from the paper by Hong (PMB 1996) - Figure A4.
    - second part: beam spread at patient entrance : sig0.
    Finally:
    the total radial emittance. See eq (8) in Hong 1996.
        sigTot = sqrt( sig0^2 + sqrt(y_t)^2)
    We return sigTot^2
     */

     float y_R;
     y_R = Cfit * dmax * dmax + Dfit * dmax;
     float y_t;
     y_t = y_R * (Afit * (deff/dmax) * (deff/dmax) + Bfit * (deff/dmax));
     return (y_t * y_t + sig0 * sig0);

}

float interp(float xval, float * xList, float * yList, int lenList){
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
        
    
float centralAxisTerm(float dd, float ssd0, float deff, float zBev){
    /**
     *Compute the central axis term, from Hong 1996, eq (10)
     *   C(z') = DD(deff) * ((ssd0 + deff) / z')^2
     *
     *   zBeV: z coordinate in the beam eye view coordinate system (BEV)
     */
    
    //return dd;// * ( (ssd0 + deff) / zBev ) * ( (ssd0 + deff) / zBev );
    return dd * ( (ssd0 + deff) / zBev ) * ( (ssd0 + deff) / zBev );
}


float offAxisTerm(float variance, float distPQ){
    /*
     Compute the off-axis term, from Hong 1996 eq (12)
        O_xyz = [ 1/(2*pi*sigTot^2) ] * exp( - ((x' - xbev)^2 + (y' - ybev)^2 ) / (2*sigTot^2) )
     */

    return (1/(2*PI*(variance)))*exp( -0.5*(distPQ * distPQ) / (variance));
    
    
}

float offAxisTerm2Sigmas(float varX,float varY, float xrshift,float yrshift, float xBev, float yBev){
    /*
     Compute the off-axis term, from Hong 1996 eq (12) extended to 2 dimensions:
     
        O_xyz = [ (1/(4*pi*sigTotX^2)) + 1/(4*pi*sigTotY^2) ] * exp( - ( ((x' - xbev)^2 / (2*sigTotX^2))  + ((y' - ybev)^2 ) / (2*sigTotY^2)) )
     */
     
    PyErr_SetString(PyExc_ValueError,"Beamlet simulation with 2 sigmas not yet implemented\n"); 
//     float xshiftBev =  xrshift;
    float yshiftBev;
    yshiftBev = yrshift;
    float deltaX =  xBev;
    float deltaY = yBev;
    
    float factor = 1/(2*PI*sqrt(varX)*sqrt(varY));
    float xRatio = (deltaX * deltaX) / (2*varX);
    float yRatio = (deltaY * deltaY) / (2*varY);
    
    return factor * exp( -(xRatio + yRatio) );
    
}

float dotProduct(float * p1,float * p2){
    /*
    Dot product
    */
    float dotProd;
    dotProd = (p1[0]*p2[0] + p1[1] * p2[1] + p1[2] * p2[2]);
    return dotProd;

}
 
float distance( float * p1, float *p2){
    /*
    Euclidean distance
    */
    return sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]));

}



void findPointIndicesInGrid( int * ind, float * pointP1, float * spacing, float * x0y0z0,float * rotMat){    
    /*
    Convert coordinates of a rotated point to indices in a grid  
    */    
    
    float p1IECF[3];
    getRotatedPoint( p1IECF , pointP1, rotMat);
    
    ind[0] = round( (p1IECF[1]-x0y0z0[1])/spacing[1]);// y IEC Fixed is along the frame axis
    ind[1] = round( (p1IECF[2]-x0y0z0[2])/spacing[2]);// z IEC Fixed -> rows
    ind[2] = round( (p1IECF[0]-x0y0z0[0])/spacing[0]);
}


void getRotatedPoint( float * newP , float * P, float * rotMat){
    /*
    Rotate a point coordinates based on a rotation matrix 
    */
    
    newP[0] = rotMat[0] * P[0] + rotMat[1] * P[1] + rotMat[2] * P[2];
    newP[1] = rotMat[3] * P[0] + rotMat[4] * P[1] + rotMat[5] * P[2];
    newP[2] = rotMat[6] * P[0] + rotMat[7] * P[1] + rotMat[8] * P[2];

}
