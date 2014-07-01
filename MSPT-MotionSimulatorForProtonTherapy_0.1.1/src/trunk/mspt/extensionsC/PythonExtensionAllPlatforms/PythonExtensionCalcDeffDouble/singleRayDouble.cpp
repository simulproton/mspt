/*
 ########################################################################
#
#  singleRayDouble.cpp
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# April, 30 2013.
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
  
 
 This code is an implementation of the method described by Robert L. Siddon in "Fast calculation of the exact radiological path for three-dimensionnal CT array", Am. Assoc. Phys. Med. 1985.
 
 It traces a ray from point1 (start) to point2 (end) and assign to all the voxel crossed their radiological depth.
 
 radiological_depth_grid (shape: Nx cols ,Ny rows,Nz frames): matrix that will be filled with radiological depths along the ray. It should be initialized with 0 when radiological depth has not been computed. If a value other than 0 is already present it means the user will have to devide each value by the number of rays traced.
 
 stpPwr_grid (shape: Nx cols ,Ny rows,Nz frames): matrix storing the electron density of each voxel
 
 Points are defined in 3D : (x,y,z)
 
 point1 should be the source
 point2 the "target" point 
 
 */
#include "singleRayDouble.h"

int singleRaytraceDouble(VOLUME * grid,
             POINTXYZ point1,POINTXYZ point2 , INDICES * indRay){
    /*
     Main function of the process
     
     stpPwr_grid : grid (VOLUME) storing the electron density (radiol. depth. for X-Rays) or relative stopping power(for protons)
     radiological_depth_grid : grid (VOLUME) that will be filled with the radiological depths.
     point1 : source of the ray
     point2 : target of the ray
     indRay : list of indices encountered by the ray
     
     */
     
     if (point1.x == point2.x) point2.x = point2.x + 1e-4;
     if (point1.y == point2.y) point2.y = point2.y + 1e-4;
     if (point1.z == point2.z) point2.z = point2.z + 1e-4;
     
     
    double alphaMinMaxX[2];
    double alphaMinMaxY[2];
    double alphaMinMaxZ[2];
    double alphaMinMax[2];
    
    double coordPlanesX[2];// ind 0: first plane, ind 1: last plane
    double coordPlanesY[2];
    double coordPlanesZ[2];
    
    int nX = grid->nbCols + 1;
    int nY = grid->nbRows + 1;
    int nZ = grid->nbFrames + 1;
    
    
    double xSpacing = (grid->inc).x;
    double ySpacing = (grid->inc).y;
    double zSpacing = (grid->inc).z;
    
    int iMinMax[2];
    int jMinMax[2];
    int kMinMax[2];
    
    double startX = (grid->start).x;
    double startY = (grid->start).y;
    double startZ = (grid->start).z;
    
    
    double d12 =  euclideanDistanceDouble( point1,point2);
    
 
 
    // Step 1 compute coord first and last planes:
    coordPlanesX[0] = coordPlane(1,startX - 0.5*xSpacing,xSpacing);
     coordPlanesX[1] = coordPlane(nX,coordPlanesX[0] , xSpacing );
    
     coordPlanesY[0] = coordPlane(1,startY- 0.5*ySpacing,ySpacing);
     coordPlanesY[1] = coordPlane(nY,coordPlanesY[0] , ySpacing );

    coordPlanesZ[0] = coordPlane(1,startZ- 0.5*zSpacing,zSpacing);
     coordPlanesZ[1] = coordPlane(nZ,coordPlanesZ[0] , zSpacing );

    // Step 2 Find alpha_min alpha_max for each axis direction (x,y,z):
    int sameX,sameY,sameZ;
    sameX = computeAlpha(1, point1.x, point2.x,coordPlanesX[0], xSpacing , &alphaMinMaxX[0]);// if sameX == -1 , X2-X1 == 0
    sameX = computeAlpha(nX, point1.x, point2.x,coordPlanesX[0], xSpacing , &alphaMinMaxX[1]);// if sameX == -1 , X2-X1 == 0
    sortAlphaMinMax(alphaMinMaxX);
    
    sameY = computeAlpha(1, point1.y, point2.y,coordPlanesY[0], ySpacing , &alphaMinMaxY[0]);// if sameY == -1 , Y2-Y1 == 0
    sameY = computeAlpha(nY, point1.y, point2.y,coordPlanesY[0], ySpacing , &alphaMinMaxY[1]);// if sameY == -1 , Y2-Y1 == 0
    sortAlphaMinMax(alphaMinMaxY);
    
    sameZ = computeAlpha(1, point1.z, point2.z,coordPlanesZ[0], zSpacing , &alphaMinMaxZ[0]);// if sameZ == -1 , Z2-Z1 == 0
    sameZ = computeAlpha(nZ, point1.z, point2.z,coordPlanesZ[0], zSpacing , &alphaMinMaxZ[1]);// if sameZ == -1 , Z2-Z1 == 0
    sortAlphaMinMax(alphaMinMaxZ);

    //Step 3 Find overal alpha_min alpha_max
    findAlphaMinMax(alphaMinMaxX,alphaMinMaxY,alphaMinMaxZ,alphaMinMax);
        
    //Step 4 find i_min, i_max,j_min, j_max,k_min, k_max:
    findIndMinMax(coordPlanesX[0],coordPlanesX[1],xSpacing,point1.x, point2.x, nX, alphaMinMax , iMinMax);
    findIndMinMax(coordPlanesY[0],coordPlanesY[1],ySpacing,point1.y, point2.y, nY, alphaMinMax , jMinMax);
    findIndMinMax(coordPlanesZ[0],coordPlanesZ[1],zSpacing,point1.z, point2.z, nZ, alphaMinMax , kMinMax);

    
    //Steap 5 Calculate parametric sets
    double * listAlphaX, * listAlphaY, *listAlphaZ;
    int lenListX, lenListY, lenListZ;
    lenListX = iMinMax[1]-iMinMax[0]+1;
    lenListY = jMinMax[1]-jMinMax[0]+1;
    lenListZ = kMinMax[1]-kMinMax[0]+1;

    listAlphaX = (double *)malloc(sizeof(double) * lenListX);
    if(!listAlphaX){
        perror("Malloc alpha X error... \n");
        return FAILURE;
    }
    
    listAlphaY = (double *)malloc(sizeof(double) * lenListY);
    if(!listAlphaY){
        perror("Malloc alpha Y error... \n");
       return FAILURE;
    }
    listAlphaZ = (double *)malloc(sizeof(double) * lenListZ);
    if(!listAlphaZ){
        perror("Malloc alpha Z error... \n");
        return FAILURE;
    }
    
    if( calculateParametricSets(iMinMax ,point1.x, point2.x, coordPlanesX[0], xSpacing, listAlphaX,alphaMinMax ) != sameX){
        perror("Error SameX should be equal to value returned by calculateParamtricSet\n");
        return FAILURE;
    }
    if( calculateParametricSets(jMinMax ,point1.y, point2.y, coordPlanesY[0], ySpacing, listAlphaY,alphaMinMax) != sameY){
        perror("Error SameY should be equal to value returned by calculateParamtricSet\n");
        return FAILURE;
    }
    if( calculateParametricSets(kMinMax ,point1.z, point2.z, coordPlanesZ[0], zSpacing, listAlphaZ,alphaMinMax) != sameZ){
        perror("Error SameZ should be equal to value returned by calculateParamtricSet\n");
       return FAILURE;
    }

    
    // Step 6 merge alpha sets:
    int lenAlpha = ( 2 + lenListX + lenListY + lenListZ);
    double * listAlpha = (double *) malloc(sizeof(double) * lenAlpha);
    if(!listAlpha){
        perror("Error - memory allocation list alphas\n");
        return FAILURE;
    }
    mergeSets(alphaMinMax, listAlphaX, lenListX, listAlphaY, lenListY , listAlphaZ, lenListZ, listAlpha,d12, ((xSpacing+ySpacing+zSpacing)/3.0));
   
   
   // Step 7 calculate voxel indices

    POINTXYZ plane1;
    plane1.x = coordPlanesX[0];
    plane1.y = coordPlanesY[0];
    plane1.z = coordPlanesZ[0];
    
    POINTXYZ spacing;
    spacing.x = xSpacing;
    spacing.y = ySpacing;
    spacing.z = zSpacing;
    int nbIndices = lenAlpha - 1;
    int ** indList = (int **)malloc(sizeof(int *) * nbIndices);
    if(!indList){
        perror("Error - memory allocation list ind\n");
        return FAILURE;
    }   
 
    for (int i = 0 ; i <  nbIndices; i++){
        indList[i] = (int *) malloc(sizeof(int ) * 3); 
        if(! (indList[i])){
            perror("Error - memory allocation list ind\n");
            return FAILURE;
        }
    }
    int ** indMinMax = (int **)malloc(sizeof(int *) * 3);
    if(! indMinMax){
        perror("Error - memory allocation indMinMax\n");
        return FAILURE;
    }
    indMinMax[0] = iMinMax;
    indMinMax[1] = jMinMax;
    indMinMax[2] = kMinMax;
    
    calculateVoxelIndices(point1, point2, plane1, spacing, listAlpha, lenAlpha, indList, indMinMax);
    (*indRay).nbVoxels = 0;
    int i,j,k;
    for ( int m= 1; m < lenAlpha; m++){
        if ((listAlpha[m] - listAlpha[m-1]) > 0){
            if (indList[m-1][0] > 0 && indList[m-1][1] > 0 && indList[m-1][2] > 0 ){
                i = indList[m-1][0]-1;
                j = indList[m-1][1]-1;
                k = indList[m-1][2]-1;
                if ( i < nX-1 && j < nY-1 && k < nZ-1){
                    (*indRay).nbVoxels++;
                }
            }
        }
    }
    
    //Fill the list of indices
    (*indRay).indList = (int **)malloc(sizeof(int *) * (*indRay).nbVoxels);
    if(!indList){
        perror("Error - memory allocation list ind\n '(*indRay).indList = (int **)malloc(sizeof(int *) * (*indRay).nbVoxels)'\n");
        return FAILURE;
    }   
    for (int i = 0 ; i <  (*indRay).nbVoxels; i++){
        (*indRay).indList[i] = (int *) malloc(sizeof(int ) * 3); 
        if(! ((*indRay).indList[i])){
            perror("Error - memory allocation list ind\n");
            return FAILURE;
        }
    }
    
    (*indRay).lenVoxels = (double *)malloc(sizeof(double) * (*indRay).nbVoxels);
    if(!( (*indRay).lenVoxels )){
        perror("Error - memory allocation list dist\n '(*indRay).lenVoxels = (double *)malloc(sizeof(double) * (*lenVoxels).nbVoxels)'\n");
        return FAILURE;
    } 
    
    int  count = 0;
    for ( int m= 1; m < lenAlpha; m++){
        if ((listAlpha[m] - listAlpha[m-1]) > 0){
           if (indList[m-1][0] > 0 && indList[m-1][1] > 0 && indList[m-1][2] > 0 ){
                i = indList[m-1][0]-1;
                j = indList[m-1][1]-1;
                k = indList[m-1][2]-1;
                if ( i < nX-1 && j < nY-1 && k < nZ-1){
                    // Voxels indices
                    (*indRay).indList[count][0]=k;
                    (*indRay).indList[count][1]=j;
                    (*indRay).indList[count][2]=i;
                    // Length of the ray in the voxel
                    (*indRay).lenVoxels [count] = d12 * (listAlpha[m] - listAlpha[m-1]);
                    count ++;
                }
            }
        }
    }
    
    
    free(listAlphaX);
    free(listAlphaY);
    free(listAlphaZ);
    free(listAlpha);
    free(indMinMax);
    for (int i = 0 ; i <  nbIndices; i++)
        free(indList[i]); 
    free(indList);
    return SUCCESS;
}

