/*
 ########################################################################
#
#  calcRadiolDepth.cpp
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# April, 30 2013.
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
  
 
 This code is an implementation of the method described by Robert L. Siddon in "Fast calculation of the exact radiological path for three-dimensionnal CT array", Am. Assoc. Phys. Med. 1985.
 
 It traces a ray from point1 (start) to point2 (end) and assign to all the voxel crossed their radiological depth.
 
 radiological_depth_grid (shape: Nx cols ,Ny rows,Nz frames): matrix that will be filled with radiological depths along the ray. It should be initialized with 0 when radiological depth has not been computed. If a value other than 0 is already present it means the user will have to divide each value by the number of rays traced.
 
 stpPwr_grid (shape: Nx cols ,Ny rows,Nz frames): matrix storing the electron density of each voxel
 
 Points are defined in 3D : (x,y,z)
 
 point1 should be the source
 point2 the "target" point 
 
 */
#include "calcRadiolDepth.h"


int calcRadiolDepth(VOLUME *stpPwr_grid,VOLUME *radiological_depth_grid,
             POINTXYZ point1,POINTXYZ point2 , int * tabCount){
    /*
     Main function of the process
     
     stpPwr_grid : grid (VOLUME) storing the electron density (radiol. depth. for X-Rays) or relative stopping power(for protons)
     radiological_depth_grid : grid (VOLUME) that will be filled with the radiological depths.
     point1 : source of the ray
     point2 : target of the ray
     tabCount : count how many rays are going through each voxels in order to perform an average afterwards.
     
     */
     
     
     
    if ( absolute(point1.x - point2.x) <= 1e-5 ){
        point2.x = point2.x + 1e-3;
    }
    if (absolute(point1.y - point2.y) <= 1e-5 ){
        point2.y = point2.y + 1e-3;
    }
    if (absolute(point1.z - point2.z) <= 1e-5 ){
        point2.z = point2.z + 1e-3;
    }
     
    float alphaMinMaxX[2];
    float alphaMinMaxY[2];
    float alphaMinMaxZ[2];
    float alphaMinMax[2];
    
    float coordPlanesX[2];// ind 0: first plane, ind 1: last plane
    float coordPlanesY[2];
    float coordPlanesZ[2];
    
    int nX = stpPwr_grid->nbCols + 1;
    int nY = stpPwr_grid->nbRows + 1;
    int nZ = stpPwr_grid->nbFrames + 1;
    
    
    
    float xSpacing = (stpPwr_grid->inc).x;
    float ySpacing = (stpPwr_grid->inc).y;
    float zSpacing = (stpPwr_grid->inc).z;
    
    int iMinMax[2];
    int jMinMax[2];
    int kMinMax[2];
    
    float startX = (stpPwr_grid->start).x;
    float startY = (stpPwr_grid->start).y;
    float startZ = (stpPwr_grid->start).z;
    
    
    float d12 =  euclideanDistance( point1,point2);
    
 
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

    
    //Step 5 Calculate parametric sets
    float * listAlphaX, * listAlphaY, *listAlphaZ;
    int lenListX, lenListY, lenListZ;
    lenListX = iMinMax[1]-iMinMax[0]+1;
    lenListY = jMinMax[1]-jMinMax[0]+1;
    lenListZ = kMinMax[1]-kMinMax[0]+1;

    listAlphaX = (float *)malloc(sizeof(float) * lenListX);
    if(!listAlphaX){
        perror("Malloc alpha X error... \n");
        return FAILURE;
    }
    
    listAlphaY = (float *)malloc(sizeof(float) * lenListY);
    if(!listAlphaY){
        perror("Malloc alpha Y error... \n");
       return FAILURE;
    }
    listAlphaZ = (float *)malloc(sizeof(float) * lenListZ);
    if(!listAlphaZ){
        perror("Malloc alpha Z error... \n");
        return FAILURE;
    }
    
    if( absolute(calculateParametricSets(iMinMax ,point1.x, point2.x, coordPlanesX[0], xSpacing, listAlphaX,alphaMinMax ) - sameX) >  1e-5){
        perror("Error SameX should be equal to value returned by calculateParamtricSet\n");
        return FAILURE;
    }
    if( absolute(calculateParametricSets(jMinMax ,point1.y, point2.y, coordPlanesY[0], ySpacing, listAlphaY,alphaMinMax) - sameY) > 1e-5){
        perror("Error SameY should be equal to value returned by calculateParamtricSet\n");
        return FAILURE;
    }
    if( absolute(calculateParametricSets(kMinMax ,point1.z, point2.z, coordPlanesZ[0], zSpacing, listAlphaZ,alphaMinMax) - sameZ) > 1e-5){
        perror("Error SameZ should be equal to value returned by calculateParamtricSet\n");
       return FAILURE;
    }


    // Step 6 merge alpha sets:
    int lenAlpha = ( 2 + lenListX + lenListY + lenListZ);
    float * listAlpha = (float *) malloc(sizeof(float) * lenAlpha);
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

    
    
    // Step 8 Calculate radiological path
    if( -1 == calculateRadiologicalPath(stpPwr_grid, radiological_depth_grid, indList,  listAlpha,lenAlpha, d12,tabCount)){
        perror("Error while calculating the radiological depth.\n");
        exit(1);
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


float coordPlane(int idx, float coordPlane1, float spacing){
    /*
     Function that computes and return the coordinate of a plane for a given index (idx)
     For instance X_plane(i) = X_plane(1) + (i-1)dx
     Y_plane(j) = Y_plane(1) + (j-1)dy
    Z_plane(k) = Z_plane(1) + (k-1)dz
    
    see eq(3) in Siddon 1985
    
    idx: given index
    coordPlane1 : coordinate of the first plane
    spacing : spacing between planes
    */
   return (coordPlane1 + (idx-1.0) * spacing);
    
    
}

int computeAlpha(int idx, float coord1, float coord2,float coordPlane1, float spacing , float * ptAlpha){
    /*
     Function that computes alpha value for a given index (idx).
     For instance:
     if X2 - X1 != 0:
     alpha_x(1) = [ X_plane(1) - X1] / [X2 - X1]
     else:
         alpha_x(1) is undefined
    see eq(4) in Siddon 1985
    idx : given index
    coord1 : coordinate of the first point (x or y or z value)
    coord2 : coordinate of the second point (x or y or z value) 
    coordPlane1 : coordinate of the first plane
    spacing: spacing between planes for which we compute alpha values
    ptAlpha : pointer where we will store computed alpha 
    Note: for undefined alpha values we will set alpha to FLT_MAX, and return -1,
    return 0 if alpha is not undefined
    
    Computes alpha value for intersection between ray and the idx_th plan in the given direction (X, Y or Z). The direction is given by the function call in the calling function.
    
    */
    if(!ptAlpha){
        perror("ptAlpha is Null is compute alpha\n");
        exit(1);
    }
    
    float deltaCoord = coord2 - coord1;

    if (absolute(deltaCoord) <= 1e-5){
        (*ptAlpha) = FLT_MAX;
        return -1;
    }
    float coordPl = coordPlane(idx,coordPlane1,spacing);
    (*ptAlpha) = (coordPl - coord1) / deltaCoord; 
    return 0;
}


void sortAlphaMinMax(float * alphasMinMax){
    /*
    Function that makes sure that alpha min is the smallest value et alpha max the highest
    */
    
    float min = alphasMinMax[0];
    float max = alphasMinMax[1];
    if (max < min){
        alphasMinMax[1] = min;
        alphasMinMax[0] = max;
    }
}



void findAlphaMinMax(float * alphasX, float * alphasY, float * alphasZ, float * alphaMinMax){
    /*
     Function that finds alpha_min and alpha_max from alpha values. See eq(5) in Siddon 1985.
     alpha_min = max{ 0 , min[ alpha_x(1) , alpha_x(Nx)],min[ alpha_y(1) , alpha_y(Ny)],min[ alpha_z(1) , alpha_z(Nz)]   }
     
     alpha_max = min{ 1 , max[ alpha_x(1) , alpha_x(Nx)],max[ alpha_y(1) , alpha_y(Ny)],max[ alpha_z(1) , alpha_z(Nz)]   }
    alphasX : [alpha_x(1) , alpha_x(Nx) ] : array of alpha values fo x coord
    alphasY : [ alpha_y(1) , alpha_y(Ny)] : array of alpha values fo y coord
    alphasZ : [ alpha_z(1) , alpha_z(Nz)] : array of alpha values fo z coord
    alphaMinMax : [ alpha_min , alpha_max] : array to fill will alpha min and max
    */
    if(!alphasX || !alphasY || !alphasZ || !alphaMinMax){
        perror("alphasX || alphasY || alphasZ || alphaMinMax  null in findAlphaMinMax\n");
        exit(1);
    }
    float minX,maxX,minY,maxY,minZ,maxZ;
    minX = alphasX[0] < 0 ? 0 : alphasX[0];
    maxX = alphasX[1] < 0 ? 0 : alphasX[1];
    minY = alphasY[0] < 0 ? 0 : alphasY[0];
    maxY = alphasY[1] < 0 ? 0 : alphasY[1];
    minZ = alphasZ[0] < 0 ? 0 : alphasZ[0];
    maxZ = alphasZ[1] < 0 ? 0 : alphasZ[1];
    if(minX > maxX){
        perror("minX > maxX in find alpha min max\n");
        exit(1);
    }
    if(minY > maxY){
        perror("minY > maxY in find alpha min max\n");
        exit(1);
    }
    if(minZ > maxZ){
        perror("minZ > maxZ in find alpha min max\n");
        exit(1);
    }
    
    float tabMin[4];
    float tabMax[4];
    tabMin[0] = 0;
    tabMin[1] = minX == FLT_MAX ? 0 : minX;
    tabMin[2] = minY == FLT_MAX ? 0 : minY;
    tabMin[3] = minZ == FLT_MAX ? 0 : minZ;
    
    tabMax[0] = 1;
    tabMax[1] = maxX == FLT_MAX ? 1 : maxX;
    tabMax[2] = maxY == FLT_MAX ? 1 : maxY;
    tabMax[3] = maxZ == FLT_MAX ? 1 : maxZ;
    
    alphaMinMax[0] = tabMin[0];
    alphaMinMax[1] = tabMax[0];
    for(int i = 1 ; i < 4 ; i++){
        if(tabMin[i] > alphaMinMax[0]) alphaMinMax[0] = tabMin[i];
        if(tabMax[i] < alphaMinMax[1]) alphaMinMax[1] = tabMax[i];
    }
    
    if ( alphaMinMax[0] < 0) alphaMinMax[0] = 0;
    if ( alphaMinMax[1] > 1) alphaMinMax[1] = 1;
    
    
}


void findIndMinMax(float coordPlane1,float coordPlaneN,float spacing,float coord1, float coord2, int nbPlanes, float * alphaMinMax , int * indMinMax){
    /* Finds minimum index of intersected planes for given coordinate planes.
     see eq(6) in Siddon 1985
     if (X2 - X1) >= 0:
        i_min = N_x - [Xplane(N_x) - alpha_min (X2 - X1) - X1] / dx 
        i_max = 1 + [X1 + alpha_max (X2 - X1) - Xplane(1)] / dx
    if (X2 - X1) <= 0:
        i_min = N_x - [Xplane(N_x) - alpha_max (X2 - X1) - X1] / dx 
        i_max = 1 + [X1 + alpha_min (X2 - X1) - Xplane(1)] / dx
        
        coordPlane1: coord first plane
        coordPlaneN: coord last plane
        spacing : plane spacing
        coord1 : coord P1
        coord2 : coord P2
        nbPlanes : number of planes
        alphaMinMax : overall alpha min and max
        indMinMax : table[2] to fill with index min index max
     */
    if( ! alphaMinMax || ! indMinMax) {
        perror("AlphaMinMax or IndMinMax Null in findIndMinMax\n");
        exit(1);
    } 
    float tmpVal;
    if(coord2 - coord1 >= 0){
        tmpVal = float(nbPlanes) - (coordPlaneN - alphaMinMax[0] * (coord2 - coord1) - coord1) / spacing ;
        indMinMax[0] = roundFloat(tmpVal);
        tmpVal =  1.0 + (coord1 + alphaMinMax[1] * (coord2 - coord1) - coordPlane1) / spacing;
        indMinMax[1]  = roundFloat(tmpVal);
    }
    else{
        tmpVal = float(nbPlanes) - ((coordPlaneN - alphaMinMax[1] * (coord2 - coord1) - coord1) / spacing);
        indMinMax[0] = roundFloat(tmpVal);
        tmpVal = 1.0 + ((coord1 + alphaMinMax[0] * (coord2 - coord1) - coordPlane1) / spacing);
        indMinMax[1]  = roundFloat(tmpVal);
    }
    if(indMinMax[0] <= 0 || indMinMax[1] <= 0 ||indMinMax[0] >  nbPlanes || indMinMax[1] >  nbPlanes ){ // indices start at 1 and not 0!!
        printf("Error indices min (%i) max(%i)  , nbPlanes :%i , alphas min(%f) max(%f), coord1-2 (%f,%f), coordPlane1: %f, coorPlaneN:%f, spacing:%f \n",indMinMax[0],indMinMax[1],nbPlanes, alphaMinMax[0], alphaMinMax[1],coord1,coord2,coordPlane1,coordPlaneN,spacing);
        exit(1);
    }
    if( 0  > (indMinMax[1] - indMinMax[0] + 1 )){
        printf("0  > (indMinMax[1] - indMinMax[0] + 1 ) : %i / %i\n",indMinMax[0],indMinMax[1]);
        perror("0  > (indMinMax[1] - indMinMax[0] + 1 ) \n");
        exit(1);
    } 
}



int calculateParametricSets(int * indMinMax , float coord1, float coord2, float coordPlane1, float spacing, float * alphaSet, float *alphaMinMax){
    /*
     Calculate parametric sets (intersection points of the ray and parallel planes). 
     As explained in eq(8) in Siddon 1985:
     alpah_x = {alpha_x(imin) , ... , alpha_x(imax)}  if (X2-X1)> 0
     alpah_x = {alpha_x(imax) , ... , alpha_x(imin)}  if (X2-X1)< 0
     alpha_x(i)=[Xplane(i) - X1]/(X2 - X1)
                = alpha_x(i-1) + [dx / (X2-X1)]
    similar for Y and Z
    
    Note the function will return 0 if (X2-X1) != 0, -1 otherwise
    indMinMax : index min , index max
    coord1 : coord P1
    coord2 : coord P2
    coordPlane1 :coord first plane
    spacing : plane spacing
    alphaSet : set of alpha values to be filled: it should be allocated before the function is called. Its length should be : abs(indMinMax[1]-indMinMax[0])+1  
     
     */
    if(!indMinMax || !alphaSet){
        perror("IndMinMax or AlphaSet Null in calculate parametric sets.\n");
        exit(1);
    }
    int nbAlphas = indMinMax[1]-indMinMax[0]+1;
    if (absolute(coord2 - coord1) <= 1e-5 ){
        printf("abs(coord2 - coord1) <= 1e-5 - calculate parametricSet\n ");
      for( int i = 0 ; i < nbAlphas ; i++)  
            alphaSet[i] = FLT_MAX;
        
      return -1;  
    } 
    
    
    if((coord2 - coord1) > 0){
        for(int i = 0 ; i < nbAlphas ; i ++){
        
            alphaSet[i] = (coordPlane(indMinMax[0]+i, coordPlane1,spacing) - coord1) / (coord2 - coord1);
            if(alphaSet[i] < alphaMinMax [0]) alphaSet[i] = alphaMinMax [0];
            else{
                if(alphaSet[i] > alphaMinMax [1]) alphaSet[i] = alphaMinMax [1];
            }
        }
    }
    else{
        for(int i= 0 ; i< nbAlphas ; i ++){
                
                alphaSet[i] =  (coordPlane(indMinMax[1]-i, coordPlane1,spacing) - coord1) / (coord2 - coord1);
                if(alphaSet[i] < alphaMinMax [0]) alphaSet[i] = alphaMinMax [0];
                else{
                    if(alphaSet[i] > alphaMinMax [1]) alphaSet[i] = alphaMinMax [1];
                }
        }
    }
    return 0;
}




void mergeSets(float * alphaMinMax, float * alphaX, int lenAlphaX, float * alphaY, int lenAlphaY , float * alphaZ, int lenAlphaZ, float * alpha, float d12 , float spacing){
    /*
     Merge alpha sets, such that : alpha = { alphaMin, merge_by_ascending_order{alphaX,alphaY,alphaZ}, alphaMax} = {alpha(0), ... , alpha(n)}
     see eq(8) in Siddon 1985.
     
     alphaMinMax : alpha min and alpha max values
     alphaX : alpha X values
     lenAlphaX : number of alpha X values
     alphaY : alpha Y values
     lenAlphaY : number of alpha Y values
     alphaZ : alpha Z values
     lenAlphaZ : number of alpha Z values
     alpha: table to fill with alpha values merged: it should be allocated with (lenAlphaX+lenAlphaY+lenAlphaZ+2) float values
     d12 : distance P1 P2 : this is only used at the end when we verify that the set is in acsending order.
     spacing : grid average spacing.this is only used at the end when we verify that the set is in acsending order.
     */
    
    if (!alphaMinMax || !alphaX ||!alphaY || !alphaZ || !alpha){
        perror("alphaMinMax or alphaX or alphaY or alphaZ or alpha null in mergeSets\n");
        exit(1);
    }
    
    
    if (alphaX[0] == FLT_MAX) lenAlphaX = 0;
    if (alphaY[0] == FLT_MAX) lenAlphaY = 0;
    if (alphaZ[0] == FLT_MAX) lenAlphaZ = 0;
    
    float * tab0 = (float *) malloc(sizeof(float) * (lenAlphaX + lenAlphaY));
    if(!tab0){
        perror("Error - memory allocation tab0 in merge alphas\n");
        exit(1);
    }
    float * tab1 = (float *) malloc(sizeof(float) * (lenAlphaX + lenAlphaY+lenAlphaZ) );
    if(!tab1){
        perror("Error - memory allocation tab1 in merge alphas\n");
        exit(1);
    } 
    merge2Arrays(tab0, alphaX, lenAlphaX, alphaY, lenAlphaY);

    float * tab3 = alpha + 1;
    merge2Arrays(tab3, tab0, lenAlphaX+lenAlphaY, alphaZ, lenAlphaZ);
    
    alpha[0] = alphaMinMax[0] < alphaMinMax[1] ? alphaMinMax[0] : alphaMinMax[1] ;
    alpha[lenAlphaX+lenAlphaY+lenAlphaZ+1] = alphaMinMax[0] < alphaMinMax[1] ? alphaMinMax[1] : alphaMinMax[0];
    
    
    free(tab0);
    free(tab1);
}

void merge2Arrays(float * dest, float * array1, int lenArray1, float * array2, int lenArray2){
    /*
     Function used to merge and sort 2 arrays in increasing order.
     dest: array to be filled. It should be allocated priori to the call with a length of lenArray1 + lenArray2
     array1 : first array
     lenArray1: length of the first array
     array2: second array
     lenArray2: length of the second array
     
     */
    if( !dest || !array1 || !array2){
        perror("In merge 2 arrays: dest or array1 or array2 Null");
        exit(1);
    }
    int countDest = 0;
    int countA1 = 0;
    int countA2 = 0;
    int lenNewArray = lenArray1 + lenArray2;
    while(countDest < lenNewArray ){
        if( countA1 < lenArray1 && countA2 < lenArray2){
            if (array1[countA1] < array2[countA2]){
                dest[countDest] = array1[countA1];
                countA1++;
                countDest++;
            } 
            else{
                dest[countDest] = array2[countA2];
                countA2++;
                countDest++;            
            }
        
        }else{
            if(countA1 < lenArray1){
                dest[countDest] = array1[countA1];
                countA1++;
                countDest++;
            }
            else{
                dest[countDest] = array2[countA2];
                countA2++;
                countDest++;   
            }
        }
    }

}



void computeVoxelsLengths(float * alpha, int lenAlpha, float * lengthsVoxels, float d12){ // can be skipped and done directly when computing the radiological path
    /*
     2 adjacent terms in alpha list correspond to the intersection of a particular voxel and the ray. 
     There are  lenAlpha-1 voxels intersected. For m = 1 , ... ,lenAlpha, 
     l(m) = d12 * [ alpha(m) - alpha(m-1) ] , where d12 is the euclidean distance between P1 and P2.
     See eq(10 and 11) in Siddon 1985
     
     
     alpha: alpha values
     lenAlpha: number of alpha values 
     lengthsVoxels: table where lengths will be stored. It should be allocated with a (lenAlpha-1) float values
     d12: distance between P1 and P2
     */
    if(!alpha || !lengthsVoxels){
        perror("Alpha of LengthsVoxels Null in compute Voxel lengths\n");
        exit(1);
    }
    for(int i=1 ; i < lenAlpha; i ++){
        lengthsVoxels[i-1] = (alpha[i] - alpha[i-1]) * d12;
    }
}

float euclideanDistance( POINTXYZ p1, POINTXYZ p2){
    /* Computes  and returns the euclidean distance between 2 points P1 and P2 */
    return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z)  );
    
    
}

void calculateVoxelIndices(POINTXYZ p1, POINTXYZ p2, POINTXYZ plane1, POINTXYZ spacing, float * alpha, int lenAlpha, int ** indList, int ** indMinMax){
/* 
      Calculates voxel indices (for 1 coordinate, i.e. x or y or z depending on input parameters).
      i(m) = 1 + [X1 + alpha_mid(X2-X1) - Xplane(1)] / dx
      alpha_mid = ( alpha(m) + alpha(m-1)) / 2              m in [1 , lenAlpha-1]
      see eq(13 and 12) in Siddon 1985
      
      p1: coordinates P1
      p2: coordinates P2
      plane1: coordinates first planes
      spacing: spacing between planes
      alpha: list of alpha values
      lenAlpha: length of alpha list
      indList: list of indices
      indMinMax: list of indMin and indMax
    */ 
    
    float p1x = p1.x;
    float p1y = p1.y;
    float p1z = p1.z;

    float p2x = p2.x;
    float p2y = p2.y;
    float p2z = p2.z;
    
    float plane1x = plane1.x;
    float plane1y = plane1.y;
    float plane1z = plane1.z;
    
    float spacingx = spacing.x;
    float spacingy = spacing.y;
    float spacingz = spacing.z;
    

    float alphaMid, tmpVal;
    for( int m= 1; m < lenAlpha; m++){

            alphaMid = 0.5*(alpha[m]+alpha[m-1]);
            tmpVal = 1.0 + (( p1x + alphaMid * ( p2x - p1x) - plane1x) / spacingx) ;
            indList[m-1][0] = (int)floor(tmpVal);

            tmpVal =  1.0 + (( p1y + alphaMid * ( p2y - p1y) - plane1y) / spacingy);
            indList[m-1][1] = (int)floor(tmpVal);

            tmpVal = 1.0 + (( p1z + alphaMid * ( p2z - p1z) - plane1z) / spacingz); 
            indList[m-1][2] = (int)floor(tmpVal);

    }
   
}

int roundFloat( float value){

    float tmp = floor((float)value);
    int ret;
    if ( ((float)value - tmp) > 0.5)
		ret = (int) tmp + 1 ; 
	else
		ret = (int) tmp ;
    return ret;

}


float calculateRadiologicalPath(VOLUME * densGrid, VOLUME * radGrid, int ** listInd, float * listAlpha,int  lenAlpha,float d12, int * countPtr){
    /*
     Function that calculates the radiological path along the ray. See eq(14) in Siddon 1985
 
    densGrid: 3D grid of density
    radGrid: 3D grid of radiological depth
    indList : list of indices 2D array
    alpha: alpha values
    lenAlpha: length of alpha
    d12: distance P1 - P2
     */
    if( !densGrid || !radGrid || !listInd || !listAlpha ){
        perror("densGrid || radGrid || indList || alpha  null in calculateRadiologicalPath\n");
        exit(1);
    }
    
    int i,j,k;
    int nbRows = (densGrid)->nbRows;
    int nbCols = (densGrid)->nbCols;
    int nbFrames  = (densGrid)->nbFrames;
    float radiolDepth = 0;
    float voxelDensity;
    int index;

    for ( int m= 1; m < lenAlpha; m++){
            if ((listAlpha[m] - listAlpha[m-1]) > 0){
            if (listInd[m-1][0] > 0 && listInd[m-1][1] > 0 && listInd[m-1][2] > 0 ){
                i = listInd[m-1][0]-1;
                j = listInd[m-1][1]-1;
                k = listInd[m-1][2]-1;
                if ( i < nbCols && j < nbRows && k < nbFrames){ 
                    index = i + nbCols * j + k * nbRows * nbCols;
                    voxelDensity = densGrid->matrix[index];
                    if(voxelDensity < 0){
                        printf("Voxel density(%i,%i,%i): %f\n",i,j,k,voxelDensity);
                        printf("In cal radiol depth : density < 0!\n");
                        return -1;
                    }
        
                    if ( radGrid->matrix[index] < 0){
                        radGrid->matrix[index] = radiolDepth ;
                    }
                    else{
                        radGrid->matrix[index] += radiolDepth ;
                     }
                    
                    radiolDepth = radiolDepth + d12 * voxelDensity * (listAlpha[m] - listAlpha[m-1]);
                     countPtr[index] += 1;
                }
            }
       }
    }
    
    return radiolDepth;
}


float absolute(float val){
    if (val < 0) return (-val);
    return val;

}


