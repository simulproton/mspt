'''Package containing python extensions in C wrapped in C++ functions:

    * PythonExtensionBeamletCalc: Extension **_beamletSimulation_PythonExt** to compute the dose distribution of a single pencil beam using float values.
        
        
       **static PyObject * beamletSimul(PyObject *self, PyObject *args)** : API Python - C
    
        BeamletSimul simulates the effect of a single pencil beam on a 3D volume.
        
        The simulation is based on Hong et al. paper, "Pencil beam algorithm for proton dose calculations",  Phys. Med. and Biol. 41,1996, and on Szymanowski et al. paper, "Two-dimensional pencil beam scaling: an improved proton dose algorithm for heterogeneous media", Phys. Med. and Biol. 47,2002.
            
            * Inputs:
            
                #. rhomw : density ratio to water
                #. Smw :stopping power ratio to water
                #. Xr : 3D matrix representing x coordinates in the IEC Gantry coordinate system (cm).
                #. Yr : 3D matrix representing y coordinates in the IEC Gantry coordinate system (cm).
                #. Zr : 3D matrix representing z coordinates in the IEC Gantry coordinate system (cm).
                #. dmax : maximal depth in water reached for the beam energy (cm).
                #. z_dd : "depth" coordinates for the depth dose curve (cm).
                #. dd : "dose" coordinates for the depth dose curve. (Gy.cm ^2 )
                #. sig0 : beam size at patient entrance (cm).
                #. xrshift, yrshift : x and y beam shifts (from isocenter) used to define target position. Provided in the IEC gantry coordinate system.(cm)
                #. sad : source-axis distance: distance between the source and the isocenter.(cm)
                #. Deff : 3D matrix representing the radiological depth.(cm)
                #. ssd0 : source surface distance
                #. indicesRay : indices of the voxels along the ray
                #. x0y0z0IECF: beam target point in the IEC fixed coordinate system.
                #. spacingIECF : voxel spacing in the IEC fixed coord system
                #. invRotMat : Matrix to go from the IECG to IECF coordinate system
      
            * Output: 
                3D numpy array dose distribution.
    
       **float radialEmmittance( float deff, float dmax, float sig0)**
            
            Computes the radial emittance (beam's width): 
     
                * first part depth dependent: (beam spread with depth) due to multiple scattering inside patient body.
                
                The beam spread parameter is calculated using these values as:
                    
                    y_t = y_R*(Afit*(t/R)^2 + Bfit*(t/R))
    
                where y_R = Cfit*R^2 + Dfit*R, R = dmax (Bragg peak depth)
    
                This is from the paper by Hong (PMB 1996) - Figure A4.
    
                * second part: beam spread at patient entrance : sig0.
                
                Finally:
            
                    the total radial emittance. See eq (8) in Hong 1996.
        
                        sigTot = sqrt( sig0^2 + sqrt(y_t)^2)
    
            The function returns sigTot^2
           
           
                * Inputs:
                
                    * deff : radiological depth 
                    * dmax : beam range in water
                    * sig0 : beam's spread at the patient entrance
                    
                * Output:
                
                    sigTot^2
                     
           **float centralAxisTerm(float dd, float ssd0, float deff, float z)**:
           
                Compute the central axis term, from Hong 1996, eq (10)
                
                C(z') = DD(deff) * ((ssd0 + deff) / z')^2
                
                zBeV: z coordinate in the beam eye view coordinate system (BEV)
            
           **float offAxisTerm(float variance, float distPQ)**:
            
                Compute the off-axis term, from Hong 1996 eq (12)
            
                O_xyz = [ 1/(2*pi*sigTot^2) ] * exp( - ((x' - xbev)^2 + (y' - ybev)^2 ) / (2*sigTot^2) )
            
           **float offAxisTerm2Sigmas(float varX,float varY, float xrshift,float yrshift, float xBev, float yBev)**
            
                Compute the off-axis term, from Hong 1996 eq (12) extended to 2 dimensions:
     
                O_xyz = [ (1/(4*pi*sigTotX^2)) + 1/(4*pi*sigTotY^2) ] * exp( - ( ((x' - xbev)^2 / (2*sigTotX^2))  + ((y' - ybev)^2 ) / (2*sigTotY^2)) )
            
            
                Note: Not implemented yet
            
           **void getRotatedPoint( float * newP , float * P, float * rotMat)**:
            
                Rotate a point coordinates based on a rotation matrix 
            
            
           **void findPointIndicesInGrid( int * ind, float * pointP1, float * spacing, float * x0y0z0,float * rotMat)**

                Convert coordinates of a rotated point to indices in a grid  

           **float interp(float xval, float * xList, float * yList, int lenList)**:
            
                Piecewise linear interpolation function: interpolates a function defined by (xList,yList) \
                and returns 'yval' corresponding to 'xval' for the interpolated function.
            
           **float dotProduct(float * p1,float * p2)**:
            
                Computes the dot product between two vectors
            
           **float distance( float * p1, float *p2)**:
            
                Euclidean distance



    * PythonExtensionBeamletCalcDouble: Extension **_beamletSimulation_PythonExtDouble** to compute the dose distribution of a single pencil beam using double values.
    
        Same as PythonExtensionBeamletCalc for double values.
    
    
    * PythonExtensionBuildSpotSet: 4 extenstions 
    
        #.  Extension **_fillCTArray_PythonExt** to fill the CT array. 
        
           **static PyObject * fillCTGrid(PyObject *self, PyObject *args)** : API Python - C
                
            Function to fill the CT grid based on a new grid shape
 
            * Inputs:
            
                #. ctGrid : 3D numpy array representing the CT
                #. newShape : grid new shape
                #. indStart : index in the input CT grid where to start taking values.

            * Output: 
            
                3D numpy array representing a new CT grid.
                    
        #.  Extension **_fillDensityArray_PythonExt** to fill the Density array. 
        
           **static PyObject * fillDensityGrid(PyObject *self, PyObject *args)** : API Python - C
            
            Function to fill the density grid based using the CT grid and a conversion table
            
            * Inputs:
                
                #. ctGrid : 3D numpy array representing the CT
                #. conversionTable:(2D Array) 2 rows, n columns: row 0 : CT values, row 1 :density values 
    
 
            * Output: 
                
                    3D numpy array representing a new Density grid.
                    
        #.  Extension **_fillDoseArray_PythonExt** to fill the expected dose array. 
        
           **static PyObject * fillDoseGrid(PyObject *self, PyObject *args)** : API Python - C
   
            FillDoseGrid aims to fill a 3D array with the planned dose at each voxel. The originalDoseGrid and newDoseGrid don't always have the same \
            size. That is why we need to associate to each coordinate in newDoseGrid with the coordinates in originalDoseGrid and set the voxel value value in newDoseGrid with the voxel value in originalDoseGrid.\
            A resampling of the dose grid can be performed of the internal variable 'interpol' is set to 1 in the code.
            
            The goal of the function is to find in originalDoseGrid the position of each spot in newDoseGrid in order \
            to set the voxel of newDoseGrid to the corresponding dose value.
 
        
        
            * Inputs:
            
                #. dimNewGrid:Dimensions of the 3D array to fill with the dose contained in originalDoseGrid. Its size is m x n x q
                #. originalDoseGrid: 3D array containing the planned dose extracted from RD dicom file. Its size can be different than m x n x q, it can be smaller or greater.  
                #. imagePosition: coordinates (in mm) of the top left corner on the first slice of originalDoseGrid
                #. imageEndPosition: coordinates (in mm) of the bottom right corner on the last slice of originalDoseGrid
                #. imageOrientation: axis information in orginialDoseGrid : 1D array with 9 elements
                #. spacingInfoCT : spacing information (mm) in newDoseGrid: x:spacingInfoCT[0] ,y:spacingInfoCT[1] ,z:spacingInfoCT[2] ,
                #. spacingInfoRD : spacing information (mm) in originalDoseGrid (from RD File) : x:spacingInfoRD[0] ,y:spacingInfoRD[1] ,z:spacingInfoRD[2] 
                #. imagePositionRD : x,y,z coordinates of the first top left voxel of the dose data.
        
            * Output:
            
                3D numpy array representing a new Dose grid.
                
           **float trilinearInterpolation(float * valNeigh ,  float  pointNeigh[8][3] , float * pointRD)**
            
                Computes the tri-linear interpolation of the dose value, at pointRD,  between 8 neighbors.
            
           **float computeBarycenter(float valN1, float valN2, float coordP, float coordN1, float coordN2)**
            
                Computes the barycenter (dose value) to P between the dose values at N1 and N2.
            
           **void createTransformationMatrix(float * imagePosition, float * imageOrientation, float * pixelSpacing,float * transformMatrix)**
            
                Returns the transformation matrix used to go from one coordinates in mm to indices in a Dicom Array
                
                    * Image Position: refers to coordinates (x0,y0,z0) of the top left corner in the image 
                    * Image Orientation: unit vectors in the 3 directions with respect to the patient coordinate system
                    * Pixel Spacing: pixel spacing in the 3 directions
                    * transformMatrix: tha matrix to fill
     
                Formula: look at page 410 from `Dicom Documentation <http://medical.nema.org/Dicom/2011/11_03pu.pdf>`_  - conversion between dicom voxel indices and coordinates. 
            
           **void createTransformationMatrixIndToCoord(float * imagePosition, float * imageOrientation, float * pixelSpacing,float * transformMatrix)**
            
                Returns the transformation matrix used to go from indices to coordinates in mm in a Dicom Array
            
           **void indexToCoordFromImagePosition( float * point, float * matrix,int *indices)**

                Returns (x, y, z) coordinates (pixel location)for given indices (frame,row,column) provided in an image  \
                coordinate system (mm).
            
            
           **void  coordinatesToIndexFromImagePosition( float * point, float * invMatrix,int *indices)**
            
                Returns (frame, row, column) indices (pixel location)for a given point (x,y,z) provided in an image  \
                coordinate system (mm).

            
           **bool invertMatrix(float m[16], float invOut[16])**
            
                Functions that inverts a 4x4 matrix
            
            
           **float absolute(float val)**
                
                returns the absolute value.
                
                
        #. Extension **_fillStopPwrArray_PythonExt** to fill the relative stopping power array. 
        
           **static PyObject * fillStopPwrGrid(PyObject *self, PyObject *args)** : API Python - C
            
            Function that converts the density grid into relative stopping power grid, using a conversion table
 
            * Inputs:
                
                #. densityGrid : 3D numpy array in which each voxel is represented by its density
                #. conversionTable : table used for the conversion. Array of 3 rows and n columns:
                
                    * Row 0 : lowest bound of the density values for each media, 
                    * Row 1 : upper bound of the density values for each madia , 
                    * Row 3: corresponding relative mass Stopping Power value. 
                    * Each column corresponds to 1 type of medium: 'air','adipose,','muscle','bone' ...\
                        sorted by increasing density
 
            * Output: 
                
                3D numpy array in which the relative stopping power is stored.
                
           **float scalingFactor(float density)**
            
                Returns the stopping power scaling factor. It is calculated from a static table. 
  
    
    * PythonExtensionBuildSpotSetDouble
    
        Same as PythonExtensionBuildSpotSet for double values.
        
    * PythonExtensionCalcDeff: 3 extenstions and 2 utility files
    
        #.  Extension **_calcRadiolDepth_PythonExt** to fill the  entire radiological depth array. 
        
           **static PyObject * calcRadDeff(PyObject *self, PyObject *args)** : API Python - C
            
                This function is the python extension. It calculates the radiological depth of all the voxels in a volume.
                
                * Inputs:
            
                    #. densityGrid : 3D numpy array representing the volume: density or relative stopping power (depending if this is photons or protons)
                    #. startVec : numpy array (1,3) : Coordinates of the first voxel (0,0,0)
                    #. incVec : numpy array (1,3): Increment vector (spacing) when incrementing (+1,+1,+1) the voxels indices.
                    #. sourceVec : numpy array (1,3): Coordinates of the source of the ray

                * Output:
                
                    Numpy array (3D) in which voxels are represented by their radiological depth

            
        #. Extension **_fewRays_RadDepth_PythonExt** to fill the radiological depth array only along few rays. Therefore \
        only some voxels will be filled.
        
           **static PyObject * calcRadDeff(PyObject *self, PyObject *args)** : API Python - C

                This function is the python extension. It finds the voxels of a volume encountered by a ray and their radiological depth \
                for few target spots:
            
                * Inputs:
                    
                    #. grid : 3D numpy array representing the volume: density or relative stopping power (depending if this is photons or protons) 
                    #. startVec : numpy array (1,3) : Coordinates of the first voxel (0,0,0)
                    #. incVec : numpy array (1,3): Increment vector (spacing) when incrementing (+1,+1,+1) the voxels indices.
                    #. sourceVec : numpy array (1,3): Coordinates of the source of the ray
                    #. targetVec : numpy array (1,3): Coordinates of the target of the ray // This should be the "central" target
                    #. auxilaryTargetVecs : numpy array (nbTargets,3): Coordinates of the targets around the "central" target
    
                * Output:
                    
                    #. Numpy array (3D) in which voxels encountered by the ray have a value of radiological depth, and -1 otherwise

        #. Extension **_singleRay_PythonExt** to fill the radiological depth array only along a single ray.
        
       **static PyObject * calcRadDeff(PyObject *self, PyObject *args)** : API Python - C

            This function is the python extension. It calculates the voxels of a volume encountered by a ray:

                * Inputs:

                    #. grid : 3D numpy array representing the volume 
                    #. startVec : numpy array (1,3) : Coordinates of the first voxel (0,0,0)
                    #. incVec : numpy array (1,3): Increment vector (spacing) when incrementing (+1,+1,+1) the voxels indices.
                    #. sourceVec : numpy array (1,3): Coordinates of the source of the ray
                    #. targetVec : numpy array (1,3): Coordinates of the target of the ray
    
                * Output:

                    Numpy array (nbVoxels,4) in which is stored a list of voxel: for each voxel (frame,row,column) 
                    indices are stored along with the length of the ray going through the voxel.



        #. Utility **calcRadiolDepth**
        
            #. **int calcRadiolDepth(VOLUME *stpPwr_grid,VOLUME *radiological_depth_grid, POINTXYZ point1,POINTXYZ point2 , int * tabCount)**

                Main function of the process
                 
                This code is an implementation of the method described by Robert L. Siddon in "Fast calculation of \
                the exact radiological path for three-dimensionnal CT array", Am. Assoc. Phys. Med. 1985.\
                It traces a ray from point1 (start) to point2 (end) and assign to all the voxel crossed their \
                radiological depth.

     
                    * Inputs:
                    
                        #. stpPwr_grid : grid (VOLUME) storing the electron density (radiol. depth. for X-Rays) or \
                        relative stopping power(for protons).It should be initialized with 0 when radiological depth has \
                        not been computed. If a value other than 0 is already present it means the user will have \
                        to divide each value by the number of rays traced.
                        #. radiological_depth_grid : grid (VOLUME) that will be filled with the radiological depths.
                        #. point1 : source of the ray
                        #. point2 : target of the ray
                        #. tabCount : count how many rays are going through each voxels in order to perform an average afterwards.

                    * Output:
                    
                        Returns Success or Failure.
                        
                    
            #. **float coordPlane(int idx, float coordPlane1, float spacing)**:
                
                Function that computes and return the coordinate of a plane for a given index (idx)
                For instance X_plane(i) = X_plane(1) + (i-1)dx
                Y_plane(j) = Y_plane(1) + (j-1)dy
                Z_plane(k) = Z_plane(1) + (k-1)dz

                see eq(3) in Siddon 1985
                
                * Inputs:
                    #. idx: given index
                    #. coordPlane1 : coordinate of the first plane
                    #. spacing : spacing between planes
           

            #. **void sortAlphaMinMax(float * alphasMinMax)**
                
                 Function that makes sure that alpha min is the smallest value et alpha max the highest
            
            #. **int computeAlpha(int idx, float coord1, float coord2,float coordPlane1, float spacing , float * ptAlpha)**
                     
                     Function that computes alpha value for a given index (idx).
                     
                     For instance:
                     if X2 - X1 != 0: \
                     alpha_x(1) = [ X_plane(1) - X1] / [X2 - X1] \
                     else: \
                         alpha_x(1) is undefined
                     
                     see eq(4) in Siddon 1985
                     
                     * Inputs:
                     
                        #. idx: given index
                        #. coord1: coordinate of the first point (x or y or z value)
                        #. coord2: coordinate of the second point (x or y or z value) 
                        #. coordPlane1: coordinate of the first plane
                        #. spacing: spacing between planes for which we compute alpha values
                        #. ptAlpha: pointer where we will store computed alpha
                    
                        Note: for undefined alpha values we will set alpha to FLT_MAX, and return -1,
                        return 0 if alpha is not undefined

                        Computes alpha value for intersection between ray and the idx_th plan in the given direction (X, Y or Z). The direction is given by the function call in the calling function.


            #. **void findAlphaMinMax(float * alphasX, float * alphasY, float * alphasZ, float * alphaMinMax)**

                Function that finds alpha_min and alpha_max from alpha values. See eq(5) in Siddon 1985.
                
                alpha_min = max{ 0 , min[ alpha_x(1) , alpha_x(Nx)],min[ alpha_y(1) , alpha_y(Ny)],min[ alpha_z(1) , alpha_z(Nz)]   }

                    * Inputs:
                        #. alpha_max = min{ 1 , max[ alpha_x(1) , alpha_x(Nx)],max[ alpha_y(1) , alpha_y(Ny)],max[ alpha_z(1) , alpha_z(Nz)]   }
                        #. alphasX : [alpha_x(1) , alpha_x(Nx) ] : array of alpha values fo x coord
                        #. alphasY : [ alpha_y(1) , alpha_y(Ny)] : array of alpha values fo y coord
                        #. alphasZ : [ alpha_z(1) , alpha_z(Nz)] : array of alpha values fo z coord
                        #. alphaMinMax : [ alpha_min , alpha_max] : array to fill will alpha min and max
            
            
            #. **void findIndMinMax(float coordPlane1,float coordPlaneN,float spacing,float coord1, float coord2, int nbPlanes, float * alphaMinMax , int * indMinMax)**

                Finds minimum index of intersected planes for given coordinate planes. see eq(6) in Siddon 1985
                
                * if (X2 - X1) >= 0:
                
                    * i_min = N_x - [Xplane(N_x) - alpha_min (X2 - X1) - X1] / dx 
                    * i_max = 1 + [X1 + alpha_max (X2 - X1) - Xplane(1)] / dx
                    
                * if (X2 - X1) <= 0:
                
                    * i_min = N_x - [Xplane(N_x) - alpha_max (X2 - X1) - X1] / dx 
                    * i_max = 1 + [X1 + alpha_min (X2 - X1) - Xplane(1)] / dx

                * Inputs:

                    #. coordPlane1: coord first plane.
                    #. coordPlaneN: coord last plane.
                    #. spacing: plane spacing.
                    #. coord1: coord P1.
                    #. coord2: coord P2.
                    #. nbPlanes: number of planes.
                    #. alphaMinMax: overall alpha min and max.
                    #. indMinMax: table[2] to fill with index min index max.



            #. **int calculateParametricSets(int * indMinMax , float coord1, float coord2, float coordPlane1, float spacing, float * alphaSet, float *alphaMinMax)**

                Calculate parametric sets (intersection points of the ray and parallel planes). 
                 
                 As explained in eq(8) in Siddon 1985:
                 
                     alpah_x = {alpha_x(imin) , ... , alpha_x(imax)}  if (X2-X1)> 0
                     alpah_x = {alpha_x(imax) , ... , alpha_x(imin)}  if (X2-X1)< 0
                     alpha_x(i)=[Xplane(i) - X1]/(X2 - X1)
                            = alpha_x(i-1) + [dx / (X2-X1)]
                similar for Y and Z

                Note the function will return 0 if (X2-X1) != 0, -1 otherwise
                
                * Inputs:
                
                    #. indMinMax : index min , index max
                    #. coord1 : coord P1
                    #. coord2 : coord P2
                    #. coordPlane1 :coord first plane
                    #. spacing : plane spacing
                    #. alphaSet : set of alpha values to be filled: it should be allocated before the function is called. Its length should be : abs(indMinMax[1]-indMinMax[0])+1  

           
           
            #. **void mergeSets(float * alphaMinMax, float * alphaX, int lenAlphaX, float * alphaY, int lenAlphaY , float * alphaZ, int lenAlphaZ, float * alpha, float d12, float spacing)**
                  Merge alpha sets, such that : alpha = { alphaMin, merge_by_ascending_order \
                  {alphaX,alphaY,alphaZ}, alphaMax} = {alpha(0), ... , alpha(n)} \
                    see eq(8) in Siddon 1985.

                * Inputs:
                    
                     #. alphaMinMax : alpha min and alpha max values
                     #. alphaX : alpha X values
                     #. lenAlphaX : number of alpha X values
                     #. alphaY : alpha Y values
                     #. lenAlphaY : number of alpha Y values
                     #. alphaZ : alpha Z values
                     #. lenAlphaZ : number of alpha Z values
                     #. alpha: table to fill with alpha values merged: it should be allocated with (lenAlphaX+lenAlphaY+lenAlphaZ+2) float values
                     #. d12 : distance P1 P2 : this is only used at the end when we verify that the set is in acsending order.
                     #. spacing : grid average spacing.this is only used at the end when we verify that the set is in acsending order.
           
            #. **void merge2Arrays(float * dest, float * array1, int lenArray1, float * array2, int lenArray2)**
            
                Function used to merge and sort 2 arrays in increasing order.

                * Inputs:
                     #. dest: array to be filled. It should be allocated priori to the call with a length of lenArray1 + lenArray2
                     #. array1 : first array
                     #. lenArray1: length of the first array
                     #. array2: second array
                     #. lenArray2: length of the second array                    
            #. **void computeVoxelsLengths(float * alpha, int lenAlpha, float * lengthsVoxels, float d12)**
                 
                 2 adjacent terms in alpha list correspond to the intersection of a particular voxel and the ray. \
                 There are  lenAlpha-1 voxels intersected. For m = 1 , ... ,lenAlpha, \
                 l(m) = d12 * [ alpha(m) - alpha(m-1) ] , where d12 is the euclidean distance between P1 and P2. \
                 See eq(10 and 11) in Siddon 1985

                * Inputs:
                
                     #. alpha: alpha values
                     #. lenAlpha: number of alpha values 
                     #. lengthsVoxels: table where lengths will be stored. It should be allocated with a (lenAlpha-1) float values
                     #. d12: distance between P1 and P2
            
            #. **float euclideanDistance( POINTXYZ p1, POINTXYZ p2)**
            
                Computes  and returns the euclidean distance between 2 points P1 and P2
            
            #. **void calculateVoxelIndices(POINTXYZ p1, POINTXYZ p2, POINTXYZ plane1, POINTXYZ spacing, float * alpha, int lenAlpha, int** indList, int** indMinMax)**
            
                Calculates voxel indices (for 1 coordinate, i.e. x or y or z depending on input parameters).\
                i(m) = 1 + [X1 + alpha_mid(X2-X1) - Xplane(1)] / dx \
                alpha_mid = ( alpha(m) + alpha(m-1)) / 2              m in [1 , lenAlpha-1] \
                see eq(13 and 12) in Siddon 1985

                * Inputs:
                
                      #. p1: coordinates P1
                      #. p2: coordinates P2
                      #. plane1: coordinates first planes
                      #. spacing: spacing between planes
                      #. alpha: list of alpha values
                      #. lenAlpha: length of alpha list
                      #. indList: list of indices
                      #. indMinMax: list of indMin and indMax
            
            
            #. **float calculateRadiologicalPath(VOLUME * densGrid, VOLUME * radGrid, int** listInd, float * listAlpha,int  lenAlpha,float d12, int * countPtr)**
            
                Function that calculates the radiological path along the ray. See eq(14) in Siddon 1985
                
                * Inputs:

                    #. densGrid: 3D grid of density
                    #. radGrid: 3D grid of radiological depth
                    #. indList : list of indices 2D array
                    #. alpha: alpha values
                    #. lenAlpha: length of alpha
                    #. d12: distance P1 - P2
            
            #. Structures:
                
                * Struct. POINTXYZ:
                
                    * float x
                    * float y 
                    * float z 

                * Struct. VOLUME:
                
                    * POINTXYZ start 
                    * POINTXYZ inc 
                    * int nbCols 
                    * int nbRows 
                    * int nbFrames 
                    * float *matrix
                    



        #. Utility **singleRay**
        
            **int singleRaytrace(VOLUME * grid, POINTXYZ point1,POINTXYZ point2 , INDICES * indRay)**
        
                * Inputs:
                
                    #. stpPwr_grid : grid (VOLUME) storing the electron density (radiol. depth. for X-Rays) or relative stopping power(for protons)
                    #. radiological_depth_grid : grid (VOLUME) that will be filled with the radiological depths.
                    #. point1 : source of the ray
                    #. point2 : target of the ray
                    #. indRay : list of indices encountered by the ray
    
    * PythonExtensionCalcDeffDouble: 3 extenstions and 2 utility files
    
        Same as PythonExtensionCalcDeff for double  values

'''


