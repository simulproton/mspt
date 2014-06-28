########################################################################
#
# ctToDensityConv.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# November 2013
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
import numpy as np
import os

pathToFile =(os.getcwd()+'/RefData/RayStation_Generic_CT_to_density.txt')



def getConversionTableCTToMassDens():
    '''Loads the data (text file) used to convert CT numbers to mass density. \
    Then creates a list as: [(CT Num0, mass dens0 ) , (CT Num1, mass dens1 ) ,...]
    
    :returns: The created list.
    
    '''
    listConversion = []
    data = np.loadtxt(pathToFile,skiprows=1)
    for value in zip(data[:,0],data[:,1]):
        listConversion.append(value)
    return listConversion
