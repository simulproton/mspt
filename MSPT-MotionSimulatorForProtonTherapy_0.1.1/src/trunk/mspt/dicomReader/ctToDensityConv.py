########################################################################
#
# ctToDensityConv.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# November 2013
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
########################################################################
import numpy as np
import os

pathToFile =(os.getcwd()+'/RefData/CTToDensData/ct_to_density_MSPT.txt')



def getConversionTableCTToMassDens(name = 'MSPT'):
    '''Loads the data (text file) used to convert CT numbers to mass density. \
    Then creates a list as: [(CT Num0, mass dens0 ) , (CT Num1, mass dens1 ) ,...]
    
    :param name: Name of the file where is contained the CT to density conversion data. It must be stored in \
    the folder '/RefData/CTToDensData/'. By default it is 'MSPT'.
    
    .. note:: 
        
        The CT to density file name must follow this pattern 'ct_to_density_NAMEDATA.txt' where NAMEDATA is the name \
        chosen by the user. For example by default in MSPT, NAMEDATA is 'MSPT' the default variable provided to \
        this function. To create this file proceed as follows:
        
            #. Create a new tab delimited text file named 'ct_to_density_NAMEDATA.txt'
            #. In the first line, enter: 1st column: 'CT_Number' , 2nd column: 'Mass_Density_g.cm3'
            #. Then for each row enter: 1st column: the CT number (in increasing order) , 2nd column: the corresponding mass density.
            #. Place the file in the folder '/RefData/CTToDensData/'
            #. Set the MSPT configuration 'nameCTToDensTable' variable 
    
    :returns: The created list.
    
    '''
    
    listConversion = []
    if name is None:
        filename = pathToFile
    else:
        filename = os.getcwd()+'/RefData/CTToDensData/ct_to_density_%s.txt'%name
    
    print "\t>>>>>>>> Importing CT to density table %s <<<<<<<<"%name
    if os.path.isfile(filename):
        data = np.loadtxt(pathToFile,skiprows=1)
        for value in zip(data[:,0],data[:,1]):
            listConversion.append(value)
        print "CT to density table %s imported"%name
    else:
        raise ValueError('No file %s was found to obtaine the conversion data from CT number to density'%filename)
    print "\t>>>>>>>> End Importing CT to density table <<<<<<<<"
    return listConversion
