########################################################################
#
# Makefile
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# June 2013
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
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# do something
INCLUDES = -I/usr/local/lib/python2.7/lib/python2.7/site-packages/numpy/core/include/numpy/ -I/usr/local/lib/python2.7/include/python2.7/

FLAGS_LINK = -shared 

FLAGS_COMP = -fpic
endif


ifeq ($(UNAME), Darwin)
# do something
INCLUDES = -I/Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/numpy/core/include/numpy/ -I/Library/Frameworks/Python.framework/Versions/Current/include/python2.7/

FLAGS_LINK = -bundle -flat_namespace -undefined suppress 

FLAGS_COMP =  -O3 -Wall

endif


CPP = g++

# ---- Link --------------------------- 
SPOTSETDOUBLE: _fillDoseArray_PythonExtDouble.so _fillDensityArray_PythonExtDouble.so _fillCTArray_PythonExtDouble.so _fillStopPwrArray_PythonExtDouble.so


_fillDoseArray_PythonExtDouble.so: fillDoseArray_PythonExtDouble.o     
			$(CPP) $(FLAGS_LINK) -o _fillDoseArray_PythonExtDouble.so fillDoseArray_PythonExtDouble.o

_fillDensityArray_PythonExtDouble.so: fillDensityArray_PythonExtDouble.o 
			$(CPP) $(FLAGS_LINK) -o _fillDensityArray_PythonExtDouble.so fillDensityArray_PythonExtDouble.o 

_fillCTArray_PythonExtDouble.so: fillCTArray_PythonExtDouble.o   
			$(CPP) $(FLAGS_LINK) -o _fillCTArray_PythonExtDouble.so fillCTArray_PythonExtDouble.o

_fillStopPwrArray_PythonExtDouble.so: fillStopPwrArray_PythonExtDouble.o 
			$(CPP) $(FLAGS_LINK) -o _fillStopPwrArray_PythonExtDouble.so fillStopPwrArray_PythonExtDouble.o 


# ---- gcc C compile ------------------

fillDoseArray_PythonExtDouble.o: fillDoseArray_PythonExtDouble.cpp fillDoseArray_PythonExtDouble.h  
			$(CPP) $(FLAGS_COMP) -c fillDoseArray_PythonExtDouble.cpp $(INCLUDES)

fillDensityArray_PythonExtDouble.o: fillDensityArray_PythonExtDouble.cpp fillDensityArray_PythonExtDouble.h 
			$(CPP) $(FLAGS_COMP) -c fillDensityArray_PythonExtDouble.cpp $(INCLUDES)

fillCTArray_PythonExtDouble.o: fillCTArray_PythonExtDouble.cpp fillCTArray_PythonExtDouble.h
			$(CPP) $(FLAGS_COMP) -c fillCTArray_PythonExtDouble.cpp $(INCLUDES)

fillStopPwrArray_PythonExtDouble.o: fillStopPwrArray_PythonExtDouble.cpp fillStopPwrArray_PythonExtDouble.h     
			$(CPP) $(FLAGS_COMP) -c fillStopPwrArray_PythonExtDouble.cpp $(INCLUDES)
 

clean:
	rm -f *.o *~
	echo clean done...

cleanSO:
	rm -f *.so
	echo .SO clean done... 