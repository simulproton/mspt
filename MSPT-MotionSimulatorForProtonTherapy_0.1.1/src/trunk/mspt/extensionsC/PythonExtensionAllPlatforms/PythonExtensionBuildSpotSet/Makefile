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

FLAGS_COMP = -O3 -Wall

endif


CPP = g++

# ---- Link --------------------------- 
SPOTSET: _fillDoseArray_PythonExt.so _fillDensityArray_PythonExt.so _fillCTArray_PythonExt.so _fillStopPwrArray_PythonExt.so


_fillDoseArray_PythonExt.so: fillDoseArray_PythonExt.o     
			$(CPP) $(FLAGS_LINK) -o _fillDoseArray_PythonExt.so fillDoseArray_PythonExt.o

_fillDensityArray_PythonExt.so: fillDensityArray_PythonExt.o 
			$(CPP) $(FLAGS_LINK) -o _fillDensityArray_PythonExt.so fillDensityArray_PythonExt.o 

_fillCTArray_PythonExt.so: fillCTArray_PythonExt.o   
			$(CPP) $(FLAGS_LINK) -o _fillCTArray_PythonExt.so fillCTArray_PythonExt.o

_fillStopPwrArray_PythonExt.so: fillStopPwrArray_PythonExt.o 
			$(CPP) $(FLAGS_LINK) -o _fillStopPwrArray_PythonExt.so fillStopPwrArray_PythonExt.o 


# ---- gcc C compile ------------------


fillDoseArray_PythonExt.o: fillDoseArray_PythonExt.cpp fillDoseArray_PythonExt.h  
			$(CPP) $(FLAGS_COMP) -c fillDoseArray_PythonExt.cpp $(INCLUDES)

fillDensityArray_PythonExt.o: fillDensityArray_PythonExt.cpp fillDensityArray_PythonExt.h 
			$(CPP) $(FLAGS_COMP) -c fillDensityArray_PythonExt.cpp $(INCLUDES)

fillCTArray_PythonExt.o: fillCTArray_PythonExt.cpp fillCTArray_PythonExt.h
			$(CPP) $(FLAGS_COMP) -c fillCTArray_PythonExt.cpp $(INCLUDES)

fillStopPwrArray_PythonExt.o: fillStopPwrArray_PythonExt.cpp fillStopPwrArray_PythonExt.h     
			$(CPP) $(FLAGS_COMP) -c fillStopPwrArray_PythonExt.cpp $(INCLUDES)
 

clean:
	rm -f *.o *~
	echo clean done...

cleanSO:
	rm -f *.so
	echo .SO clean done... 