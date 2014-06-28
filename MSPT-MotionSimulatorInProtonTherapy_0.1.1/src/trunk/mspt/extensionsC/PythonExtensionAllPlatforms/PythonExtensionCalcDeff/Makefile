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

FLAGS_COMP = -O3 -Wall

endif

CPP = g++


# ---- Link --------------------------- 
DEFF: _calcRadiolDepth_PythonExt.so _singleRay_PythonExt.so _fewRays_RadDepth_PythonExt.so


_calcRadiolDepth_PythonExt.so: calcRadiolDepth_PythonExt.o calcRadiolDepth.o      
			$(CPP) $(FLAGS_LINK) -o _calcRadiolDepth_PythonExt.so calcRadiolDepth.o calcRadiolDepth_PythonExt.o
			
_singleRay_PythonExt.so:  singleRay_PythonExt.o singleRay.o calcRadiolDepth.o
			g++ $(FLAGS_LINK) -o _singleRay_PythonExt.so singleRay.o singleRay_PythonExt.o calcRadiolDepth.o
			

_fewRays_RadDepth_PythonExt.so: fewRays_RadDepth_PythonExt.o calcRadiolDepth.o      
			$(CPP) $(FLAGS_LINK) -o _fewRays_RadDepth_PythonExt.so calcRadiolDepth.o fewRays_RadDepth_PythonExt.o





# ---- gcc C compile ------------------

calcRadiolDepth.o: calcRadiolDepth.cpp calcRadiolDepth.h
			$(CPP) $(FLAGS_COMP) -c calcRadiolDepth.cpp

calcRadiolDepth_PythonExt.o: calcRadiolDepth_PythonExt.cpp calcRadiolDepth.h calcRadiolDepth_PythonExt.h    
			$(CPP) $(FLAGS_COMP)  -c calcRadiolDepth_PythonExt.cpp $(INCLUDES)

singleRay_PythonExt.o: singleRay_PythonExt.cpp singleRay_PythonExt.h singleRay.h calcRadiolDepth.h
			g++  $(FLAGS_COMP) -c singleRay_PythonExt.cpp $(INCLUDES)

singleRay.o: singleRay.cpp singleRay.h calcRadiolDepth.h  
			g++  $(FLAGS_COMP) -c singleRay.cpp $(INCLUDES)

fewRays_RadDepth_PythonExt.o: fewRays_RadDepth_PythonExt.cpp calcRadiolDepth.h fewRays_RadDepth_PythonExt.h    
			$(CPP) $(FLAGS_COMP)  -c fewRays_RadDepth_PythonExt.cpp $(INCLUDES)


clean:
	rm -f *.o *~
	echo clean done... 
	
cleanSO:
	rm -f *.so
	echo .SO clean done... 