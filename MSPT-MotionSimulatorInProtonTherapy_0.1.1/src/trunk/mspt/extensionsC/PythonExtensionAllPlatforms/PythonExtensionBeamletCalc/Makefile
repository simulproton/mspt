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

#ifeq (,$(findstring CYGWIN_NT,$(UNAME))
## do something
#INCLUDES = -I/cygdrive/c/python27/lib/site-packages/numpy/core/include/numpy/  -I/cygdrive/c/python27/include/
#
#FLAGS_LINK = -shared 

#FLAGS_COMP = 
#endif


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
BEAMLET: _beamletSimulation_PythonExt.so


_beamletSimulation_PythonExt.so: beamletSimulation_PythonExt.o     
			$(CPP) $(FLAGS_LINK) -o _beamletSimulation_PythonExt.so beamletSimulation_PythonExt.o


# ---- gcc C compile ------------------

beamletSimulation_PythonExt.o: beamletSimulation_PythonExt.cpp beamletSimulation_PythonExt.h    
			$(CPP) $(FLAGS_COMP) -c beamletSimulation_PythonExt.cpp $(INCLUDES)

 

clean:
	rm -f *.o *~
	echo clean done...
	
cleanSO:
	rm -f *.so
	echo .SO clean done... 