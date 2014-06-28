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

# Main math makefile

export CPLUS_INCLUDE_PATH=/usr/local/lib/python2.7/site-packages/numpy/core/include/numpy:/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/:/include/python2.7/Scientific:/usr/include/python2.7:/usr/include/python2.7/numpy:$CPLUS_INCLUDE_PATH



PROG = pythonExt
TRGTS = DEFF SPOTSET BEAMLET DEFFDOUBLE SPOTSETDOUBLE BEAMLETDOUBLE

$(PROG): $(TRGTS)

all: $(TRGTS)


DEFF:
	cd ./PythonExtensionCalcDeff ; make DEFF
	@echo "Compilation PythonExtensionCalcDeff done"
	@echo "**********************************************************"
	@echo 
	@echo 
SPOTSET:
	cd ./PythonExtensionBuildSpotSet ; make SPOTSET
	@echo "Compilation PythonExtensionBuildSpotSet done"
	@echo "**********************************************************" 
	@echo 
	@echo 

BEAMLET:
	cd ./PythonExtensionBeamletCalc ; make BEAMLET
	@echo "Compilation PythonExtensionBeamletCalc done"
	@echo "**********************************************************" 
	@echo 
	@echo 
	
	
DEFFDOUBLE:
	cd ./PythonExtensionCalcDeffDouble ; make DEFFDOUBLE
	@echo "Compilation PythonExtensionCalcDeffDouble done"
	@echo "**********************************************************" 
	@echo 
	@echo 
SPOTSETDOUBLE:
	cd ./PythonExtensionBuildSpotSetDouble ; make SPOTSETDOUBLE
	@echo "Compilation PythonExtensionBuildSpotSetDouble done"
	@echo "**********************************************************" 
	@echo 
	@echo 
BEAMLETDOUBLE:
	cd ./PythonExtensionBeamletCalcDouble ; make BEAMLETDOUBLE
	@echo "Compilation PythonExtensionBeamletCalcDouble done"
	@echo "**********************************************************" 
	@echo 
	@echo 

clean:
	rm -f *.o *~
	cd ./PythonExtensionCalcDeff ; make clean
	cd ./PythonExtensionBuildSpotSet ; make clean
	cd ./PythonExtensionBeamletCalc ; make clean
	cd ./PythonExtensionCalcDeffDouble ; make clean
	cd ./PythonExtensionBuildSpotSetDouble ; make clean
	cd ./PythonExtensionBeamletCalcDouble ; make clean


cleanSO:
	rm -f *.so
	cd ./PythonExtensionCalcDeff ; make cleanSO
	cd ./PythonExtensionBuildSpotSet ; make cleanSO
	cd ./PythonExtensionBeamletCalc ; make cleanSO
	cd ./PythonExtensionCalcDeffDouble ; make cleanSO
	cd ./PythonExtensionBuildSpotSetDouble ; make cleanSO
	cd ./PythonExtensionBeamletCalcDouble ; make cleanSO


# End of the main math makefile

