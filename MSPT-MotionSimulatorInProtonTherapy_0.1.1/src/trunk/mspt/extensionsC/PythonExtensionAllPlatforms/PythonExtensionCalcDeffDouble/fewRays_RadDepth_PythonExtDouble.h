/*
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
*/
#include "Python.h" // MUST be the first library included!!
#include "arrayobject.h" // from python numpy library 
#include <math.h>
#include "calcRadiolDepthDouble.h"



/* ==== Prototypes =================================== */

// .... Python callable functions ..................

static PyObject *fewRaysRadDepthDouble(PyObject *self, PyObject *args);

int * fillDimensionTable( PyArrayObject * array);

