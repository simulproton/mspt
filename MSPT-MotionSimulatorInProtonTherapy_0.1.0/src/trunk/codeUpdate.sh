#!/bin/sh
svn up
echo "svn up Done"

cd mspt/extensionsC/PythonExtensionAllPlatforms/
make
echo "C++ code updated"
make clean

cd ../../..
echo "Code update done..." 