########################################################################
#
# zipFiles.py
# Proton Therapy Simulator Project
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee
# July 2013
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



import zipfile, os


def createZipArchiveForFiles(listFiles, zipname, rmCommonPath = True):
    '''Creates a zip archive from a list of files.
    
    :param listFiles: list of paths to files
    :param zipname: name of the archive
    :param rmCommonPath: True if you want to remove the common prefixe path to all the files when they are being archived
    
    :returns: A list [Bool , String], where bool is True if the zip process was achieved, False otherwise. The string correspond to an informative message.
    
    '''
    if listFiles == []:
        return [False , 'No files to archive... Archiving failed..']
        
    if not zipname.endswith(".zip"):
        zipname = zipname + '.zip'
        
    zipf = zipfile.ZipFile(zipname, 'w', zipfile.ZIP_DEFLATED,allowZip64 = True)
    if rmCommonPath:
        if len(listFiles) > 1:
            common = os.path.commonprefix(listFiles)
        else:
            common = os.path.split(listFiles[0])[0]
    for f in listFiles:
        if rmCommonPath:
            newName = os.path.relpath(f, start = common)
            if len(listFiles) == 1:
                newName = os.path.split(os.path.splitext(zipname)[0])[1] +'/'+ newName
            zipf.write(f,arcname = newName)
        else:
            zipf.write(f)
    zipf.close()
    return [True ,'Archive %s done'%zipname]


def getFileListForDirectory( dirPath , visitSubDir = True , listExt = None):
    '''Creates a list of files from a directory path: the files are in the directory or its subdirectories.
    
    :param dirPath: path to the directory
    :param visitSubDir: True if we list the files contained in the subdirectories
    :param listExt: list of extensions, if we want to list only files with given extensions
    
    :returns: A list of files.
    
    '''
    listFiles = []
    for file in os.listdir(dirPath):
        newPath = os.path.join(dirPath, file)
        if os.path.isfile(newPath):
            if  listExt is None:
                listFiles.append(newPath)
            else:
                flag = 0
                for ext in listExt:
                    if not ext.startswith('.'):
                        ext = '.'+ext
                    if os.path.splitext(newPath)[1] == ext:
                        flag = 1
                        break
                if flag == 1:
                    listFiles.append(newPath)
        elif os.path.isdir(newPath) and visitSubDir:
            subList = getFileListForDirectory( newPath , visitSubDir  , listExt)
            listFiles.extend(subList)
    return listFiles



def archiveDataFromListPath(listPath , zipname, goToSubDir = True ,listExt = None, rmCommonPath = True):
    '''From given directory and file path create a zip archive.
    
    :param listPath: list of path of files and directories
    :param zipname: name of the archive
    :param goToSubDir: if True, when the path is a directory list the files in the subdirectories
    :param listExt: list of extensions we want add in the zip archive. Files with other extensions won't be added
    :param rmCommonPath: True if you want to remove the common prefixe path to all the files when they are being archived
    
    :returns: Informative text on the success or the failure.
    
    '''
    #Make sure that all extensions start with a ".".
    if listExt is not None:
        newListExt = []
        for item in listExt:
            if not item.startswith('.'):
                newExt = '.'+item
                newListExt.append(newExt)
            else:
                newListExt.append(item)
    else:
        newListExt = listExt
    #Create a list of all the files that need to be added to the archive        
    listFiles = []
    for path in listPath:
        if os.path.isdir(path):
            listFiles.extend(getFileListForDirectory( path , goToSubDir , newListExt))
        elif os.path.isfile(path):
            if newListExt is not None:
                extFile = os.path.splitext(path)[1] 
                if extFile in newListExt:
                    listFiles.append(path)
            else:
                listFiles.append(path)
    [success,text] = createZipArchiveForFiles(listFiles, zipname, rmCommonPath)     
    print text
