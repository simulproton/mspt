########################################################################
#
# array2DViewer.py
# 
# Created by Paul Morel, LIGM, Universite Paris-Est Marne La Vallee, France
# paul.morel@univ-mlv.fr
# June 2012
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
#   Code used to save a 2D numpy array to a PNG or Tiff image
#   A file is also stored with each image file: ".met". It is a file readable
#   with cPickle library and is intended to store the data used to convert actual intensities
#   to the intensity range used in the image, i.e. 0-255 for PNG and 0-65535 for Tiff
########################################################################
import PIL.Image as Image
import numpy as np
import cPickle
import os
import matplotlib as mpl



def XServer_is_running():
    '''Code used to test if the X11 server is running. This is used to know how configure the displays in matplotlib
    
    :returns: True or False
    
    '''
    from subprocess import Popen, PIPE
    p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0


if not 'DISPLAY' in os.environ.keys():
            mpl.use("Agg")
            print "No display"
else: #Some problems appeared for the display on some servers ([ "c-iibi000" ])  
        # when connecting by ssh without using X11 on the client and when the ssh options
        # -X or -Y are not used: 'DISPLAY' was in os.environ.keys() whereas it should not
        # So the lines right after intend to fix this problem. However, there is no problem is 
        # the server is accessed using ssh -X or -Y. No other fix has been found yet.
    import socket
    currentServer = socket.gethostname()
    if  currentServer in [ "c-iibi000" ]:
        if not  XServer_is_running():# Test if X11 is running on the server.
            mpl.use("Agg")
            print "No display"


import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
import matplotlib.image as mpimg
import matplotlib.cm as cm
import time

    
def saveDataToImagePNG( data , filename,min_max=[],bw = True):
    '''Save a 2D array to an image file: png if the data is in the range 0-255 and tiff otherwise. \ 
    The array is linearly scaled to be in the range 0-255 (png) or 0-65535 (tiff). This is determined by the maximum value stored in min_max.
    
    :param data: 2D numpy array
    :param filename: path and name of the image to be saved. The extension (.png or .tiff) will be added in the function.
    :param min_max: By default min_max = []. In this case the maximum image intensity is automatically set to 255 or 65535 depending on what\
    the maximum value is in data. If min_max = [min_Val , max_Val], the linear scaling determines the slope and the intercept to be used \
    based on these min_max values. Max_Val will determine if the image will be a png (0-255) or a tiff (0-65535).
    :param bw: True (Default) / False. True to save the image in black and white (grayscale), False to save it in color. Note that images saved \
    in tiff format are automatically saved in black-and white.
    
    '''
    
    if min_max == []:
        mini = np.min(data)
        maxi = np.max(data)
    else:
        mini = min_max[0]
        maxi = min_max[1]
          
    if mini != 0 or maxi!= 255:
        if maxi < 255:
            limit = 255 #pow(2,8)
        else: 
            limit = 65535 #pow(2,16)
        
        if mini != maxi:
            slope =  limit/ (maxi - mini)
            intercept = -slope * mini
        else:
            slope = 1.0
            if mini < 0:
                intercept = 0
            elif maxi > limit:
                intercept = limit
            else:
                intercept = 0
        if limit == 255:
            scaledData_int = np.zeros((data.shape),dtype=np.uint8)
            scaledData_int[:,:] = (scaledData_int[:,:] +  slope * (data) + intercept)
            if not bw:
                dataRange01 = scaledData_int / 255.0 
                coloredData = cm.jet(dataRange01)
                scaledData_int = np.zeros(( coloredData.shape),dtype=np.uint8)
                scaledData_int[:,:,:] = np.uint8(coloredData*255.0)[:,:,:]
                
        else:
            scaledData_int = np.zeros((data.shape),dtype=np.uint16)
            scaledData_int[:,:] = (scaledData_int[:,:] +  slope * (data) + intercept)
            if not bw:
                pass

        if limit == 255:
            img = Image.fromarray(scaledData_int)
            img.save((filename+'.png'))
        else:
            img = Image.fromarray(scaledData_int,'I;16')
            img.save((filename+'.tiff'))            

    else:
        data_int = data.astype(np.uint8)
        img = Image.fromarray(data_int)
        img.info["dpi"] = 360
        img.save((filename+'.png'))
        
        
def openPNGImageToArray(filename):
    '''Loads a PNG image and returns it as a numpy array.
    
    :param filename: The image file to open.
    
    :returns: A 2D numpy array corresponding to the input image.
    
    '''
    img = Image.open(filename)
    fileName, fileExtension = os.path.splitext(filename)
    data = np.asarray(img)
    return data
        
def displayImageFromArray( data , name = None):
    '''Displays an image from a numpy array using matplotlib. It uses "jet" colormap and \
    the figure opens and closes automatically.
    
    :param data: 2D numpy array.
    :param name: If name is not None, it is used as the title of the figure. 
    
    '''
    plt.ion()
    grid = plt.subplots(1, 1)
    if name is not None:
        plt.title(name)
   
    img = plt.subplot(111).imshow(data,cmap = 'jet')
    cbar = grid[0].colorbar(img, cmap='jet')
    
    displayPlot()

    



def displayPlot():
    '''Draws the window with the figure and the subplots. It is displayed during 10 seconds.
    
    '''
    plt.draw()
    if 'DISPLAY' in os.environ.keys():
        time.sleep(10)


def closePlot():
    '''Closes the matplotlib plotting area.
    
    '''
    plt.close('all')
 
    
def saveFigure( pathToSave,filename , ext = 'png', subplotIdx = None):
    '''Saves a matplotlib figure to an image file (png or pdf).
    
    :param pathToSave: Path where to save the image.
    :param filename: Image name.
    :param ext: File extension (png or pdf). By default: png.
    :param subplotIdx: If the figure contains several subplots and the user wants to save only 1 subplot it needs to specify it.\
    For example if the plots are in a grid 2x3 and one wants to plot the 4th subplot, subplotIdx = 234. \
    If there is a single subplot: 111. It is None by default.
    
    '''
    fig = plt.gcf()
    #format = png or pdf
    if ext in ['png','pdf']:
        if not filename.endswith('.'+str(ext)):
            filename = filename + '.' + str(ext)       
    else:
        print "Wrong extension: %s. Should be pdf or png"%str(ext)
        print "Plot will be saved as png."
        if not filename.endswith('.png'):
            filename = filename + '.png' 
    
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)
    
    plt.figure(fig.number)
    if subplotIdx is None:
        try:
            plt.savefig(pathToSave +filename,bbox_inches='tight',dpi = 300,transparent=True)
        except:
            try:
                plt.savefig(pathToSave +filename[0:-4]+'.pdf',bbox_inches='tight',dpi = 300,transparent=True)
            except:
                print "Something went wrong: image : %s could not be saved."%(pathToSave +filename)
                return
            print "Sorry, something went wrong. Image was saved as a pdf."
            filename = filename[0:-4]+'.pdf'
    else:
        extent = plt.subplot(subplotIdx).get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        try:
            fig.savefig(pathToSave +filename , bbox_inches=extent.expanded(1.5, 1.5),dpi = 360,transparent=True) 
        except:
            try:
                fig.savefig(pathToSave +filename[0:-4]+'.pdf', bbox_inches=extent.expanded(1.7, 1.7),dpi = 360,transparent=True)
            except:
                print "Something went wrong: image : %s could not be saved."%(pathToSave +filename)
                return
            print "Sorry, something went wrong. Image was saved as a pdf."
            filename = filename[0:-4]+'.pdf'
    print "Image %s saved"%((pathToSave +filename))



def storeMatrixAsImage(pathToSave,data,name, alongAxis = None,cpkl = False, bw = True,  minmax = None):
    '''Stores a 3D numpy array as a stack of images and stores the numpy array data in a cPickle structure through the function 'storeData()'.
    Slices can be stored along a given axis: 'z' dicom patient axis, 'y' dicom patient axis, 'x' dicom patient axis. \
    By default 'z' axis is chosen. 
    
    :param pathToSave: Path where to save the images.
    :param data: 3D numpy array to save.
    :param name: Name for the images to save.
    :param alongAxis: List of the axis along which the images will be saved. By default the value is None which corresponds to 'z'. \
    Possible values are 'x' (i.e. the images are stored along the column axis of the dicom image), 'y' (i.e. the images are stored \
    along the rows axis of the dicom image), 'z' (i.e. the images are stored along the frame axis of the dicom image). This can be used\
    when the patient is moving and the user wants to see the motions from the side, the top or the front of the dicom image.
    :param cpkl: True to save the 3D numpy array to a cPickle structure. This allows to have have access to the numpy array outside the simulator. \
    False otherwise. Default: False.
    :param bw: True to save the images in black and white. False otherwise. Default: True
    :param minmax: Min and Max values used to scale the image intensities. Default: None, in this case the min and max used are the min and max \
    of the numpy array.
    
    '''
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)
    if cpkl:
        filename = pathToSave + str(name)+".cpkl"
        storeData(filename,data)

    if alongAxis != 'z' and alongAxis is not None:
        if alongAxis == 'x':
            print "%s : Slices being stored along X axis"%name
            data = np.transpose(data , axes=(2,1,0) )
            name = name +"_alongX"
        elif alongAxis == 'y':
            print "%s : Slices being stored along Y axis"%name
            data = np.transpose(data , axes=(1,0,2) )
            name = name +"_alongY"
        else:
            print " Unrecognized 'alongAxis' value : %s. Slices will be sotred along Z axis."%(str(alongAxis))
            name = name +"_alongZ"
    else:
        print "%s : Slices being stored along Z axis"%name
        name = name +"_alongZ"
    path = pathToSave + str(name)+"/"
    if not os.path.exists(path):
        os.makedirs(path)
                
    countFrame = 0
    if minmax is None:
        mini = np.min(data)
        maxi = np.max(data)
    else:
        mini = minmax[0]
        maxi = minmax[1]
        
    if maxi > 255 and not bw:
        print "Warning - image saved will be stored in tiff and color in not supported yet. It will be saved in black/white."
        bw = True
    for frame in data:
        strFrame = "%04d"%countFrame
        filename = str(path) + str(name)+"-frame-"+strFrame
        saveDataToImagePNG(frame , filename,min_max=[mini,maxi],bw = bw)
        countFrame = countFrame+1    

def storeData(filename, data):
    '''Stores some data to a cPickle structure. This allows to open the data, as it would be in the simulator, outside the simulator.
    
    :param filename: Path where to save the cPickle file.
    :param data: Any python data to save.
    
    '''
    cPickle.dump(data, open(filename,'wb'), protocol=cPickle.HIGHEST_PROTOCOL)
    print "\t\tData stored to cPickle struct.: %s"%filename 
    
def loadData( filename):
    '''Loads cPickle data to retrieve the python data stored in it.
    
    :param filename: Path to the file to load.
    
    :returns: The data loaded.
    
    '''
    fileData = open(filename,"rb")
    data = cPickle.load(fileData)
    print "\tData loaded from cPickle struct.: %s"%filename 
    return data