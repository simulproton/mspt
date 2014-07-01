########################################################################
#
# plotDataGraphs.py
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
import numpy as np
import os, sys
import matplotlib
if not 'DISPLAY' in os.environ.keys():
            matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import time

globalFont =  FontProperties()
globalFont.set_size('small')
globalFont.set_style('italic')
globalFont.set_family('fantasy')
#family = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
#style  = ['normal', 'italic', 'oblique']
#size: 'xx-small, x-small, small, medium, large, x-large, xx-large or an absolute font sze: 12 ...


def drawData(listXData,listYData,listNames = None, display = True, title = None, xLabel = None, yLabel = None , posLegend = None ):
    '''Draws the several set of (X,Y) data in a pyplot figure. The first xValues array and the first yValues array\
    is plot. Then the second and so on.
    
    If display is True, the figure will be displayed.
    
    :param listXdata: List of 1D arrays of x Values
    :param listYData: List of 1D arrays of y Values
    :param listNames: List of names to be used with each set of (X,Y) data. If listNames is None (default)\
    the names are defined automatically as "Data_IDX" where IDX is their index in the lists.
    :param display: True (default) to display the window, False otherwise.
    :param title: String for the figure title. None by default.
    :param xLabel: title for x axis. None by default.
    :param yLabel: title for y axis. None by default.
    :param posLegend: Legend position: ['right', 'bottom']. By default it is set to None which is considered as 'right'.
    
    '''
    if listXData is None or listYData is None:
        print "Nothing to draw in drawData..."
        return
    else:
        if display:
            plt.ion()
            plt.clf()
            plt.show()
                
        figure = getCurrentFigure()
    
        [xMin,xMax] = findminMax(listXData)
        [yMin,yMax] = findminMax(listYData)
        xMax = xMax + xMax*10/100
        yMax = yMax + yMax*10/100
        count = 0
        plotNum = 1
        subplot = newSubplot(figure)
        clearSubplot(subplot)
        for idx,(xData,yData) in enumerate(zip(listXData,listYData)):
            if listNames is not None:
                name = listNames[idx]
            else:
                name = "Data_%0.2d"%idx
            drawPlot(subplot , xData, yData,xmin = xMin, xmax = xMax, ymin = yMin, ymax =yMax,xlabel = xLabel, ylabel=yLabel, titlePlot = name, singlePlot = True, nameSinglePlot = title, grid = True)
            title = None
            if display:
                displayPlot()
                    
                    
        makePlotLegend(figure,subplot,pos =  posLegend)
        if display:
            displayPlot()
        return (figure,subplot)
        
        
        
                
def findminMax(listData):
    '''Function that finds the overall minimum and maximum values in a list of data arrays.
    
    :param listData: List of 1D arrays.
    
    :returns: List : [min, max]
    
    '''
    if listData is not None:
        min = None
        max = None
        for data in listData:
            values = data
            localMin = np.min(values)
            localMax = np.max(values)
            
            if min is None:
                min = localMin
            else:
                if localMin < min:
                    min = localMin
                
            if max is None:
                max = localMax
            else:
                if localMax > max:
                    max = localMax          
    else:
        raise ValueError ("list data is None in findXMinMax") 
    return [min , max]


def newSubplot(fig , multiple = False , num = 1):
    '''Creates a new subplot in a figure.
    
    :param fig: Reference of the matplotlib figure.
    :param multiple: True to draw multiple subplots in the figure. False (default) otherwise. 
    :param num: Number of desired subplots (used only if multiple is True). If there are multiple subplots in the figure\
    it is set such that there are 2 rows of subplots and n columns.
    
    :returns: A reference to the subplot created and added to the figure.
    
    
    '''
    if multiple:
        if (num - 1)%2 == 1: 
            sub = fig.add_subplot( 2,1 + ((num-1) - 1)/2, num) 
        else:
            sub = fig.add_subplot( 2, 1 + (num-1) /2, num)
    else:
        sub = fig.add_subplot(111) 
    return sub

def drawPlot(subplot , xData, yData,xmin = None, xmax = None, ymin = None, ymax = None,xlabel = None, ylabel=None, titlePlot = None, singlePlot = False, nameSinglePlot = None, grid = True): 
    '''Draws the data (xData,yData) in subplot. 
    
    :param xmin: minimum value of the x axis. None by default.
    :param xmax: maximum value of the x axis. None by default.
    :param ymin: minimum value of the y axis. None by default.
    :param ymax: maximum value of the y axis. None by default.
    :param xlabel: title for the x axis. None by default.
    :param ylabel: title for the y axis. None by default.
    :param titlePlot: title for the subplot. None by default.
    :param singlePlot: True if a single plot is in the figure. False by default.
    :param nameSinglePlot:  title for the figure if there is only one plot (singlePlot is True). None by default.
    :param grid: True (default) or False if we want to display or not a grid in the subplot
    
    '''
    subplot.plot( xData, yData, label = titlePlot )
    if xmin is not None:
        subplot.set_xlim(xmin = xmin)
    if xmax is not None:
        subplot.set_xlim(xmax = xmax)
    if ymin is not None:
        subplot.set_ylim(ymin = ymin)
    if ymax is not None:
        subplot.set_ylim(ymax = ymax)
    if xlabel is not None:
        subplot.set_xlabel(xlabel)
    if ylabel is not None:
        subplot.set_ylabel(ylabel)
    
    if not singlePlot:
        if titlePlot is not None:
            subplot.set_title(titlePlot)
    else:
        if nameSinglePlot is not None:
            subplot.set_title(nameSinglePlot)
    if grid:
        subplot.grid(b='on')


def displayPlot():
    '''Draws the window with the figure and the subplots. It is displayed during 2 seconds.
    
    '''
    plt.draw()
    time.sleep(2)
 
def getCurrentFigure():
    '''Get the reference of the current figure. Equivalent to pyplot.gcf().
    
    :returns: The reference of the current figure.
    '''
    return plt.gcf()
 
def closePlot():
    '''Closes a plotting area when done.
    
    '''
    plt.close()

def clearSubplot(subplot):
    '''Clears the subplot.
    
    :param subplot: Reference to the subplot.
    
    '''
    subplot.clear()

    
def makePlotLegend(fig , subplot, pos = None):
    '''Builds a legend for the subplot in a figure.
        
    :param fig: reference to the figure
    :param subplot: reference to the subplot. 
    :param pos: Legend position: ['right', 'bottom']. By default None, which corresponds to 'right'
    
    '''
    if pos is None:
        pos = 'right'
    elif pos not in ['right', 'bottom']:
        print "In makePlotLegend: unkown argument pos = %s. pos = 'right' or 'bottom' "%str(pos)
        pos = 'right'
    if pos == 'right':
        #Put legend on the right
        box = subplot.get_position()
        subplot.set_position([box.x0, box.y0, box.width * 0.75, box.height* 0.75])
        leg = subplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    elif pos == 'bottom':
        #Put the legend under the plot:
        box = subplot.get_position()
        subplot.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
        subplot.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=5)
    else:    
        box = subplot.get_position()
        subplot.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
        subplot.legend(loc=pos, bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=5)

def savePlotImage(fig, filename , dirPath = None, ext = 'png', subplot = None):
    '''Saves a figure to image file (png or pdf).
    
    :param fig: reference to the figure to save
    :param filename: image name
    :param dirPath: directory where to save the image. By default it is None: dirPath is set to os.getcwd() + "/PlotImages/".
    :param ext: file extension (png or pdf). By default: png
    :param subplot: Reference to the subplot. If the figure contains several subplots and the user wants to save only 1 subplot it needs to specify it.
    
    '''
    #format = png or pdf
    if ext in ['png','pdf']:
        if not filename.endswith('.'+str(ext)):
            filename = filename + '.' + str(ext)       
    else:
        print "Wrong extension: %s. Should be pdf or png"%str(ext)
        print "Plot will be saved as png."
        if not filename.endswith('.png'):
            filename = filename + '.png' 
    if  dirPath is None:
        pathToSave = os.getcwd() + "/PlotImages/"
    else:
        pathToSave = dirPath
    if not os.path.exists(pathToSave):
        os.makedirs(pathToSave)
    
    plt.figure(fig.number)
    
    if subplot is None:
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
        extent = subplot.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
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

    print "Image saved: %s"%((pathToSave +str(filename)))
    
def plotHorizontalLine(subplot, yVal , xEnd = None ,color = None , label = None,font = None):
    '''Plots a horizontal line in the subplot at yVal. By default it draws the line from y axis to the right of the subplot, but\ 
    a x value (corresponding to the real xValue) can be specified if the user wants to stop it before.
    
    :param yVal: Y value where to draw the line.
    :param xEnd: X value where to end the line. By default it is None.
    :param color: Matplotlib color used to draw the line and the label if a label is provided. By default black.
    :param label: Text being displayed at the left the y axis. 
    :param font: Font properties if a label is defined. By default it is None and uses matplotlib default font
    
    '''
    if color is None:
            color = 'k'
    # Label: text being displayed
    if label is not None:
        lab = label # label cannot be used in text() because it is a keyword
        if font is None:
            font = FontProperties()  
        subplot.text(-0.5,yVal,lab,rotation=0,fontproperties=font,color = color)
    
    # Darw the line    
    if xEnd is None:
        xEnd = 1
    else:
        #Set xEnd in range 0-1:
        xMax = subplot.get_xlim()[1]
        xEnd = xEnd / xMax
    subplot.axhline(y=yVal,xmin=0, xmax=xEnd,color = color,linestyle = 'dashed')  
        
    
        
def plotVerticalLine(subplot, xVal, yEnd = None , color = None, label = None, font = None):
    '''Plots a vertical line in the subplot at xVal. By default it draws the line from x axis to the top of the subplot, but 
    a y value (corresponding to the real yValue) can be specified if the user wants to stop it before.

    :param xVal: X value where to draw the line.
    :param yEnd: Y value where to end the line. By default it is None.
    :param color: Matplotlib color used to draw the line and the label if a label is provided. By default black.
    :param label: Text being displayed at the left the y axis. 
    :param font: Font properties if a label is defined. By default it is None and uses matplotlib default font

    
    '''
    if color is None:
            color = 'k'
    # Label: text being displayed
    if label is not None:
        lab = label # label cannot be used in text() because it is a keyword
        if font is None:
            font = FontProperties()        
        subplot.text(xVal,-5.5,lab,rotation=0,fontproperties=font,color = color)
     
    # Darw the line
    if yEnd is None:
        yEnd = 1
    else:
        #Set xEnd in range 0-1:
        yMax = subplot.get_ylim()[1]
        yEnd = yEnd / yMax        
    subplot.axvline(x=xVal,ymin=0, ymax=yEnd,color = color,linestyle = 'dashed') 