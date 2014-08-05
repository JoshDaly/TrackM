#!/usr/bin/env python
###############################################################################
#                                                                             #
#    view.py                                                                  #
#                                                                             #
#    The View class handles all TrackM based visualisations.                  #
#    NOTE: This class expects a TrackM server to be running in 'daemon' mode. #
#    when it is started.                                                      #
#                                                                             #
#    Copyright (C) Josh Daly, Michael Imelfort                                #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Josh Daly, Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system imports
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
np.seterr(all='raise')
# local imports
from trackm.viewInterface import ViewInterface, Condition


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ImagePropertiesGeneral(object):
    """Standard image output properties"""
    def __init__(self,
                 showPlot,
                 imageFormat,
                 xLabel,
                 yLabel,
                 title,
                 outFile,
                 dpi,
                 labelFontSize,
                 edgeColour
                 ):
        self.showPlot      = showPlot
        self.imageFormat   = imageFormat
        self.xLabel        = xLabel
        self.yLabel        = yLabel
        self.title         = title
        self.outFile       = outFile
        self.dpi           = dpi 
        self.labelFontSize = labelFontSize
        self.edgeColour    = edgeColour 

class ImagePropertiesFrequency(ImagePropertiesGeneral):
    """options specific for frequency plot"""
    def __init__(self,
                 showPlot,
                 imageFormat, 
                 xLabel, 
                 yLabel, 
                 title, 
                 outFile, 
                 dpi, 
                 labelFontSize, 
                 markerStyle, 
                 markerSize,
                 edgeColour
                 ):
        ImagePropertiesGeneral.__init__(self, showPlot, imageFormat, xLabel, yLabel, title, outFile, dpi, labelFontSize, edgeColour)
        self.markerStyle   = markerStyle
        self.markerSize    = markerSize

        
class ImagePropertiesScatter(ImagePropertiesGeneral):
    """options specific for scatter plot """
    def __init__(self,
                 showPlot,
                 imageFormat, 
                 xLabel, 
                 yLabel, 
                 title, 
                 outFile, 
                 dpi, 
                 labelFontSize, 
                 markerStyle, 
                 markerSize,
                 edgeColour
                 ):
        ImagePropertiesGeneral.__init__(self, showPlot, imageFormat, xLabel, yLabel, title, outFile, dpi, labelFontSize, edgeColour)
        self.markerStyle   = markerStyle
        self.markerSize    = markerSize
        
class ImagePropertiesNetwork(ImagePropertiesGeneral):
    """options specific for network plot """
    def __init__(self,
                 showPlot,
                 imageFormat, 
                 xLabel, 
                 yLabel, 
                 title, 
                 outFile, 
                 dpi, 
                 labelFontSize, 
                 nodeShape,
                 nodeSize, 
                 edgeColour,
                 nodeColour, # colour of nodes, can be array of colours of len(x)
                 nodeFontColour, # colour of node labels
                 nodeFontSize, # size of node labels 
                 nodeLabels # labels on/off
                 ):
        ImagePropertiesGeneral.__init__(self, showPlot, imageFormat, xLabel, yLabel, title, outFile, dpi, labelFontSize, edgeColour)
        self.nodeColour      = nodeColour
        self.nodeFontColour  = nodeFontColour
        self.nodeFontSize    = nodeFontSize
        self.nodeLabels      = nodeLabels      
        self.nodeSize        = nodeSize
        self.nodeShape       = nodeShape


    def __str__(self):
        print "****************************************"
        print self.dpi
        print self.imageFormat


###############################################################################
###############################################################################
###############################################################################
###############################################################################


class Plotter(object):
    """Class for plotting functions
       Scatter
       Network
       Frequency
    """
    
    def __init__(self,
                 hitData):
        self.HD = hitData
    
        """
        hitData objects:
        self.hitCounts          = {} # dict to store hit data
        self.hitLengthTotals    = {} # dict to store length data
        self.workingIds         = {} # dict used Ids { gid -> bool } True == has hit
        self.identityAni        = {} # dict to store ANI distance
        """

    def makeScatterPlot(self):
        """Produces a scatter plot indicating the number of lgt events

           between each pair of genomes, as well as the cumulative length

           of DNA transferred
        """

        # objects
        IPS = ImagePropertiesScatter()
        
        """ Produce scatter plot """
        xs = []         # hits
        ys = []         # cumulative contigs lengths
        workingIds = self.HD.workingIds # master list of genome tree ids

        #Plot labels
        plt.xLabel(IPS.xLabel)
        plt.yLabel(IPS.yLabel)
        plt.title(IPS.title)
        #ax.grid(b=True,which='both')

        #set the number of tick labels
        for i in range(len(working_ids_list)):
            ticks.append(i)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        #Change label names according to taxonomy
        labels = [item.get_text() for item in ax.get_xticklabels()]
        for i,v in enumerate(labels):
            labels[i] = ordered_tax_string_lowest[i]
        ax.set_xticklabels(labels,size=args.text_size)
        ax.set_yticklabels(labels,size=args.text_size)

        #rotate x labels 90deg
        plt.xticks(rotation=90)

        #adjust margin size
        #plt.subplots_adjust(bottom=0.2)
        plt.tight_layout()

        #plt.xticks(np.arange(min(xs), max(xs)+1,1.0))
        #plt.yticks(np.arange(min(xs), max(xs)+1,1.0))
        # set the plot axes
        plt.axis([plot_border*-1,
                  len(working_ids_list)+plot_border,
                  plot_border*-1,
                  len(working_ids_list)+plot_border])

        if args.show_plot == "True":
            plt.show()
        else:
            plt.savefig(IPS.outFile,dpi=IPS.dpi,format=IPS.imageFormat)


    def makeNetworkPlot(self,
                        showPlot,
                        imageFormat, 
                        xLabel, 
                        yLabel, 
                        title, 
                        outFile, 
                        dpi, 
                        labelFontSize, 
                        nodeShape, 
                        nodeSize,
                        edgeColour,
                        nodeColour, # colour of nodes, can be array of colours of len(x)
                        nodeFontColour, # colour of node labels
                        nodeFontSize, # size of node labels 
                        nodeLabels, # labels on/off
                        ):
        """Produces a network plots with nodes representing individual genomes

           and edges representing an LGT event between the two gneomes
        """
        # objects
        IPN = ImagePropertiesNetwork(showPlot,
                                     imageFormat, 
                                     xLabel, 
                                     yLabel, 
                                     title, 
                                     outFile, 
                                     dpi, 
                                     labelFontSize, 
                                     nodeShape, 
                                     nodeSize,
                                     edgeColour,
                                     nodeColour, # colour of nodes, can be array of colours of len(x)
                                     nodeFontColour, # colour of node labels
                                     nodeFontSize, # size of node labels 
                                     nodeLabels, # labels on/off
                                     )

        """Creating the network graph"""
        G=nx.Graph()
        G.add_nodes_from(self.HD.workingIds)
        # graph variables
        edgeWidth=[]
        edgeColour=[]
        phylumCols = []
        for id_1 in self.HD.hitCounts.keys():
            for id_2 in self.HD.hitCounts[id_1]:
                #print id_1,id_2,str(self.hits[id_1][id_2])
                G.add_edge(id_1,id_2) # loop through dict, and add edges

        # edit edge properties # This needs to be edited to include the number of hits/genome size.
        #for edge in G.edges():
        #    edgeWidth.append(int(ids_dict[edge[0]][edge[1]][0]))

        #values = [phylum_cols.get(node,0.25) for node in G.nodes()]

        pos= nx.spring_layout(G)
        #nx.draw(G,pos,node_color=values,with_labels=args.labels,width=edgeWidth,font_size=12,font_color='#006400')
        nx.draw(G,pos,with_labels=IPN.nodeLabels,font_size=IPN.nodeFontSize,font_color=IPN.nodeFontColour, node_color = IPN.nodeColour, node_size = IPN.nodeSize)
        #nx.draw(G,pos)
        if showPlot:
            plt.show()
        else:
            plt.savefig(IPS.outFile,dpi=IPS.dpi,format=IPS.imageFormat)
        
        

    def makeFrequencyPlot(self,
                          showPlot,
                          imageFormat,
                          xLabel,
                          yLabel,
                          title,
                          outFile,
                          dpi,
                          labelFontSize,
                          markerStyle, 
                          markerSize,
                          edgeColour
                          ):
        """Produces a line graph showing the frequency of LGT between genomes

           per 100 comparisons relative the ANI distance between the two genomes

        NEEDS TO INCORPORATE ALL COMPARISONS, NOT JUST THE ONES THAT HAD AN LGT EVENT!!!!
        """
        
        """
        hitData objects:
        self.hitCounts          = {} # dict to store hit data
        self.hitLengthTotals    = {} # dict to store length data
        self.workingIds         = {} # dict used Ids { gid -> bool } True == has hit
        self.identityAni        = {} # dict to store ANI distance
        """
        
        # objects
        IPF = ImagePropertiesFrequency(showPlot,
                                       imageFormat,
                                       xLabel,
                                       yLabel,
                                       title,
                                       outFile,
                                       dpi,
                                       labelFontSize,
                                       markerStyle, 
                                       markerSize,
                                       edgeColour)
        roundedAllComparisonsDec = {} # dict to store no. of comparisons at each ANI, to the decimal
        roundedAllComparisonsInt = {} # dict to store no. of comparisons at each ANI, to the integer
        percentListInt           = []
        percentListDec           = []
        roundedHitsInt           = {}
        roundedHitsDec           = {}

        # create total number of comparisons dict
        for id1 in self.HD.identityAni.keys():
            for id2 in self.HD.identityAni[id1]:
                if float(self.HD.identityAni[id1][id2]) != -1.0: 
                    try:
                        roundedAllComparisonsDec[round(float(self.HD.identityAni[id1][id2]),1)] += 1
                    except KeyError:
                        roundedAllComparisonsDec[round(float(self.HD.identityAni[id1][id2]),1)] = 1 
                    
                    try:
                        roundedAllComparisonsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] += 1
                    except KeyError:
                        roundedAllComparisonsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] = 1
                    if round(float(self.HD.identityAni[id1][id2]),1) not in  percentListDec:
                        percentListDec.append(round(float(self.HD.identityAni[id1][id2]),1))
                    if math.ceil(float(self.HD.identityAni[id1][id2])) not in percentListInt:
                        percentListInt.append(math.ceil(float(self.HD.identityAni[id1][id2])))
                    
        # create total number of hits at each percentage
        for id1 in self.HD.hitCounts.keys():
            for id2 in self.HD.hitCounts[id1]:
                try: 
                    roundedHitsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] += self.HD.hitCounts[id1][id2]
                except KeyError:
                    roundedHitsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] = self.HD.hitCounts[id1][id2]
                try:
                    roundedHitsDec[round(float(self.HD.identityAni[id1][id2]),1)] += self.HD.hitCounts[id1][id2]
                except KeyError:
                    roundedHitsDec[round(float(self.HD.identityAni[id1][id2]),1)] = self.HD.hitCounts[id1][id2]                    
                    
        percentListInt.sort()
        percentListDec.sort()
        
        intList = np.arange(0,96,1)
        decList = np.arange(0,96,0.1)
        
        percentListIntY = []
        percentListDecY = []
        
        percentHitsIntY = []
        percentHitsDecY = []
        
        # append zeros to lists
        for i in range(96):
            percentListIntY.append(0)
            percentHitsIntY.append(0)
        #for i in range(96):
        #    percentListDecY.append(0)
        #    percentHitsDecY.append(0)
        
        for i,v in enumerate(intList):
            try:
                percentListIntY[i] = roundedAllComparisonsInt[v]
            except KeyError:
                pass
        
        #for i,v in enumerate(decList):
        #    try:
        #        percentListDecY[i] = roundedAllComparisonsDec[v]
        #    except KeyError:
        #        pass
        
        
        # get y values for hits
        for i,perc in enumerate(intList):
            try:
                c = roundedAllComparisonsInt[perc] / float(100)
                percentHitsIntY[perc] = roundedHitsInt[perc] / float(c)
                #print perc, roundedAllComparisonsInt[perc], roundedHitsInt[perc], percentHitsIntY[perc]
            except KeyError:
                pass
        
        # plot 
        plt.scatter(intList, percentHitsIntY, marker='|')
        plt.plot(intList, percentHitsIntY, linestyle='-')
        
        # plot no. of comparisons
        plt.plot(intList,percentListIntY,c='#FFFFFF')
        plt.fill_between(intList, percentListIntY, 1e-6, facecolor = '#C0C0C0')
        plt.axis([100,-5,0,max(percentHitsIntY)+10])
        #plt.show()
        
        #normalise hits per 100 comparisons
        percList = self.DD.roundedComparisons.keys() # list of percentages in DistanceData
        percList.sort()
        normalisedDirtyHits = []
        standardisedDirty = []
        normalisedHits = []
        standardised = []
        comparisons = []

        for perc in percList:
            normalise = 0
            normaliseDirty = 0
            # array of comparisons
            comparisons.append(self.DD.roundedComparisons[perc]/float(100))

            # create x and y coordinates for line graph
            c = self.DD.roundedComparisons[perc]/float(100)
            try:
                normalise = self.HD.roundedDistance[perc] / float(c)
            except KeyError:
                normalise = 0
            normalisedHits.append(normalise)
            try:
                normaliseDirty = self.DD.dirtyRoundedHits[perc] / float(c)
            except KeyError:
                normaliseDirty = 0
            normalisedDirtyHits.append(normaliseDirty)
            # calculate standard deviation
            try:
                hitList   = self.HD.standardDeviation[perc] # list of hits
                fullHitList = self.DD.roundedComparisons[perc] # int
                zerosToAdd = fullHitList - len(hitList)
                for i in range(zerosToAdd):
                    hitList.append(0)
                #print hitList
                stdDev = np.array(hitList)
                standardised.append(np.std(stdDev, dtype=np.float64))
                #print standardised
            except KeyError:
                standardised.append(0)

        x,y    = percList,normalisedHits # clean
        xd, yd = percList,normalisedDirtyHits # dirty
        xc, yc = percList,comparisons
        # print out data
        for i in range(len(percList)):
            print percList[i],comparisons[i],normalisedHits[i]




        # Build plot
        # comparisons
        plt.plot(xc,yc,c='#FFFFFF')
        plt.fill_between(xc, yc, 1e-6, facecolor = '#C0C0C0')
        # clean
        plt.scatter(x, y, marker='|')
        plt.plot(x,y, linestyle='-')
        # dirty
        plt.scatter(xd, yd, marker='|') # dirty
        plt.plot(xd,yd, linestyle='-')
        # add error bars
        for i in range(len(x)):
            plt.plot([x[i],x[i]],[y[i]-standardised[i],y[i]+standardised[i]],'k')
        plt.axis([100,75,0,80])
        plt.xLabel('16S distance (%)')
        plt.yLabel('LGT per 100 comparisons')
        plt.show() # plot 
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    """visualise database"""
    def __init__(self,
                 dbFileName,
                 ani=95,
                 batch=None
                 ):
        # get an interface to the DB
        VI = ViewInterface(dbFileName)
        
        # build the condition we want to select on and get the hit data
        C = Condition("ani_comp", "<=", ani)
        if batch is not None:
            C = Condition(C, "and", Condition("batch", "=", batch))
        # get hitData object
        self.hits = VI.getHitData(C)
        
        # call plotting functions
        self.PLOT = Plotter(self.hits)
        
        """
                 serverURL,         # URL of the commanding TrackM server
                 serverPort         # port of the commanding TrackM server
                 ):
        self.serverURL = serverURL
        self.serverPort = serverPort
        """
        """
        hitData objects:
        self.hitCounts          = {} # dict to store hit data
        self.hitLengthTotals    = {} # dict to store length data
        self.workingIds         = {} # dict used Ids { gid -> bool } True == has hit
        self.identityAni        = {} # dict to store ANI distance
        """
    def testSomething(self
                      ):
        count = 0 
        #for key in self.hits.hitCounts.keys():
        #    for key2 in self.hits.hitCounts[key]:
        #        print key, key2, self.hits.hitCounts[key][key2]
                
        #print self.hits.identityAni
        
        
    def scatterPlot(self,
                    showPlot,
                    imageFormat, 
                    xLabel, 
                    yLabel, 
                    title, 
                    outFile, 
                    dpi, 
                    labelFontSize, 
                    markerStyle,
                    markerSize, 
                    edgeColour
                    ):
        """Produce a scatter plot"""
        IPS = ImagePropertiesScatter(showPlot,
                                     imageFormat, 
                                     xLabel, 
                                     yLabel, 
                                     title, 
                                     outFile, 
                                     dpi, 
                                     labelFontSize, 
                                     markerStyle,
                                     markerSize, 
                                     edgeColour
                                     )    
        self.PLOT.makeScatterPlot()
    
    def networkPlot(self,
                    showPlot,
                    imageFormat, 
                    xLabel, 
                    yLabel, 
                    title, 
                    outFile, 
                    dpi, 
                    labelFontSize, 
                    nodeShape,
                    nodeSize, 
                    edgeColour,
                    nodeColour, # colour of nodes, can be array of colours of len(x)
                    nodeFontColour, # colour of node labels
                    nodeFontSize, # size of node labels 
                    nodeLabels # labels on/off
                    ):
        """Produce a network plot"""
    
        self.PLOT.makeNetworkPlot(showPlot,
                                  imageFormat, 
                                  xLabel, 
                                  yLabel, 
                                  title, 
                                  outFile, 
                                  dpi, 
                                  labelFontSize, 
                                  nodeShape, 
                                  nodeSize,
                                  edgeColour,
                                  nodeColour, # colour of nodes, can be array of colours of len(x)
                                  nodeFontColour, # colour of node labels
                                  nodeFontSize, # size of node labels 
                                  nodeLabels # labels on/off
                                  )
    
    def frequencyPlot(self,
                      showPlot,
                      imageFormat, 
                      xLabel, 
                      yLabel, 
                      title, 
                      outFile, 
                      dpi, 
                      labelFontSize, 
                      markerStyle,
                      markerSize, 
                      edgeColour
                      ):
        """Produce a frequency plot"""

        self.PLOT.makeFrequencyPlot(showPlot,
                      imageFormat, 
                      xLabel, 
                      yLabel, 
                      title, 
                      outFile, 
                      dpi, 
                      labelFontSize, 
                      markerStyle,
                      markerSize, 
                      edgeColour)
    
    def connect(self):
        """Try connect to the TrackM server"""
        pass

    def plotSomething(self):
        """Place holder"""
        # an example entry point from the main TrackM entry point
        # Ideally, you should have a separate plot function
        # defined for each type of plot.


###############################################################################
###############################################################################
###############################################################################
###############################################################################
