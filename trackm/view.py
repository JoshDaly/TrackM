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

# local imports
from trackm.viewInterface import ViewInterface, Condition


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ImagePropertiesGeneral(object):
    """Standard image output properties"""
    def __init__(self,
                 imageFormat,
                 xlabel,
                 ylabel,
                 title,
                 outFile,
                 dpi,
                 labelFontSize,
                 markerStyle,
                 markerSize,
                 edgeColour
                 ):
        self.imageFormat   = imageFormat
        self.xLabel        = xLabel
        self.yLabel        = yLabel
        self.title         = title
        self.outFile       = outFile
        self.dpi           = dpi 
        self.labelFontSize = labelFontSize
        self.markerStyle   = markerStyle
        self.markerSize    = markerSize
        self.edgeColour    = edgeColour 

class ImagePropertiesFrequency(object,
                               ImagePropertiesGeneral):
    """options specific for frequency plot"""
    def __init__(self
                 ):
        pass

        
class ImagePropertiesScatter(object,
                             ImagePropertiesGeneral):
    """options specific for scatter plot """
    def __init__(self
                 ):
        pass
        
class ImagePropertiesNetwork(object,
                             ImagePropertiesGeneral):
    """options specific for network plot """
    def __init__(self,
                 nodeColour, # colour of nodes, can be array of colours of len(x)
                 nodeFontColour, # colour of node labels
                 nodeFontSize, # size of node labels 
                 nodeLabels # labels on/off
                 ):
        ImagePropertiesGeneral.__init__(self, imageFormat, xlabel, ylabel, title, outFile, dpi, labelFontSize, markerStyle)
        self.nodeColour     = nodeColour
        self.nodeFontColour  = nodeFontColour
        self.nodeFontSize    = nodeFontSize
        self.nodeLabels      = nodeLabels      

###############################################################################
###############################################################################
###############################################################################
###############################################################################


class Plotter(object):
    """Class for plotting functions
    Scatter
    Network
    Frequency"""
    
    def __init__(self):
        self.HD = HitData()
    

    def scatterPlot(self):
        """Produces a scatter plot indicating the number of lgt events

           between each pair of genomes, as well as the cumulative length

           of DNA transferred
        """

        # objects
        IPS = ImagePropertiesScatter()
        
        """ Produce scatter plot """
        xs = []         # hits
        ys = []         # cumulative contigs lengths
        workingIds = [] # master list of genome tree ids

        #Plot labels
        plt.xlabel('Cumulative length of contigs (bp)')
        plt.ylabel('No. of hits')
        plt.title("Gut-Oral LGT events")
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
            plt.savefig(output_file,dpi=args.dpi,format=str(image_format))


    def networkPlot(self):
        """Produces a network plots with nodes representing individual genomes

           and edges representing an LGT event between the two gneomes
        """
        HD = HitData()
        """Creating the network graph"""
        G=nx.Graph()
        G.add_nodes_from(self.workingIDs)
        # graph variables
        edgeWidth=[]
        edgeColour=[]
        phylumCols = []
        for id_1 in self.HD.hits.keys():
            for id_2 in self.HD.hits[id_1]:
                #print id_1,id_2,str(self.hits[id_1][id_2])
                G.add_edge(id_1,id_2) # loop through dict, and add edges

        # edit edge properties # This needs to be edited to include the number of hits/genome size.
        #for edge in G.edges():
        #    edgeWidth.append(int(ids_dict[edge[0]][edge[1]][0]))

        #values = [phylum_cols.get(node,0.25) for node in G.nodes()]

        pos= nx.spring_layout(G)
        #nx.draw(G,pos,node_color=values,with_labels=args.labels,width=edgeWidth,font_size=12,font_color='#006400')
        nx.draw(G,pos,with_labels=True,font_size=12,font_color='#006400')
        #nx.draw(G,pos)
        plt.show()


    def frequencyPlot(self,
                      lookUpFile,
                      comparisonsFile,
                      dirtyFile):
        """Produces a line graph showing the frequency of LGT between genomes

           per 100 comparisons relative the ANI distance between the two genomes

        NEEDS TO INCORPORATE ALL COMPARISONS, NOT JUST THE ONES THAT HAD AN LGT EVENT!!!!
        """
        # objects
        DFP = DistanceFileParser()
        IDFP = IDFileParser()
        self.DD = DistanceData()

        # read in id look up file
        with open(lookUpFile,'r') as fh:
            for hit in IDFP.readFile(fh):
                self.DD.addIDS(hit[IDFP._IMG_ID], hit[IDFP._GT_ID])

        # read in comparisons file
        with open(comparisonsFile,'r') as fh:
            for hit in DFP.readFile(fh):
                 self.DD.addComparison(hit[DFP._IMG_ID_1], hit[DFP._IMG_ID_2], hit[DFP._IDENTITY_16S])
        self.DD.collapse16S() # calculate no. of comparisons at each rounded 16S distance

        # read in dirty hit file
        DHFP = DirtyHitFileParser()

        with open(dirtyFile,'r') as fh:
            for hit in DHFP.readHit(fh):
                self.DD.addDirtyHit(hit[DHFP._ID_1], hit[DHFP._ID_2])
        self.DD.getDirty16S() # creates 16S -> hits

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
        plt.xlabel('16S distance (%)')
        plt.ylabel('LGT per 100 comparisons')
        plt.show() # plot 
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    def __init__(self,
                 dbFileName
                 ):
        self.dbFileName = dbFileName
        """
                 serverURL,         # URL of the commanding TrackM server
                 serverPort         # port of the commanding TrackM server
                 ):
        self.serverURL = serverURL
        self.serverPort = serverPort
        """
        
    def testSomething(self,
                      ani=1.,             # only get hits with this ani or less
                      batch=None          # only get hits from this batch (None == all)
                      ):
        # get an interface to the DB
        VI = ViewInterface(self.dbFileName)

        # build the condition we want to select on and get the hit data
        C = Condition("ani_comp", "<=", ani)
        if batch is not None:
            C = Condition(C, "and", Condition("batch", "=", batch))

        hits = VI.getHitData(C)
        print hits
        

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
