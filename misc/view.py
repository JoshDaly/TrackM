#!/usr/bin/env python
###############################################################################
#                                                                             #
#    view.py [development version]                                            #
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
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# local imports
from trackm.exceptions import *


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HitFileParser(object):
    """class for parsing hit data files"""
    # constants for readability
    _UID     = 0  
    _ID_1    = 1
    _LGT_LEN = 2
    _ID_2    = 3
    _PERC_ID = 4
    
    def __init__(self):
        self.prepped = False
    
    def readHit(self, # this is a generator function
                fh 
                ):
        while True:
            if not self.prepped:
                # we still need to strip out the header
                for l in fh:
                    if l[0] == "u": # header line
                        self.prepped = True
                        break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[2],
                        fields[7],
                        fields[11],
                        float(fields[-1])]
                       )
            break # done! 


###############################################################################
###############################################################################
###############################################################################
###############################################################################


class HitData(object):
    """class for capturing hit data"""
    def __init__(self):
        self.hits       = {} # dict to store hit data
        self.length     = {} # dict to store length data
        self.distance   = {} # dict to store 16S distance
    
    def addLen(self,
               _ID_1,
               _ID_2,
               _LGT_LEN,
               ):   
        """create cumulative contig length object"""
        try: 
            self.length[_ID_1][_ID_2] += _LGT_LEN
        except KeyError:
            self.length[_ID_1] = {}
            self.length[_ID_1][_ID_2] = _LGT_LEN
        try: 
            self.length[_ID_2][_ID_1] += _LGT_LEN
        except KeyError:
            self.length[_ID_2] = {}
            self.length[_ID_2][_ID_1] = _LGT_LEN
                    
    def add16SDist(self,
                   _ID_1,
                   _ID_2,
                   _PERC_ID
                   ): 
        """create 16S distance object"""
        try: 
            self.distance[_ID_1][_ID_2] = _PERC_ID
        except KeyError:
            self.distance[_ID_1] = {}
            self.distance[_ID_1][_ID_2] = _PERC_ID
        try: 
            self.distance[_ID_2][_ID_1] = _PERC_ID
        except KeyError:
            self.distance[_ID_2] = {}
            self.distance[_ID_2][_ID_1] = _PERC_ID
                
    def addHit(self,
               _ID_1,
               _ID_2,
               ):
        """add an LGT instance to the data store"""
        try: 
            self.hits[_ID_1][_ID_2] += 1
        except KeyError:
            self.hits[_ID_1] = {_ID_2 : 1}
        try: 
            self.hits[_ID_2][_ID_1] += 1
        except KeyError:
            self.hits[_ID_2] = {_ID_1 : 1}
          
    def printHits(self):  
        print "*****************"
        for ida in self.hits.keys():
            for idb in self.hits[ida]:
                print ida,idb 
        print "*****************"
            
    def getIDS(self):
        return self.hits.keys()

        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    def __init__(self,
                 transfersFile,
                 chartType):
        self.transfersFile = transfersFile
        self.chartType = chartType
        
    def readFile(self):
        HFP = HitFileParser()
        HD = HitData()
        with open(self.transfersFile,'r') as fh:
            for hit in HFP.readHit(fh):
                HD.add16SDist(hit[HFP._ID_1], hit[HFP._ID_2], hit[HFP._PERC_ID])
                HD.addHit(hit[HFP._ID_1], hit[HFP._ID_2])
                HD.addLen(hit[HFP._ID_1], hit[HFP._ID_2], hit[HFP._LGT_LEN])
        HD.printHits()
        self.workingIDs = HD.getIDS() # working ids list   
        self.hits = HD.hits           # hits dictionary 
        self.distance = HD.distance   # 16S distance dictionary
        self.length = HD.length       # cumulative contig length dictionary
        
    def connect(self):
        """Try connect to the TrackM server"""
        pass

    def plotSomething(self):
        """Place holder"""
        # an example entry point from the main TrackM entry point
        # Ideally, you should have a separate plot function
        # defined for each type of plot.

        # >>>>>>>>>> REMOVE THIS WHEN YOU HAVE CODE HERE <<<<<<<<<<<<<<<<<<
        from inspect import currentframe, getframeinfo
        frameinfo = getframeinfo(currentframe())
        print "Make me plot something! I live at File: %s Line: %s" % (frameinfo.filename, frameinfo.lineno)
        # >>>>>>>>>> END <<<<<<<<<<<<<<<<<<

    def scatterPlot(self):
        """Produces a scatter plot indicating the number of lgt events
        
           between each pair of genomes, as well as the cumulative length
           
           of DNA transferred
        """
        pass
        
        
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
        for id_1 in self.hits.keys():
            for id_2 in self.hits[id_1]:
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
        
        
    def lgtFrequency(self):
        """Produces a line graph showing the frequency of LGT between genomes
        
           per 100 comparisons relative the ANI distance between the two genomes
        """
        pass
        
        
    
    
    
    
    
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################