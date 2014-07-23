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
import math

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

class DistanceFileParser(object):
    """img_16S_gut_oral_97percent.csv"""
    # constants for readability
    _IMG_ID_1     = 0
    _IMG_ID_2     = 1
    _IDENTITY_16S = 2
    
    def __init__(self):
        self.prepped = False

    def readFile(self, # this is a generator function
                 fh):
        while True:
            if not self.prepped:
                # we still need to strip out the header
                for l in fh:
                    if l[0] == "I": # header line
                        self.prepped = True
                        break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[1],
                        float(fields[2])]
                       )
            break # done! 

        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class IDFileParser(object):
    """img_ids_ordered_by_genome_tree.Mar14.csv"""
    # constants
    _GT_ID  = 0 
    _IMG_ID = 1 
    
    def __init__(self):
        self.prepped = False
    
    def readFile(self, # this is a generator function
                 fh):
        while True:
            if not self.prepped:
                # we still need to strip out the header
                for l in fh:
                    if l[0] == "g": # header line
                        self.prepped = True
                        break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[1]])
            break # done! 


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class DistanceData(object):
    """class to store genome comparison data, with 16S distance"""
    def __init__(self):
        self.idLookUp           = {} # IMG_ID -> GT_ID
        self.comparisons        = {}
        self.roundedComparisons = {}
        
    def addIDS(self,
               IMG_ID,
               GT_ID
               ):
        """create a lookup table with IMG_ID -> GT_ID"""
        self.idLookUp[IMG_ID] = GT_ID 
    
    def addComparison(self,
                      IMG_ID_1,
                      IMG_ID_2,
                      perc16S):
        """Call this function after creating the lookup table for IDs"""
        try: 
            self.comparisons[self.idLookUp[IMG_ID_1]][self.idLookUp[IMG_ID_2]] = math.ceil(perc16S)
        except KeyError:
            try:
                self.comparisons[self.idLookUp[IMG_ID_1]] = {self.idLookUp[IMG_ID_2] : math.ceil(perc16S)}
            except KeyError:
                pass
        try:
            self.comparisons[self.idLookUp[IMG_ID_2]][self.idLookUp[IMG_ID_1]] = math.ceil(perc16S)
        except KeyError:
            try:
                self.comparisons[self.idLookUp[IMG_ID_2]] = {self.idLookUp[IMG_ID_1] : math.ceil(perc16S)}
            except KeyError:
                pass
    
    def collapse16S(self):
        """create data structure: rounded 16S perc -> no. of comparisons"""
        for id_1 in self.comparisons.keys():
            for id_2 in self.comparisons[id_1]:
                try:
                    self.roundedComparisons[self.comparisons[id_1][id_2]] += 1 # add one to the count of comparisons 
                except KeyError:
                    self.roundedComparisons[self.comparisons[id_1][id_2]] = 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################


class HitData(object):
    """class for capturing hit data"""
    def __init__(self):
        self.hits               = {} # dict to store hit data
        self.length             = {} # dict to store length data
        self.distance           = {} # dict to store 16S distance
        self.roundedDistance    = {} # dict to store rounded 16S distance and total hits per 100 comparisons
        self.standardDeviation  = {} # dict to store standard deviation at each percentage
    
    def addLen(self,
               _ID_1,
               _ID_2,
               _LGT_LEN,
               ):   
        """create cumulative contig length object"""
        try: 
            self.length[_ID_1][_ID_2] += _LGT_LEN
        except KeyError:
            try:
                self.length[_ID_1][_ID_2]  = _LGT_LEN
            except KeyError:
                self.length[_ID_1] = {_ID_2 : _LGT_LEN}
        
        try: 
            self.length[_ID_2][_ID_1] += _LGT_LEN
        except KeyError:
            try:
                self.length[_ID_2][_ID_1] = _LGT_LEN
            except KeyError:
                self.length[_ID_2] = {_ID_1 : _LGT_LEN}
                    
    def add16SDist(self,
                   _ID_1,
                   _ID_2,
                   _PERC_ID
                   ): 
        """create 16S distance object, round up %"""
        try: 
            self.distance[_ID_1][_ID_2] = math.ceil(_PERC_ID)
        except KeyError:
            self.distance[_ID_1] = {_ID_2 : math.ceil(_PERC_ID)}
        try: 
            self.distance[_ID_2][_ID_1] = math.ceil(_PERC_ID)
        except KeyError:
            self.distance[_ID_2] = {_ID_1 : math.ceil(_PERC_ID)}
                
    def addHit(self,
               _ID_1,
               _ID_2,
               ):
        """add an LGT instance to the data store"""
        try: 
            self.hits[_ID_1][_ID_2] += 1
        except KeyError:
            try:
                self.hits[_ID_1][_ID_2] = 1
            except KeyError:
                self.hits[_ID_1] = {_ID_2 : 1}  

            
    def getIDS(self):
        """return array of IDs"""
        return self.hits.keys()

    def groupBy16S(self):
        """calculate the number of transfers per rounded 16S score"""
        for id_1 in self.hits.keys():
            for id_2 in self.hits[id_1]:
                try:
                    #self.roundedDistance[self.distance[id_1][id_2]] += [self.hits[id_1][id_2]]
                    self.roundedDistance[self.distance[id_1][id_2]] += self.hits[id_1][id_2]     # total hits per percentage
                    self.standardDeviation[self.distance[id_1][id_2]] += [self.hits[id_1][id_2]] # hit array per percentage
                except KeyError:
                    #self.roundedDistance[self.distance[id_1][id_2]] = [self.hits[id_1][id_2]]
                    self.roundedDistance[self.distance[id_1][id_2]] = self.hits[id_1][id_2]     # total hits per percentage
                    self.standardDeviation[self.distance[id_1][id_2]] = [self.hits[id_1][id_2]] # hit array per percentage
                  
    def calculateStD(self,perc):
        """return the standard deviation for each percentage"""
        hitList = []
        hitList = self.standardDeviation[perc]
        print perc,hitList
        stdev = np.array(hitList)
        #print  str(np.std(stdev, dtype=np.float64))
        return  np.std(stdev, dtype=np.float64) 
                    
    def numHits16S(self):
        """print stats about the number of hits in the roundedDistance dict"""
        for perc in self.roundedDistance.keys():
            totalPerPercent = 0
            length = len(self.roundedDistance[perc])
            for i in self.roundedDistance[perc]:
                totalPerPercent = totalPerPercent + i
            print "********************************************************************"
            print "%i" % (perc)
            print "Total hits: %i" % (totalPerPercent)
            print "Number of comparisons: %i" % (length)
            print "Average hits: %f" % (totalPerPercent/float(length))
            print "Min number of hits: %i" % (min(self.roundedDistance[perc]))
            print "Max number of hits: %i" % (max(self.roundedDistance[perc]))
            print "********************************************************************"
            
    def return16SPercArray(self):
        """return array of rounded 16S percentages """
        return self.roundedDistance.keys()
        
    def returnNormValues(self):
        """return array of hits normalised by 100 comparisons per rounded 16S score"""
        return 
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    def __init__(self,
                 transfersFile
                 ):
        self.transfersFile = transfersFile
        
    def readFile(self):
        """read data from csv file, and capture as object"""
        HFP = HitFileParser()
        self.HD = HitData()
        with open(self.transfersFile,'r') as fh:
            for hit in HFP.readHit(fh):
                self.HD.add16SDist(hit[HFP._ID_1], hit[HFP._ID_2], hit[HFP._PERC_ID])
                self.HD.addHit(hit[HFP._ID_1], hit[HFP._ID_2])
                self.HD.addLen(hit[HFP._ID_1], hit[HFP._ID_2], hit[HFP._LGT_LEN])
        self.workingIDs = self.HD.getIDS() # working ids list   
        self.HD.groupBy16S() # create dictionary of rounded 16S distance scores
        
        #print self.HD.roundedDistance[88]
        #print self.HD.standardDeviation[88]

        
    def connect(self):
        """Try connect to the TrackM server"""
        pass

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
                      comparisonsFile):
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
        
        #normalise hits per 100 comparisons
        percList = self.DD.roundedComparisons.keys() # list of percentages in DistanceData
        percList.sort()
        normalisedHits = []
        for perc in percList:
            try:
                #print str(perc),self.DD.roundedComparisons[perc], self.HD.roundedDistance[perc]
                x = self.DD.roundedComparisons[perc]/float(100)
                normalised = self.HD.roundedDistance[perc] / float(x)
            except KeyError: # no hits at that percentage
                x = 1
                normalised = 0 
            normalisedHits.append(normalised)
             
            # calculate standard deviation
            try:
                print perc, self.HD.calculateStD(perc)
            except KeyError:
                pass
            
            
        x,y = percList,normalisedHits
        #print x
        #print y
        
        
        
        
        # Build plot
        plt.scatter(x, y, marker='|', s=1000)
        plt.plot(x,y, linestyle='-')
        plt.axis([100,75,0,100])
        plt.xlabel('16S distance (%)')
        plt.ylabel('HGT per 100 comparisons')
        plt.show() # plot 
    
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################