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
import random
from colorsys import hsv_to_rgb as htr
np.seterr(all='raise')
from Bio import SeqIO
from Bio.Seq import Seq


# local imports
from trackm.viewInterface import ViewInterface, Condition
from trackm.cb2cols import Cb2Cols as CB2

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

class MetaDataFileParser(object):
    """class for parsing img metadata files"""
    _GT_ID      = 0
    _STATUS     = 1
    _PHYLUM     = 2
    _BODY_SITE  = 3
    _SEQ_CENTRE = 4
    _SEQ_PLAT   = 5
        
    def __init__(self):
        self.prepped = False
        
    def readFile(self, # this is a generator function
                 fh
                 ):
        while True:
            if not self.prepped:
                # we still need to strip out the header
                for l in fh:
                    if l[0] == "g":
                        self.prepped = True
                        break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield ([fields[0],
                        fields[3],
                        fields[7],
                        fields[57],
                        fields[6],
                        fields[-1]
                        ])
            break # done! 

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class MetaData(object):
    """class for storing metadata"""
    def __init__(self):
        self.metadata = {} # gid -> {phylum, bodysite, seqCentre, seqPlatform}

    def getMetaData(self, 
                    metadataFile
                    ):
        MFP = MetaDataFileParser() # call generator
        with open(metadataFile, 'r') as fh:
            for meta in MFP.readFile(fh):
                self.metadata[meta[MFP._GT_ID]] = [meta[MFP._PHYLUM],meta[MFP._BODY_SITE], meta[MFP._SEQ_PLAT], meta[MFP._SEQ_CENTRE], meta[MFP._STATUS]]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CleanHitData(object):
    def __init__(self,):
        self.hitCounts          = {} # dict to store hit data
        self.hitLengthTotals    = {} # dict to store length data
        self.workingIds         = {} # dict used Ids { gid -> bool } True == has hit
        
    def getData(self,
                fastaFile,
                hitData):
        """read in fasta file"""
        """accession = >pid_sqid"""
        # parse fasta file using biopython
        for accession,sequence in SeqIO.to_dict(SeqIO.parse(fastaFile,"fasta")).items():
            pid  = int(accession.split("_")[0])
            try:
                gid1 = hitData.pidLookUp[pid][0] # maintains order of gids
                gid2 = hitData.pidLookUp[pid][1] 
                sqid = accession.split("_")[1]
                length = len(sequence)
                
                # add to working ids
                self.workingIds[gid1] = True
                self.workingIds[gid2] = True 
                
                # build hit hash
                try:
                    self.hitCounts[gid1][gid2] += 1 
                except KeyError: 
                    try:
                        self.hitCounts[gid1][gid2] = 1
                    except KeyError:
                        self.hitCounts[gid1] = {gid2: 1}
                
                # build contig length hash
                try:
                    self.hitLengthTotals[gid1][gid2] += length 
                except KeyError: 
                    try:
                        self.hitLengthTotals[gid1][gid2] = length
                    except KeyError:
                        self.hitLengthTotals[gid1] = {gid2: length}
            except KeyError:
                print "key %s not found" % (pid)
        

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
                 hitData,
                 subSet,
                 metadata,
                 colourBy,
                 status,
                 cleanHitData
                 ):
        self.HD             = hitData
        self.MD             = MetaData() # calls metadata object
        self.MD.getMetaData(metadata)
        self.subSet         = subSet
        self.colourBy       = colourBy
        self.status         = status
        self.cleanHitData   = cleanHitData
        
        
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
                        weighted,
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
        
        # graph variables
        edgeWidth     = []
        edgeColour    = []
        phylumCols    = {}
        bodySiteCols  = {}
        seqPlatCols   = {}
        seqCentreCols = {}
        vales         = []
        # Get ColorBrewer colous
        cb2 = CB2()
        col_set = "qualSet1"
        ColBrewColours = cb2.maps[col_set].values()[0:10]
        

        if self.subSet == 0:
            print "                          Using all of the data!                               "
            print "*******************************************************************************"
            
            # building colour set for metadata
            phylumDict     = {}
            bodySiteDict   = {}
            seqPlatDict    = {}
            seqCentreDict  = {}
            workingIdsDict = {}
            workingIdsHits = []
            count          = 1 
            colours = self.ReturnColours(1, 100) # list of rgb colours
            
            #####################################################
            """Build clean data structure""" 
            #####################################################
            
            if len(self.cleanHitData.workingIds) > 0: 
                
                # only ids with hits, and finished genomes
                for id in self.cleanHitData.workingIds:
                    if self.status == "finished":
                        try:
                            if self.MD.metadata[str(id)][4] == "Finished": # only finished genomes! 
                                if self.cleanHitData.workingIds[id]:
                                    workingIdsHits.append(id)
                        except KeyError:
                            pass 
                    else:
                        if self.cleanHitData.workingIds[id]:
                            workingIdsHits.append(id) # all genomes!
                
            
            #####################################################
            #Build dirty data structure 
            #####################################################
            else: 
                # only ids with hits, and finished genomes
                for id in self.HD.workingIds:
                    if self.status == "finished":
                        try:
                            if self.MD.metadata[str(id)][4] == "Finished": # only finished genomes! 
                                if self.HD.workingIds[id]:
                                    workingIdsHits.append(id)
                        except KeyError:
                            pass 
                    else:
                        if self.HD.workingIds[id]:
                            workingIdsHits.append(id) # all genomes!   
            #######################
            # Add nodes to all plots
            #######################
            G.add_nodes_from(workingIdsHits) # add all nodes
            print "length of working Ids %d" % len(workingIdsHits)
                        
            #####################################
            """ colour by phylum"""
            #####################################
            for id in workingIdsHits:
                
                workingIdsDict[id] = 0
                try: 
                    phylumCols[str(id)] = phylumDict[self.MD.metadata[str(id)][0]]
                except KeyError:
                    i = 1 * count 
                    #phylumDict[self.MD.metadata[str(id)][0]] = colours[i]
                    #phylumCols[str(id)] = colours[i]
                    #count += 10 # adjust to change the colour range
                    if self.MD.metadata[str(id)][0] == "Bacteroidetes":
                        phylumDict[self.MD.metadata[str(id)][0]] =  ColBrewColours[0] #"#ea00ff"
                        phylumCols[str(id)] = ColBrewColours[0] #"#ea00ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Lentisphaerae":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[1] #"#ffb700" 
                        phylumCols[str(id)] = ColBrewColours[1] #"#ffb700"
                    elif self.MD.metadata[str(id)][0] == "Firmicutes":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[2] #"#0047ff"
                        phylumCols[str(id)] = ColBrewColours[2] #"#0047ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Spirochaetes":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[3] #"#14ff00"
                        phylumCols[str(id)] = ColBrewColours[3] #"#14ff00"
                        
                    elif self.MD.metadata[str(id)][0] == "Synergistetes":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[4] #"#6600CC" 
                        phylumCols[str(id)] = ColBrewColours[4] #"#6600CC"
                        
                    elif self.MD.metadata[str(id)][0] == "Actinobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[5] #"#ffff00" 
                        phylumCols[str(id)] = ColBrewColours[5] #"#ffff00"
                        
                    elif self.MD.metadata[str(id)][0] == "Tenericutes":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[6] #"#006600" 
                        phylumCols[str(id)] = ColBrewColours[6] #"#0080ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Fusobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = ColBrewColours[7] #"#00e0ff" 
                        phylumCols[str(id)] = ColBrewColours[7] #"#00e0ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Proteobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                        phylumCols[str(id)] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                        
                        
            for phylum in phylumDict.keys():
                print "phylum: %s Colour: %s" % (phylum, phylumDict[phylum])
            
            print "There are %d different phylums represented" % (len(phylumDict))
            print "Number of interacting genomes %d" % len(phylumCols)
            print "*******************************************************************************"
            phylumValues = [phylumCols.get(node,0.25) for node in G.nodes()]
                  
            #####################################
            """ colour by body site"""
            #####################################
            for id in workingIdsHits:
                workingIdsDict[id] = 0
                try: 
                    bodySiteCols[str(id)] = bodySiteDict[self.MD.metadata[str(id)][1]]
                except KeyError:
                    i = 1 * count 
                    #bodySiteDict[self.MD.metadata[str(id)][1]] = colours[i]
                    #bodySiteCols[str(id)] = colours[i]
                    #count += 7 # adjust to change the colour range
                    if self.MD.metadata[str(id)][1] == "plant":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[2] #"#ea00ff"
                        bodySiteCols[str(id)] = ColBrewColours[2] #"#ea00ff"
                        
                    elif self.MD.metadata[str(id)][1] == "Eye":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = '#00CCCC' #ColBrewColours[8] #"#006600" 
                        bodySiteCols[str(id)] = '#00CCCC' #ColBrewColours[8] #"#ffb700"
                        
                    elif self.MD.metadata[str(id)][1] == "Airways":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[0] #"#0047ff"
                        bodySiteCols[str(id)] = ColBrewColours[0] #"#0047ff"
                        
                    elif self.MD.metadata[str(id)][1] == "internal_organs":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[3] #"#14ff00"
                        bodySiteCols[str(id)] = ColBrewColours[3] #"#14ff00"
                        
                    elif self.MD.metadata[str(id)][1] == "Gastrointestinal tract":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[4] #"#6600CC" 
                        bodySiteCols[str(id)] = ColBrewColours[4] #"#6600CC"
                        
                    elif self.MD.metadata[str(id)][1] == "Blood":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[5] #"#ffff00" 
                        bodySiteCols[str(id)] = ColBrewColours[5] #"#ffff00"
                        
                    elif self.MD.metadata[str(id)][1] == "skin":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[6] #"#0080ff" 
                        bodySiteCols[str(id)] = ColBrewColours[6] #"#0080ff"
                        
                    elif self.MD.metadata[str(id)][1] == "Urogenital tract":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[7] #"#00e0ff" 
                        bodySiteCols[str(id)] = ColBrewColours[7] #"#00e0ff"
                        
                    elif self.MD.metadata[str(id)][1] == "Ear":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                        bodySiteCols[str(id)] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                    
                    elif self.MD.metadata[str(id)][1] == "Oral":
                        bodySiteDict[self.MD.metadata[str(id)][1]] = ColBrewColours[1] #"#ffb700"
                        bodySiteCols[str(id)] = ColBrewColours[1] #"#ff1e00"
                        
                        
            for bodysite in bodySiteDict.keys():
                print "bodysite: %s Colour: %s" % (bodysite, bodySiteDict[bodysite])
                
            print "There are %d different bodysites represented" % (len(bodySiteDict))
            print "Number of interacting genomes %d" % len(bodySiteCols)
            print "*******************************************************************************"
            bodysiteValues = [bodySiteCols.get(node,0.25) for node in G.nodes()]
            
            #####################################
            """colour by sequencing platform"""
            #####################################
            for id in workingIdsHits:
                workingIdsDict[id] = 0
                try: 
                    seqPlatCols[str(id)] = seqPlatDict[self.MD.metadata[str(id)][2]]
                except KeyError:
                    i = 1 * count 
                    #seqPlatDict[self.MD.metadata[str(id)][2]] = colours[i]
                    #seqPlatCols[str(id)] = colours[i]
                    #count += 5 # adjust to change the colour range
                    
                    
                    
                    if self.MD.metadata[str(id)][2].rstrip().rstrip() == "Illumina GAii" or self.MD.metadata[str(id)][2].rstrip() == "Illumina HiSeq 2000" or self.MD.metadata[str(id)][2].rstrip() == "Illumina": 
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[0] #"#FF0000"
                        seqPlatCols[str(id)] = ColBrewColours[0] #"#FF0000"
                        
                    elif self.MD.metadata[str(id)][2].rstrip() == "454" or self.MD.metadata[str(id)][2].rstrip() == "454 GS FLX Titanium" or self.MD.metadata[str(id)][2].rstrip() == "454 GS FLX"  or self.MD.metadata[str(id)][2].rstrip() == "454-GS-FLX" or self.MD.metadata[str(id)][2].rstrip() == "454-GS20" or self.MD.metadata[str(id)][2].rstrip() == "454-GS-FLX-Titanium" :
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[1] #"#FF00FF" 
                        seqPlatCols[str(id)] = ColBrewColours[2] #"#FF00FF"
                        
                    elif self.MD.metadata[str(id)][2].rstrip() == "Solexa":
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[2] #"#0000FF"
                        seqPlatCols[str(id)] = ColBrewColours[1] #"#0000FF"
                        
                    elif self.MD.metadata[str(id)][2].rstrip() == "Sanger":
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[3] #"#00FFFF"
                        seqPlatCols[str(id)] = ColBrewColours[3] #"#00FFFF" 
                    
                    elif self.MD.metadata[str(id)][2].rstrip() == "454; ABI3730":
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[4] #"#00FF00"
                        seqPlatCols[str(id)] = ColBrewColours[4] #"#00FF00"
                    
                    elif self.MD.metadata[str(id)][2].rstrip() == "454; Illumina" or self.MD.metadata[str(id)][2].rstrip() == "Illumina;454":
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[5] #"#FFFF00"
                        seqPlatCols[str(id)] = ColBrewColours[5] #"#FFFF00"
                        
                    elif self.MD.metadata[str(id)][2].rstrip() == "Unknown":
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = ColBrewColours[6] #"#006600"
                        seqPlatCols[str(id)] = ColBrewColours[6] #"#006600"
                    
                    else:
                        print " %s" % (self.MD.metadata[str(id)][2].rstrip())
                        seqPlatDict[self.MD.metadata[str(id)][2].rstrip()] = '#00CCCC' #ColBrewColours[8] #"#C6C6C6"
                        seqPlatCols[str(id)] = '#00CCCC' #ColBrewColours[8] #"#C6C6C6"
                        
            for seqPlat in seqPlatDict.keys():
                print "Sequencing platform: %s Colour: %s" % (seqPlat.rstrip(), seqPlatDict[seqPlat])
                
            print "There are %d different sequencing platforms represented" % (len(seqPlatDict))
            print "Number of interacting genomes %d" % len(seqPlatCols)
            print "*******************************************************************************"
            seqPlatValues = [seqPlatCols.get(node,0.25) for node in G.nodes()]
            
            #####################################
            """colour by sequencing centre"""
            #####################################
            for id in workingIdsHits:
                workingIdsDict[id] = 0
                try: 
                    seqCentreCols[str(id)] = seqCentreDict[self.MD.metadata[str(id)][3]]
                except KeyError:
                    i = 1 * count
                    #print  str(id), self.MD.metadata[str(id)][3]
                    #seqCentreDict[self.MD.metadata[str(id)][3]] = colours[i]
                    #seqCentreCols[str(id)] = colours[i]
                    #count += 1 # adjust to change the colour range
                    
                    if self.MD.metadata[str(id)][3] == "DOE Joint Genome Institute":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[0] #"#FF0000"
                        seqCentreCols[str(id)] = ColBrewColours[0] #"#FF0000"
                        
                    elif self.MD.metadata[str(id)][3] == "Washington University in St. Louis":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[1] #"#FF00FF" 
                        seqCentreCols[str(id)] = ColBrewColours[1] #"#FF00FF"
                        
                    elif self.MD.metadata[str(id)][3] == "Baylor College of Medicine":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[2] #"#0000FF"
                        seqCentreCols[str(id)] = ColBrewColours[2] #"#0000FF"
                        
                    elif self.MD.metadata[str(id)][3] == "Broad Institute":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[3] #"#00FFFF"
                        seqCentreCols[str(id)] = ColBrewColours[3] #"#00FFFF" 
                    
                    elif self.MD.metadata[str(id)][3] == "J. Craig Venter Institute":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[4] #"#00FF00"
                        seqCentreCols[str(id)] = ColBrewColours[4] #"#00FF00"
                    
                    elif self.MD.metadata[str(id)][3] == "Iowa State Univ":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[5] #"#FFFF00"
                        seqCentreCols[str(id)] = ColBrewColours[5] #"#FFFF00"
                        
                    elif self.MD.metadata[str(id)][3] == "Sanger Institute":
                        seqCentreDict[self.MD.metadata[str(id)][3]] = ColBrewColours[6] #"#006600"
                        seqCentreCols[str(id)] = ColBrewColours[6] #"#006600"
                    
                    else:
                        seqCentreDict[self.MD.metadata[str(id)][3]] = '#00CCCC' # ColBrewColours[8] #"#C6C6C6"
                        seqCentreCols[str(id)] = '#00CCCC' #ColBrewColours[8] #"#C6C6C6"
                        
                    
            for seqCent in seqCentreDict.keys():
                print "Sequencing Centre: %s Colour: %s" % (seqCent, seqCentreDict[seqCent])
                
            print "There are %d different sequencing centres represented" % (len(seqCentreDict))
            print "Number of interacting genomes %d" % len(seqCentreCols)
            print "*******************************************************************************"
            seqCentreValues = [seqCentreCols.get(node,0.25) for node in G.nodes()]
            
            
            
            #####################################
            """Add edges between interacting genomes"""
            #####################################
            
            interPhyla              = 0
            intraPhyla              = 0
            SpirochaeteInter        = 0 
            SpirochaeteIntra        = 0
            LentisphaeraeInter      = 0
            LentisphaeraeIntra      = 0
            FirmicutesInter         = 0
            FirmicutesIntra         = 0
            BacteroidetesInter      = 0
            BacteroidetesIntra      = 0
            SynergistetesInter      = 0
            SynergistetesIntra      = 0
            ActinobacteriaInter     = 0
            ActinobacteriaIntra     = 0
            TenericutesInter        = 0
            TenericutesIntra        = 0
            FusobacteriaInter       = 0
            FusobacteriaIntra       = 0
            ProteobacteriaInter     = 0
            ProteobacteriaIntra     = 0

            
            
            if len(self.cleanHitData.workingIds) > 0:
                print "##############"
                print "Using user-specified hit data"
                for id_1 in self.cleanHitData.hitCounts.keys():
                    for id_2 in self.cleanHitData.hitCounts[id_1]:
                        # Determine intra- vs inter- phyla transfers
                        ##############################################################################
                        if self.MD.metadata[str(id_1)][0] == self.MD.metadata[str(id_2)][0]:
                            intraPhyla += self.cleanHitData.hitCounts[id_1][id_2]
                            if self.MD.metadata[str(id_1)][0] == 'Spirochaetes':
                                SpirochaeteIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Lentisphaerae':
                                LentisphaeraeIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Firmicutes':
                                FirmicutesIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Bacteroidetes':
                                BacteroidetesIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Synergistetes':
                                SynergistetesIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Actinobacteria':
                                ActinobacteriaIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Tenericutes':
                                TenericutesIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Fusobacteria':
                                FusobacteriaIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Proteobacteria':
                                ProteobacteriaIntra += self.cleanHitData.hitCounts[id_1][id_2]
                            
                        else: 
                            interPhyla += self.cleanHitData.hitCounts[id_1][id_2]
                            if self.MD.metadata[str(id_1)][0] == 'Spirochaete' or self.MD.metadata[str(id_1)][1] == 'Spirochaetes':
                                SpirochaeteInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Lentisphaerae' or self.MD.metadata[str(id_1)][1] == 'Lentisphaerae':
                                LentisphaeraeInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Firmicutes' or self.MD.metadata[str(id_1)][1] == 'Firmicutes':
                                FirmicutesInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Bacteroidetes' or self.MD.metadata[str(id_1)][1] == 'Bacteroidetes':
                                BacteroidetesInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Synergistetes' or self.MD.metadata[str(id_1)][1] == 'Synergistetes':
                                SynergistetesInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Actinobacteria' or self.MD.metadata[str(id_1)][1] == 'Actinobacteria':
                                ActinobacteriaInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Tenericutes' or self.MD.metadata[str(id_1)][1] == 'Tenericutes':
                                TenericutesInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Fusobacteria' or self.MD.metadata[str(id_1)][1] == 'Fusobacteria':
                                FusobacteriaInter += self.cleanHitData.hitCounts[id_1][id_2]
                            elif self.MD.metadata[str(id_1)][0] == 'Proteobacteria' or self.MD.metadata[str(id_1)][1] == 'Proteobacteria':
                                ProteobacteriaInter += self.cleanHitData.hitCounts[id_1][id_2]
                        
                        ##############################################################################
                        ### REMOVE THIS ####
                        #print "\t".join([id_1,id_2,str(self.cleanHitData.hitCounts[id_1][id_2])])
                        ### REMOVE THIS ####
                        try: 
                            goku = workingIdsDict[id_1]
                            picolo = workingIdsDict[id_2]
                            if self.cleanHitData.hitCounts[id_1][id_2] >0:
                                if weighted:
                                    G.add_edge(id_1,
                                               id_2,
                                               weight = self.cleanHitData.hitCounts[id_1][id_2])
                                else: 
                                    G.add_edge(id_1,
                                               id_2)
                        except KeyError:
                            pass
            else: 
                for id_1 in self.HD.hitCounts.keys():
                    for id_2 in self.HD.hitCounts[id_1]:
                        if self.HD.hitCounts[id_1][id_2] >0:
                            # Determine intra- vs inter- phyla transfers
                            if self.MD.metadata[str(id_1)][0] == self.MD.metadata[str(id_2)][0]:
                                intraPhyla += self.HD.hitCounts[id_1][id_2]
                                if self.MD.metadata[str(id_1)][0] == 'Spirochaetes':
                                    SpirochaeteIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Lentisphaerae':
                                    LentisphaeraeIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Firmicutes':
                                    FirmicutesIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Bacteroidetes':
                                    BacteroidetesIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Synergistetes':
                                    SynergistetesIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Actinobacteria':
                                    ActinobacteriaIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Tenericutes':
                                    TenericutesIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Fusobacteria':
                                    FusobacteriaIntra += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Proteobacteria':
                                    ProteobacteriaIntra += self.HD.hitCounts[id_1][id_2]
                            else: 
                                interPhyla += self.HD.hitCounts[id_1][id_2]
                                if self.MD.metadata[str(id_1)][0] == 'Spirochaete' or self.MD.metadata[str(id_1)][1] == 'Spirochaetes':
                                    SpirochaeteInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Lentisphaerae' or self.MD.metadata[str(id_1)][1] == 'Lentisphaerae':
                                    LentisphaeraeInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Firmicutes' or self.MD.metadata[str(id_1)][1] == 'Firmicutes':
                                    FirmicutesInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Bacteroidetes' or self.MD.metadata[str(id_1)][1] == 'Bacteroidetes':
                                    BacteroidetesInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Synergistetes' or self.MD.metadata[str(id_1)][1] == 'Synergistetes':
                                    SynergistetesInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Actinobacteria' or self.MD.metadata[str(id_1)][1] == 'Actinobacteria':
                                    ActinobacteriaInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Tenericutes' or self.MD.metadata[str(id_1)][1] == 'Tenericutes':
                                    TenericutesInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Fusobacteria' or self.MD.metadata[str(id_1)][1] == 'Fusobacteria':
                                    FusobacteriaInter += self.HD.hitCounts[id_1][id_2]
                                elif self.MD.metadata[str(id_1)][0] == 'Proteobacteria' or self.MD.metadata[str(id_1)][1] == 'Proteobacteria':
                                    ProteobacteriaInter += self.HD.hitCounts[id_1][id_2]
                            G.add_edge(id_1,id_2) # loop through dict, and add edges
            #################################
            # print out inter vs intra phyla counts
            #################################
            print "############################################################################"
            print "Intra phyla hits: %d" % (intraPhyla)
            print "Inter phyla hits: %d" % (interPhyla)
            print "Spirochaete"
            print "Intra: %d" % (SpirochaeteIntra)
            print "Inter: %d" % (SpirochaeteInter)
            print "Lentisphaerae"
            print "Intra: %d" % (LentisphaeraeIntra)
            print "Inter: %d" % (LentisphaeraeInter)
            print "Firmicutes"
            print "Intra: %d" % (FirmicutesIntra)
            print "Inter: %d" % (FirmicutesInter)
            print "Bacteroidetes"
            print "Intra: %d" % (BacteroidetesIntra)
            print "Inter: %d" % (BacteroidetesInter)
            print "Synergistetes"
            print "Intra: %d" % (SynergistetesIntra)
            print "Inter: %d" % (SynergistetesInter)
            print "Actinobacteria"
            print "Intra: %d" % (ActinobacteriaIntra)
            print "Inter: %d" % (ActinobacteriaInter)
            print "Tenericutes"
            print "Intra: %d" % (TenericutesIntra)
            print "Inter: %d" % (TenericutesInter)
            print "Fusobacteria"
            print "Intra: %d" % (FusobacteriaIntra)
            print "Inter: %d" % (FusobacteriaInter)
            print "Proteobacteria"
            print "Intra: %d" % (ProteobacteriaIntra)
            print "Inter: %d" % (ProteobacteriaInter)
            print "############################################################################"
            
        else: 
            print "                          Subsetting data to %d                                "% (self.subSet)
            print "*******************************************************************************"
            
            subSetIDs = []
            workingIdsHits = []
            # get random subset of nodes
            for id in self.HD.workingIds:
                if self.HD.workingIds[id]:
                    workingIdsHits.append(id)
                          
            subSetIDs = random.sample(workingIdsHits,self.subSet)
            
            G.add_nodes_from(subSetIDs)
            subSetIDsDict = {}
            phylumDict    = {}
            colours = self.ReturnColours(1, 100) # list of rgb colours
            count = 1 
            
            #####################################
            """colour by phylum"""
            #####################################
            
            #if self.colourBy == "sequencingPlatform":
            
            for id in subSetIDs:
                subSetIDsDict[id] = 0
                try: 
                    phylumCols[str(id)] = phylumDict[self.MD.metadata[str(id)][0]]
                except KeyError:
                    i = 1 * count 
                    #phylumDict[self.MD.metadata[str(id)][0]] = colours[i]
                    #phylumCols[str(id)] = colours[i]
                    #count += 10
                    if self.MD.metadata[str(id)][0] == "Spirochaetes":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#ea00ff"
                        phylumCols[str(id)] = "#ea00ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Lentisphaerae":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#ffb700" 
                        phylumCols[str(id)] = "#ffb700"
                        
                    elif self.MD.metadata[str(id)][0] == "Firmicutes":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#0047ff"
                        phylumCols[str(id)] = "#0047ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Bacteroidetes":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#14ff00"
                        phylumCols[str(id)] = "#14ff00"
                        
                    elif self.MD.metadata[str(id)][0] == "Synergistetes":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#6600CC" 
                        phylumCols[str(id)] = "#6600CC"
                        
                    elif self.MD.metadata[str(id)][0] == "Actinobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#ffff00" 
                        phylumCols[str(id)] = "#ffff00"
                        
                    elif self.MD.metadata[str(id)][0] == "Tenericutes":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#0080ff" 
                        phylumCols[str(id)] = "#0080ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Fusobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#00e0ff" 
                        phylumCols[str(id)] = "#00e0ff"
                        
                    elif self.MD.metadata[str(id)][0] == "Proteobacteria":
                        phylumDict[self.MD.metadata[str(id)][0]] = "#ff1e00"
                        phylumCols[str(id)] = "#ff1e00"
                    
            for phylum in phylumDict.keys():
                print "phylum: %s Colour: %s" % (phylum, phylumDict[phylum])
            
            print "There are %d different phylums represented" % (len(phylumDict))
            print "Length of Phylum Colour list %d" % len(phylumCols)

            values = [phylumCols.get(node,0.25) for node in G.nodes()]
            
            #####################################
            """colour by body site"""
            #####################################
            
            #elif self.colourBy == "sequencingPlatform":
                
            values = [phylumCols.get(node,0.25) for node in G.nodes()]
            #####################################
            """colour by sequencing platform"""
            #####################################
            
            #elif self.colourBy == "sequencingPlatform":
                
            values = [phylumCols.get(node,0.25) for node in G.nodes()] 
            #####################################
            """colour by sequencing centre"""
            #####################################
            
            #elif self.colourBy == "sequencingPlatform":


            values = [phylumCols.get(node,0.25) for node in G.nodes()]
            #####################################
            """Add edges between interacting genomes"""
            #####################################
            for id_1 in self.HD.hitCounts.keys():
                try: 
                    fred = subSetIDsDict[id_1] # check if in dictionary
                    for id_2 in self.HD.hitCounts[id_1]:
                        try: 
                            jimmy = subSetIDsDict[id_2] # check if in dictionary
                            if self.HD.hitCounts[id_1][id_2] >0:
                                if weighted:
                                    G.add_edge(id_1,
                                               id_2,
                                               weight=self.HD.hitCounts[id_1][id_2]) # loop through dict, and add edges
                                else:
                                    G.add_edge(id_1,id_2) # loop through dict, and add edges
                                "Edge added between %s and %s" % (id_1,id_2)
                        except KeyError:
                            pass
                except KeyError:
                    pass
                
        # edit edge properties # This needs to be edited to include the number of hits/genome size.
        #for edge in G.edges():
        #    edgeWidth.append(int(ids_dict[edge[0]][edge[1]][0]))
        
        # Set seed
        #random.seed(1)
        
        pos= nx.spring_layout(G,iterations=500)
        #pos= nx.circular_layout(G)
        fig = plt.figure(figsize=(21,10),dpi=300)
        
        
        #nx.draw(G,pos,node_color=values,with_labels=args.labels,width=edgeWidth,font_size=12,font_color='#006400')
        ###################
        # Colour by phylum
        ##################
        # Defaults:
        # alpha=0.4
        # node_size=10
        
        plt.subplot(2,2,1,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        #nx.draw(G,pos,node_color=phylumValues,with_labels=IPN.nodeLabels,width=edgeWidth,font_size=12,font_color='#006400')
        nx.draw_networkx_nodes(G,pos,node_color = phylumValues, linewidths=0, alpha=0.45, node_size=10)
        nx.draw_networkx_edges(G, pos, edge_color = "#373737",width=0.3) # 373737
        plt.tick_params(axis='both',
                        which='both',      
                        bottom='off',      
                        top='off',         
                        labelbottom='off',
                        labelleft='off')
        IPN.nodeLabels
        ###################
        # Colour by bodysite
        ################## 
        plt.subplot(2,2,2,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,pos,node_color = bodysiteValues, linewidths=0, alpha=0.45, node_size=10)
        nx.draw_networkx_edges(G, pos, edge_color = "#373737",width=0.3)
        plt.tick_params(axis='both',
                        which='both',      
                        bottom='off',      
                        top='off',         
                        labelbottom='off',
                        labelleft='off')
        
        ###################
        # Colour by seqPlat
        ##################
        plt.subplot(2,2,3,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,pos,node_color = seqPlatValues, linewidths=0, alpha=0.45 ,node_size= 10)
        nx.draw_networkx_edges(G, pos, edge_color = "#373737",width=0.3)
        plt.tick_params(axis='both',
                        which='both',      
                        bottom='off',      
                        top='off',         
                        labelbottom='off',
                        labelleft='off')
        
        ###################
        # Colour by seqCentre
        ##################
        plt.subplot(2,2,4,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,pos,node_color = seqCentreValues, linewidths=0, alpha=0.45, node_size=10)
        nx.draw_networkx_edges(G, pos, edge_color = "#373737",width=0.3)
        plt.tick_params(axis='both',
                        which='both',      
                        bottom='off',      
                        top='off',         
                        labelbottom='off',
                        labelleft='off')
        
        
        #nx.draw(G,
        #        pos,
        #        with_labels =IPN.nodeLabels,
        #        font_size       =IPN.nodeFontSize,
        #        font_color      =IPN.nodeFontColour, 
        #        node_color      = values, 
        #        node_size       = IPN.nodeSize,
        #        alpha           = 0.5, 
        #        edge_color      = "#E8E8E8"
        #        )

        
        if showPlot:
            plt.show()
        else:
            outputFile = IPN.outFile+".%s" % (IPN.imageFormat)
            print "Writing to file: %s" % (outputFile)
            plt.savefig(outputFile,dpi=IPN.dpi,format=IPN.imageFormat)

    #def plotNetwork(self,G,pos):
        
        

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

        total = 0 # total number of hits  
        # create total number of hits at each percentage
        for id1 in self.HD.hitCounts.keys():
            for id2 in self.HD.hitCounts[id1].keys():
                total += self.HD.hitCounts[id1][id2]
                try: 
                    roundedHitsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] += self.HD.hitCounts[id1][id2]
                except KeyError:
                    roundedHitsInt[math.ceil(float(self.HD.identityAni[id1][id2]))] = self.HD.hitCounts[id1][id2]
                try:
                    roundedHitsDec[round(float(self.HD.identityAni[id1][id2]),1)] += self.HD.hitCounts[id1][id2]
                except KeyError:
                    roundedHitsDec[round(float(self.HD.identityAni[id1][id2]),1)] = self.HD.hitCounts[id1][id2]
                    
        print "total no. of hits in hitObject %d" % total                    
                    
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
        for i in range(960):
            percentListDecY.append(0)
            percentHitsDecY.append(0)
        
        for i,v in enumerate(intList):
            try:
                percentListIntY[i] = roundedAllComparisonsInt[v]
            except KeyError:
                pass
        
        """
        for i,v in enumerate(decList):
            try:
                percentListDecY[i] = roundedAllComparisonsDec[v]
            except KeyError:
                pass
        
        
        # get y values for hits
        for i,perc in enumerate(decList):
            try:
                c = roundedAllComparisonsDec[i] / float(100)
                percentHitsDecY[i] = roundedHitsDec[i]/float(c)
                print perc, percentHitsDecY
                print perc, roundedAllComparisonsDec[perc], roundedHitsDec[i], percentHitsDecY[i]
            except KeyError:
                pass
        """
        
        # get y values for hits
        for i,perc in enumerate(intList):
            print i
            print perc
            try:
                c = roundedAllComparisonsInt[perc] / float(100)
                percentHitsIntY[perc] = roundedHitsInt[perc] / float(c)
                print perc, roundedAllComparisonsInt[perc], roundedHitsInt[perc], percentHitsIntY[perc]
            except KeyError:
                pass
        
        
        print "*******************************************************************************"
        print "                          Building frequency plot                              "
        print "*******************************************************************************"
        
        # plot 
        plt.scatter(intList, percentHitsIntY, marker='|')
        plt.plot(intList, percentHitsIntY, linestyle='-')
        
        #plt.scatter(decList, percentHitsDecY, marker='|')
        #plt.plot(decList, percentHitsDecY, linestyle='-')
        
        # plot no. of comparisons
        plt.plot(intList,percentListIntY,c='#FFFFFF')
        plt.fill_between(intList, percentListIntY, 1e-6, facecolor = '#C0C0C0')
        plt.axis([100,-5,0,max(percentHitsIntY)+10])
        
        #plt.plot(decList,percentListDecY,c='#FFFFFF')
        #plt.fill_between(decList, percentListDecY, 1e-6, facecolor = '#C0C0C0')
        #plt.axis([100,-5,0,max(percentHitsDecY)+10])
        
        #plt.xLabel(IPF.xLabel)
        #plt.yLabel(IPF.yLabel)
        plt.title(IPF.title)
        if IPF.showPlot:
            plt.show()
        else:
            outputFile = IPF.outFile+".%s" % (IPF.imageFormat)
            print "                          Writing to file: %s                              " % (outputFile)
            plt.savefig(outputFile,dpi=IPF.dpi,format=IPF.imageFormat)
        
        
        
        
        """
        
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
        """
    def ReturnColours(self,start,count):
        colour_array = []
        for i in range(int(start),int(start)+int(count)):
            colour_array.append('#%02x%02x%02x' % tuple(np.array(htr(float(i)/100.,1,1))*255))
        return colour_array


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    """visualise database"""
    def __init__(self,
                 dbFileName,
                 metadata,
                 subSet,
                 colourBy,
                 status,
                 cleanFastaFile,
                 ani=95,
                 batches=[],       # list of batch numbers to process
                 ):
        
        # get an interface to the DB
        VI = ViewInterface(dbFileName)
        
        # build the condition we want to select on and get the hit data
        C = Condition("ani_comp", "<=", ani)
        if len(batches) > 0:
            bc = Condition("batch", "=", batches[0])
            for batch_num in range(1, len(batches)):
                bc = Condition(bc, "or", Condition("batch", "=", batches[batch_num]))
            C = Condition(C, "and", bc)
        
        # get hitData object
        self.hits = VI.getHitData(C)
         
        # Build clean hit data object
        self.cleanHits = CleanHitData()
        self.cleanHits.getData(cleanFastaFile, self.hits)
        
        # call plotting functions
        self.PLOT = Plotter(self.hits, 
                            subSet,
                            metadata,
                            colourBy,
                            status,
                            self.cleanHits)

        
        
        
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
        pass
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
                    weighted,
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
                                  weighted,
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
