#!/usr/bin/env python
###############################################################################
#                                                                             #
#    TrackM.py                                                         #
#                                                                             #
#    Description!!                                                            #
#                                                                             #
#    Copyright (C) Josh Daly                                                  #
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

__author__ = "Josh Daly"
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
import operator
import random

from multiprocessing import Pool
from subprocess import Popen, PIPE

import os
import numpy as np
np.seterr(all='raise')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

from colorsys import hsv_to_rgb as htr



###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TemplateClass():
    """Utilities wrapper"""
    def __init__(self): pass

    def sayHi(self):
        print('write some "REAL" code you bum!')
      
    def ReturnColours(self,start,count):
        colour_array = []
        for i in range(int(start),int(start)+int(count)):
            colour_array.append('#%02x%02x%02x' % tuple(np.array(htr(float(i)/100.,1,1))*255))
        return colour_array
        
        
    def doWork(self,args):
        """ Main wrapper"""
        
        """make scatter plot from a tab delimited file 
        img_id_a        contig_a        img_id_b        contig_b        hits    length
        """
        # variables
        ids_list = [] # list to store uniquq ids
        unordered_ids_list = [] # contained unordered img Ids 
        #sorted_ids_list = [] # contains IMG IDs ordered according to sorted_master
        genome_tree_ids_list = [] # Ordered list of genome tree IDs
        genome_tree_ids_dict = {} # {genome_tree_id:ladder_position}
        ids_dict = {} # dict to store {id:[[id,hits,length],...]}
        genome_tree_img_dict = {} # dict to store {id:tax_string}
        #master_tax = [] #list to store img_ids taxa information
        #sorted_master_tax = [] # Master sorted img IDs by tax string
        #img_sorted_dict = {} # contains img_id -> x coordinate
        working_ids_list = [] # sorted img_ids 
        plot_border = args.plot_border 
        ordered_tax_string_lowest = [] # store lowest taxonomical information in order
        ordered_tax_string_phylum = [] # store phylum-level information in order
        subs_ordered_tax_string_phylum = [] #subset of unique entries in ordered_tax_string_lowest
        ticks = [] # array of numbers, for producing ticks on figure
        unique_taxa = {} #dictionary to check for unique taxa for plotting rectangles
        colours = [] #list of colours produced from input colour file
        Lengths = []
        Hits = []
        
        #-----
        
        #change arguments into variables
        output_file = args.output_file + ".%s" % (args.image_format)
        image_format = args.image_format
       
        #-----
        
        
        """
        1/ create a dictionary containing img_ids -> ladder_position
        A00000033       637000039       123    taxa_string
        """ 
        with open(args.genome_tree_file,"r") as fh:
            for l in fh:
                genome_tree_id = l.split('\t')[0].rstrip()
                img_id = l.split('\t')[1].rstrip()
                ladder_position = l.split('\t')[2].rstrip()
                taxa_string = l.split('\t')[3].rstrip()
                genome_tree_img_dict[img_id] = [ladder_position,genome_tree_id,taxa_string] 
        #print genome_tree_img_dict
          
        """  
        with open(args.master_tax_file, "r") as fh:
            for l in fh:
                genome_tree_id = l.split('\t')[0].rstrip()
                img_id = l.split('\t')[1].rstrip()
                genome_tree_img_dict[genome_tree_id] = {img_id:genome_tree_ids_dict[genome_tree_id]}
                #tax_string = l.split('\t')[2].rstrip()
                #master_ids_list.append(img_id)
                master_tax.append(tax_string)
                tax_dict[img_id] = tax_string
            sorted_master_tax = sorted(tax_dict.iteritems(), key=operator.itemgetter(1))
        # make a dictionary containing img ID and order
        for i, v in enumerate(sorted_master_tax):
            img_sorted_dict[v[0]] = int(i)
        """
        #-----
        
        #read through file line by line
        """ genome_tree_id_a    img_id_a    bodysite_a    contig_a    genome_tree_id_b    img_id_b    bodysite_b    contig_b    hits    length"""
        
        with open(args.input_file, "r") as fh:
            # capture the file header
            header = fh.readline()
                           
            #read through the file line by line
            for l in fh:
                # assign each column in line to a variable
                genome_tree_a = l.split('\t')[0]
                id_a = l.split('\t')[1]
                body_site_a = l.split("\t")[2]
                contig_a = l.split("\t")[3]
                genome_tree_b = l.split('\t')[4]
                id_b = l.split('\t')[5]
                body_site_b = l.split("\t")[6]
                contig_b = l.split("\t")[7]
                hits = l.split('\t')[8].rstrip()
                length = l.split('\t')[9].rstrip()
                Hits.append(int(hits))
                Lengths.append(int(length))
                
                # check to see if key is in hash
                try:
                    ids_dict[id_a][id_b] = [int(hits),int(length),body_site_a,body_site_b,contig_a,contig_b]
                except KeyError:
                    ids_dict[id_a] = {id_b:[int(hits),int(length),body_site_a,body_site_b,contig_a,contig_b]}
                    #unordered_ids_list.append(img_sorted_dict[id_a])
                    
                try:                
                    ids_dict[id_b][id_a] = [int(hits),int(length),body_site_a,body_site_b,contig_b,contig_a]  
                except KeyError: 
                    ids_dict[id_b] = {id_a:[int(hits),int(length),body_site_a,body_site_b,contig_b,contig_a]}
                    #unordered_ids_list.append(img_sorted_dict[id_b])  
                    
            # make ids_list a numpy array
            ids_list= np.array(ids_dict.keys())
            for key in ids_dict.keys():
                unordered_ids_list.append(int(genome_tree_img_dict[key][0]))    
            # ordered img_ids  
            working_ids_list = ids_list[np.argsort(unordered_ids_list)]
        
        
            
        
        #-----
   
        #print out the tax string
        for i in working_ids_list:
            g = genome_tree_img_dict[i][2].split(';')[5] #genus
            f = genome_tree_img_dict[i][2].split(';')[4] #family
            o = genome_tree_img_dict[i][2].split(';')[3] #order
            c = genome_tree_img_dict[i][2].split(';')[2] #class
            p = genome_tree_img_dict[i][2].split(';')[1] #phylum
            
            ordered_tax_string_phylum.append(p)
            #print str(i) + "\t" +str(p) + "\t" + str(g) +"\t"+ str(genome_tree_img_dict[i][0])
            if len(g) == 4:
                if len(f) == 4:
                    if len(o) == 4:
                        if len(c) == 4:
                            if len(p) == 4:
                                print "error"
                            else:
                                #print p
                                ordered_tax_string_lowest.append(p)
                        else:
                            #print c
                            ordered_tax_string_lowest.append(c)
                    else:
                        #print o
                        ordered_tax_string_lowest.append(o) 
                else:
                    #print f
                    ordered_tax_string_lowest.append(f)   
            else:
                #print g
                ordered_tax_string_lowest.append(g)
        
        #print ordered_tax_string_lowest
        #-----
        #randomise list
        #random.shuffle(ids_list)
        
        #-----
        """ Create scatter plot data 
        """
        
        #create scatter plot variable, x and y
        xs = [] 
        ys = [] 
        val = []
        heatmap_colours = [] # list of colours to be used for the heatmap (scatterplot)
        #BuPu Hex colours
        BuPu = "#edf8fb" "#bfd3e6" "#9ebcda" "#8c96c6" "#8c6bb1" "#88419d" "#6e016b"
        #OrRd Hex colours
        OrRd = "#fef0d9"    "#fdd49e"    "#fdbb84"    "#fc8d59"  "#ef6548"  "#d7301f"    "#990000"
        #blue    #red    #green    #yellow
        "#2600FF"    "#FF0019"  "#51FF00"   "#FFFF00" 
        #Blue
        Blue = ["#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594"]
        #Green
        Green = ["#edf8e9", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32"]
        #Orange
        Orange = ["#feedde","#fdd0a2","#fdae6b","#fd8d3c","#f16913", "#d94801", "#8c2d04"]
        #Red
        Red = ["#fee5d9","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d"]
        #Black
        Black = ["#f7f7f7","#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525"]
        #purple
        Purple = ["#FFCCCC","#FF99FF","#FF66FF","#FF33FF","#FF00FF","#CC00CC","#990099"]
        
        # hits breaks
        hits_breaks = [0.01,0.1,0.2,0.3,0.4,0.5]
        #length breaks
        length_breaks = [0.01,0.05,0.1,0.15,0.25,0.5]
        
        # variables for heatmap
        hits_max = int(max(Hits))
        Length_max = int(max(Lengths))
        
        #-----
        
        # loop over dictionary, to produce x and y
        for x in range(len(working_ids_list)): 
            try:
                fred = ids_dict[working_ids_list[x]]
                for y in range(x+1, len(working_ids_list)):
                    try:
                        # Upper triangle - Hits data
                        val.append(fred[working_ids_list[y]][0])
                        xs.append(x)
                        ys.append(y)
                        
                        #-----
                        # colour based on body site
                        if fred[working_ids_list[y]][2] == "gut":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Red[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Red[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Red[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Red[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Red[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Red[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Red[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Blue[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Blue[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Blue[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Blue[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Blue[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Blue[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Blue[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                        if fred[working_ids_list[y]][2] == "oral":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Blue[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Blue[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Blue[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Blue[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Blue[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Blue[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Blue[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Green[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Green[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Green[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Green[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Green[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Green[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Green[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[0])
                        if fred[working_ids_list[y]][2] == "both":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[0]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[1]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[2]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[3]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[4]) and fred[working_ids_list[y]][0] <= (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][0] > (hits_max * hits_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                        
                        
                        """
                        if fred[working_ids_list[y]][2] == "gut":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#2600FF")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#FF0019")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                                
                        elif fred[working_ids_list[y]][2] == "oral":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#FF0019")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#51FF00")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                                
                        elif fred[working_ids_list[y]][2] == "both":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#FFFF00")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#FFFF00")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                             
                        
                        #-----
                        # scripts to make different coloured heatmaps for hits/length data
                       
                        # Add colour depending on value of hits
                        if fred[working_ids_list[y]][0] <= (hits_max * 0.1):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.1) and fred[working_ids_list[y]][0] <= (hits_max * 0.15):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.15) and fred[working_ids_list[y]][0] <= (hits_max * 0.2):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.2) and fred[working_ids_list[y]][0] <= (hits_max * 0.3):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.3) and fred[working_ids_list[y]][0] <= (hits_max * 0.4):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.4) and fred[working_ids_list[y]][0] <= (hits_max * 0.5):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][0] > (hits_max * 0.5):
                            heatmap_colours.append("#")
                        """
                        # Lower triangle - Length data
                        val.append(fred[working_ids_list[y]][1])
                        xs.append(y)
                        ys.append(x)
                        
                        
                        if fred[working_ids_list[y]][2] == "gut":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Red[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Red[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Red[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Red[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Red[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Red[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Red[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Blue[0])
                                elif fred[working_ids_list[y]][1] > (Length_max* length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Blue[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Blue[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Blue[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Blue[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Blue[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Blue[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][1] > (Length_max* length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max* length_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                        if fred[working_ids_list[y]][2] == "oral":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Blue[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Blue[1])
                                elif fred[working_ids_list[y]][1] > (Length_max* length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Blue[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Blue[3])
                                elif fred[working_ids_list[y]][1] > (Length_max* length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Blue[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Blue[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Blue[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Green[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Green[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Green[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Green[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Green[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Green[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Green[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                        if fred[working_ids_list[y]][2] == "both":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                                    
                            elif fred[working_ids_list[y]][3] == "oral":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max *length_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                                    
                            elif fred[working_ids_list[y]][3] == "both":
                                if fred[working_ids_list[y]][1] <= (Length_max * length_breaks[0]):
                                    heatmap_colours.append(Purple[0])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[0]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[1]):
                                    heatmap_colours.append(Purple[1])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[1]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[2]):
                                    heatmap_colours.append(Purple[2])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[2]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[3]):
                                    heatmap_colours.append(Purple[3])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[3]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[4]):
                                    heatmap_colours.append(Purple[4])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[4]) and fred[working_ids_list[y]][1] <= (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[5])
                                elif fred[working_ids_list[y]][1] > (Length_max * length_breaks[5]):
                                    heatmap_colours.append(Purple[6])
                        
                        
                        """
                        if fred[working_ids_list[y]][2] == "gut":
                            if fred[working_ids_list[y]][3] == "gut":
                                if fred[working_ids_list[y]][1] <= (length_max * 0.01):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.01) and fred[working_ids_list[y]][1] <= (length_max * 0.1):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.1) and fred[working_ids_list[y]][1] <= (length_max * 0.15):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.15) and fred[working_ids_list[y]][1] <= (length_max * 0.2):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.2) and fred[working_ids_list[y]][1] <= (length_max * 0.3):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.3) and fred[working_ids_list[y]][1] <= (length_max * 0.5):
                                    heatmap_colours.append("#")
                                elif fred[working_ids_list[y]][1] > (length_max * 0.5):
                                    heatmap_colours.append("#")
                        """
                        
                        
                        
                        
                        
                        """
                        # colour based on body site
                        if fred[working_ids_list[y]][2] == "gut":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#2600FF")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#FF0019")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                                
                        elif fred[working_ids_list[y]][2] == "oral":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#FF0019")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#51FF00")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                                
                        elif fred[working_ids_list[y]][2] == "both":
                            if fred[working_ids_list[y]][3] == "gut":
                                heatmap_colours.append("#FFFF00")
                            elif fred[working_ids_list[y]][3] == "oral":
                                heatmap_colours.append("#FFFF00")
                            elif fred[working_ids_list[y]][3] == "both":
                                heatmap_colours.append("#FFFF00")
                        
                        # Add color depending on value of length
                        if fred[working_ids_list[y]][1] <= (length_max * 0.01):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.01) and fred[working_ids_list[y]][1] <= (length_max * 0.1):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.1) and fred[working_ids_list[y]][1] <= (length_max * 0.15):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.15) and fred[working_ids_list[y]][1] <= (length_max * 0.2):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.2) and fred[working_ids_list[y]][1] <= (length_max * 0.3):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.3) and fred[working_ids_list[y]][1] <= (length_max * 0.5):
                            heatmap_colours.append("#")
                        elif fred[working_ids_list[y]][1] > (length_max * 0.5):
                            heatmap_colours.append("#")
                        """    
                    except KeyError:
                        pass
            except KeyError:
                pass
        
           
        #-----
        val = np.sqrt(val)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # paint a rectangle over plot
        currentAxis = plt.gca()
        
        #-----
        """ Taxonomy(phylum-level) grids """
        # Create dictionary with unique taxa strings, and ordered subset of unique taxa strings
        for i,v in enumerate(ordered_tax_string_phylum):
            if v in unique_taxa:
                unique_taxa[v][1] += 1
            else:
                unique_taxa[v] = [i,1]
                subs_ordered_tax_string_phylum.append(v)    
        #-----
        """
        # produce a random list of colours of len(working_ids_list)
        with open(args.colour_file,"r") as fh:
            for l in fh:
                colours.append(str(l.rstrip()))
        #randomise list of colours
        random.shuffle(colours)
        #print colours
        """
        
        #-----
        # produce rectangles for each unique taxa string
        for i in subs_ordered_tax_string_phylum:
            #y-axis rectangle
            currentAxis.add_patch(Rectangle((0-plot_border, unique_taxa[i][0]), len(working_ids_list)+plot_border-1, unique_taxa[i][1]-1, edgecolor = "grey",facecolor="none",alpha=0.1))
            #currentAxis.add_patch(Rectangle((0-plot_border, unique_taxa[i][0]), len(working_ids_list)+plot_border, unique_taxa[i][1], edgecolor = colours[v],facecolor="none" ,alpha=0.45))
            #x-axis rectangle
            currentAxis.add_patch(Rectangle((unique_taxa[i][0],0-plot_border), unique_taxa[i][1]-1, len(working_ids_list)+plot_border-1,edgecolor = "grey",facecolor="none",alpha=0.1))
            #currentAxis.add_patch(Rectangle((unique_taxa[i][0],0-plot_border), unique_taxa[i][1], len(working_ids_list)+plot_border,edgecolor = colours[v],facecolor="none",alpha=0.45))    
            
        #-----
        """ Produce scatter plot """
        ax.scatter(xs,
                   ys,
                   s=1.8,
                   marker='s',
                   alpha=1,
                   c=heatmap_colours,
                   edgecolors = 'grey',
                   linewidths = 0.1,
                   )
        
        #draw arbitrary line to differentiate top and bottom triangles
        plt.plot([0, len(working_ids_list)-1],[0, len(working_ids_list)-1],'k-',lw=0.5,alpha=0.3)
        
        
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
            labels[i] = (ordered_tax_string_lowest[i] + "_" + genome_tree_img_dict[working_ids_list[i]][1])
            #labels[i] = genome_tree_img_dict[working_ids_list[i]][1]
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
            

    def demoStuff(self):

        """
        # parse a file
        try:
            with open(filename, "r") as fh:
                for line in fh:
                    print line
        except:
            print "Error opening file:", filename, exc_info()[0]
            raise
        """

        """
        fig = plt.figure()

        #-----
        # make a 3d plot
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(points[:,0],
                   points[:,1],
                   points[:,2],
                   #edgecolors='none',
                   #c=colors,
                   #s=2,
                   #marker='.'
                   )

        #-----
        # make a 2d plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(points[:,0],
                points[:,1],
                '*g')

        #-----
        # show figure
        plt.show()
        # or save figure
        plt.savefig(filename,dpi=300,format='png')

        #-----
        # clean up!
        plt.close(fig)
        del fig
        """

        return 0
