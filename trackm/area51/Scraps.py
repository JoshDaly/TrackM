#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = ""
__status__ = "Development"

###############################################################################

import argparse
import sys

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    #global variables
    genome_tree_ids_list = []
    genome_tree_ids_dict = {}
    genome_tree_img_dict = {}
    
    with open(args.genome_tree_file, "r") as fh:
        for l in fh:
            genome_tree_ids_list.append(l.rstrip()) #ordered list of genome tree IDs
    for i, v in enumerate(genome_tree_ids_list):
        genome_tree_ids_dict[v] = int(i) #{genome_tree_id:ladder_position}
    
    with open(args.master_tax_file, "r") as fh:
        for l in fh:
            genome_tree_id = l.split('\t')[0].rstrip()
            img_id = l.split('\t')[1].rstrip()
            tax_info = l.split('\t')[2].rstrip()
            try:
                genome_tree_img_dict[genome_tree_id] = [img_id,genome_tree_ids_dict[genome_tree_id],tax_info]
            except KeyError:
                pass

    for key in genome_tree_img_dict:
            print "\t".join([str(key),
                            str(genome_tree_img_dict[key][0]),
                            str(genome_tree_img_dict[key][1]),
                            str(genome_tree_img_dict[key][2])
                            ])
            
    #-----
    #print header
        print "\t".join(["genome_tree_id_a",
                         "img_id_a",
                         "bodysite_a",
                         "contig_a",
                         "genome_tree_id_b",
                         "img_id_b",
                         "bodysite_b",
                         "contig_b",
                         "hits",
                         "length"
                         ])
        
        for x in range(len(working_ids_list)):
            try:
                jimmy = ids_dict[working_ids_list[x]] # [id_b:[],id_c:[]]...
                for y in range(x+1, len(working_ids_list)):
                    try:
                        print "\t".join([str(genome_tree_img_dict[working_ids_list[x]][1]),
                                         str(working_ids_list[x]),
                                         str(jimmy[working_ids_list[y]][2]),
                                         str(jimmy[working_ids_list[y]][4]),
                                         str(genome_tree_img_dict[working_ids_list[y]][1]),
                                         str(working_ids_list[y]),
                                         str(jimmy[working_ids_list[y]][3]),
                                         str(jimmy[working_ids_list[y]][5]),
                                         str(jimmy[working_ids_list[y]][0]),
                                         str(jimmy[working_ids_list[y]][1])
                                         ])
                    except KeyError:
                        pass
            except KeyError:
                pass         
    

    """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""

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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_tree_file', help="Input file")
    parser.add_argument('master_tax_file', help="Input file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
