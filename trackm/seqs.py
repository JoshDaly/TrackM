#!/usr/bin/env python
###############################################################################
#                                                                             #
#    seqs.py                                                                  #
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

# local imports
from trackm.viewInterface import ViewInterface, Condition


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Seqs(object):
    def __init__(self,
                 dbFileName,
                 ani=95,
                 batches=[]
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
            
        seqs = VI.getHitData(C)
        
    
    


