#!/usr/bin/env python
###############################################################################
#                                                                             #
#    viewInterface.py                                                         #
#                                                                             #
#    Read only Interface for accessing data in a TrackM DB                    #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
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

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system includes
import sys
import numpy as np

# local includes
from dancingPeasant.interface import Interface
from trackm.db import TrackMDB

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ViewInterface(Interface):
    """Use this interface for querying the TrackM DB"""

    def __init__(self,
                 dbFileName,        # file name to connect to
                 verbosity=-1       # turn off all DP chatter
                 ):
        Interface.__init__(self, dbFileName, verbosity)
        self.db = TrackMDB(verbosity=verbosity)

    def getOutstandingPairs(self, aniCut=95., batch=None):
        """Get all the outstanding pairs based on batch and ani cutoff

        returns a list of tuples of type:
        (job_id, gPath1, gPath2, gid1, gid2, batch, ani_comp)
        """
        self.connect()
        # get all the outstanding pairs
        if batch is None:
            all_pairs = self.select('pairs', ["pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch"], "ani_comp = '-1'")
        else:
            all_pairs = self.select('pairs', ["pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch"], "ani_comp = '-1' and batch = '%s'" % batch)
        gids = self.select('paths', ['gid', 'path'])
        self.disconnect()
        gids = dict(gids)
        vetted_pairs = []
        seen_gids = {}
        for pair in all_pairs:
            highest_ani = np.max([pair[3], pair[4]])
            if highest_ani <= aniCut:
                # pair is OK
                vetted_pairs.append((pair[0], gids[pair[1]], gids[pair[2]], pair[1], pair[2], pair[5], highest_ani))
        return vetted_pairs

    def getGids(self):
        """Get a list of all the GIDs"""
        self.connect()
        gids = self.select('paths', ['gid', 'path'])
        self.disconnect()
        # unwrap these guys
        self.gids = dict(zip([gid[0] for gid in gids], [True] * len(gids)))
        return self.gids

    def getGenomePaths(self):
        """Get a dict of type {gid -> path}"""
        self.connect()
        paths = self.select('paths', ["gid", "path"])
        self.disconnect()
        if len(paths) == 0:
            return {}
        return dict(paths)

    def getContigHeaders(self):
        """Get a dict of type {'header' -> cid}"""
        self.connect()
        headers = self.select('contigs', ['header', 'cid'])
        self.disconnect()
        if len(headers) == 0:
            return {}
        return dict(headers)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

