#!/usr/bin/env python
###############################################################################
#                                                                             #
#    importInterface.py                                                       #
#                                                                             #
#    Interface for importind data into a TrackM DB                            #
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system includes
import sys

# local includes
from dancingPeasant.interface import Interface
from trackm.db import TrackMDB
import csv

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ImportInterface(Interface):
    """Use this interface for importing data stored in csv
    files into the TrackM DB"""

    def __init__(self,
                 dbFileName,        # file name to connect to
                 verbosity=-1       # turn off all DP chatter
                 ):
        Interface.__init__(self, dbFileName, verbosity)
        self.db = TrackMDB(verbosity=verbosity)
        self.gids = None

    def getGids(self):
        """Get a list of all the GIDs"""
        self.connect()
        gids = self.select('paths', ['gid'])
        self.disconnect()
        # unwrap these guys
        self.gids = dict(zip([gid[0] for gid in gids], [True] * len(gids)))
        return self.gids

    def importGenomes(self, genomesFileHandle):
        """import new genomes into the TrackM DB"""
        if not self.connect(createDB=True):
            # database exists
            gids = self.select('paths', ['gid'])
            self.gids = dict(zip([gid[0] for gid in gids], [True] * len(gids)))
        else:
            self.gids = {}

        # skip comment lines
        dr = csv.DictReader((row for row in genomesFileHandle if not row.startswith('#')), delimiter="\t")
        dr.fieldnames = "gid", "path"

        to_db = [(row['gid'], row['path'])
                 for row in dr
                 if row["gid"] not in self.gids]    # not going to add the same GID twice

        if len(to_db) > 0:
            self.insert("paths", ["gid", "path"], to_db)
        self.disconnect()

    def importJobs(self, jobsFileHandle):
        """Import new jobs into the TrackM DB"""
        # get existing genomes
        if self.gids is None:
            self.getGids()
        elif self.gids == {}:
            self.getGids()

        self.connect()
        # skip comment lines
        dr = csv.DictReader((row for row in jobsFileHandle if not row.startswith('#')), delimiter="\t")
        dr.fieldnames = "pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch"

        # we will set processed to 0
        to_db = [(row["pid"], row["gid_1"], row["gid_2"], row["ani_1"], row["ani_2"], row["batch"], 0)
                 for row in dr]

        # work out now if there are any pairs which have no reference in the paths table
        to_db_b = []
        for row in to_db:
            if row[1] in self.gids and row[2] in self.gids:
                to_db_b.append(row)
            else:
                sys.stderr.write("Warning: unknown GIDS in job %s\n" %str(row))

        gid_pairs = {}
        to_db = []
        existing_jobs = self.select('jobs', ['gid_1', 'gid_2'])

        if len(existing_jobs) > 0:
            for job in existing_jobs:
                print job
                gid_pairs["%s|%s" % (job[0], job[1])] = True
            print gid_pairs
            for row in to_db_b:
                this_pair = "%s|%s" % (row[1], row[2])
                if this_pair not in gid_pairs:
                    to_db.append(row)
                else:
                    sys.stderr.write("Warning: Duplicate job: (Compare: %s to %s)\n" % tuple(this_pair.split("|")))
        else:
            to_db = to_db_b

        if len(to_db) > 0:
            self.insert("jobs", ["pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch", "processed"], to_db)
        self.disconnect()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

