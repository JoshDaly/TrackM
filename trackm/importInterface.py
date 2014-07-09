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

# local includes
from dancingPeasant.exceptions import *
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
            # we've jst created this database so set all the meta counts to -1
            self.insert('ids', ['sqid', 'hid', 'cid', 'oid'], [tuple([-1]*4)])

        # skip comment lines
        dr = csv.DictReader((row for row in genomesFileHandle if not row.startswith('#')), delimiter="\t")
        dr.fieldnames = "gid", "path"

        to_db = [(row['gid'], row['path'])
                 for row in dr
                 if row["gid"] not in self.gids]    # not going to add the same GID twice

        if len(to_db) > 0:
            self.insert("paths", ["gid", "path"], to_db)
        self.disconnect()

    def importPairs(self, pairsFileHandle):
        """Import new pairs into the TrackM DB"""
        # get existing genomes
        if self.gids is None:
            self.getGids()
        elif self.gids == {}:
            self.getGids()

        self.connect()
        # skip comment lines
        dr = csv.DictReader((row for row in pairsFileHandle if not row.startswith('#')), delimiter="\t")
        dr.fieldnames = "pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch"

        # we will set ani to -1
        to_db = [(row["pid"], row["gid_1"], row["gid_2"], row["ani_1"], row["ani_2"], row["batch"], -1)
                 for row in dr]

        # work out now if there are any pairs which have no reference in the paths table
        to_db_b = []
        for row in to_db:
            if row[1] in self.gids and row[2] in self.gids:
                to_db_b.append(row)
            else:
                sys.stderr.write("Warning: unknown GIDS in pair %s\n" %str(row))

        gid_pairs = {}
        to_db = []
        existing_pairs = self.select('pairs', ['gid_1', 'gid_2'])

        if len(existing_pairs) > 0:
            for pair in existing_pairs:
                print pair
                gid_pairs["%s|%s" % (pair[0], pair[1])] = True
            print gid_pairs
            for row in to_db_b:
                this_pair = "%s|%s" % (row[1], row[2])
                if this_pair not in gid_pairs:
                    to_db.append(row)
                else:
                    sys.stderr.write("Warning: Duplicate pair: (Compare: %s to %s)\n" % tuple(this_pair.split("|")))
        else:
            to_db = to_db_b

        if len(to_db) > 0:
            self.insert("pairs", ["pid", "gid_1", "gid_2", "ani_1", "ani_2", "batch", "ani_comp"], to_db)
        self.disconnect()

    def importHits(self,
                   contigHeaders,   # dict of type {header -> cid}
                   hits,            # list of hit and other information returned from TrackM worker
                   ):
        """Add the new hits to the DB

        hits is a list structured like this:

        [pid, ani_comp, <trackm.hit.Hit object>, <trackm.hit.Hit object>, ...]
        """
        self.connect()
        sqid, hid, cid, oid = self.getHighestIds()
        pid = hits[0]
        ani_comp = hits[1]
        to_db = []
        new_contigs = []
        new_seqs = []
        for i in range(2, len(hits)):
            # for each new hit

            # work out if we'e seen the contigs before or make a new entry
            if hits[i].contig1 == "repeated":
                print 'REPEAT!'
            try:
                cid1 = contigHeaders[hits[i].contig1]
            except KeyError:
                # new contig!
                cid += 1
                cid1 = cid
                contigHeaders[hits[i].contig1] = cid1
                new_contigs.append((cid, "%s" % hits[i].contig1))

            try:
                cid2 = contigHeaders[hits[i].contig2]
            except KeyError:
                # new contig!
                cid += 1
                cid2 = cid
                contigHeaders[hits[i].contig2] = cid2
                new_contigs.append((cid, "%s" % hits[i].contig2))

            # we know the seqs are new
            sqid += 1
            sqid1 = sqid
            new_seqs.append((sqid1, hits[i].seq1))
            sqid += 1
            sqid2 = sqid
            new_seqs.append((sqid2, hits[i].seq2))

            # new hid and add!
            hid += 1

            # make a tuple for the add
            to_db.append(tuple([hid,
                                pid,
                                cid1,
                                hits[i].start1,
                                hits[i].len1,
                                hits[i].strand1,
                                sqid1,
                                cid2,
                                hits[i].start2,
                                hits[i].len2,
                                hits[i].strand2,
                                sqid2,
                                hits[i].identity]
                               )
                         )

        # insert the hits
        self.insert("hits", ["hid", "pid", "cid_1", "start_1", "len_1", "strand_1", "sqid_1", "cid_2", "start_2", "len_2", "strand_2", "sqid_2", "ident"], to_db)

        # insert the ani_comp into the pairs table to say that this pair has been processed
        self.updateSingle("pairs", ["ani_comp"], ["%0.2f"% ani_comp], condition="pid='%d'"%pid)

        # insert the new contigs
        self.insert('contigs', ["cid", "header"], new_contigs)

        # insert the new sequences
        self.insert('seqs', ["sqid", "seq"], new_seqs)

        # update Id counts
        self.updateIds(sqid, hid, cid, oid)

        self.disconnect()

#------------------------------------------------------------------------------
# Handling IDs

    def updateIds(self, sqid, hid, cid, oid):
        """Update the values of the highest Ids for certain primary keys"""
        disconnect = True
        try:
            self.connect()
        except DP_FileAlreadyOpenException:
            # sometimes called from within a connect block
            disconnect = False

        self.updateSingle('ids', ['sqid', 'hid', 'cid', 'oid'], tuple([sqid, hid, cid, oid]))

        if disconnect:
            self.disconnect()

    def getHighestIds(self):
        """Get the values of the highest Ids for certain primary keys"""
        disconnect = True
        try:
            self.connect()
        except DP_FileAlreadyOpenException:
            # sometimes called from within a connect block
            disconnect = False

        ids = self.select('ids', ['sqid', 'hid', 'cid', 'oid'])

        if disconnect:
            self.disconnect()

        return ids[0]


###############################################################################
###############################################################################
###############################################################################
###############################################################################

