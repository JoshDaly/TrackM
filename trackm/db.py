#!/usr/bin/env python
###############################################################################
#                                                                             #
#    db.py                                                                    #
#                                                                             #
#    Implements the TrackMDB file format and corresponding API                #
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
__version__ = "0.1.0"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Dev"

__TRACKM_DB_VERSION__ = "1.0.0"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system imports

# local imports
from dancingPeasant.baseFile import BaseFile

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TrackMDB(BaseFile):
    def __init__(self, verbosity=0):
        BaseFile.__init__(self, verbosity)


    def createNewFile(self,
                      fileName,             # name of the new file
                      force=False,          # should we check to see if this is a wise move?
                      ):
        """Create a new TrackM database file"""
        # make a basic file
        BaseFile.createNewFile(self,
                               fileName,
                               type="TrackM_DB",
                               version=__TRACKM_DB_VERSION__,
                               force=force)

        # add TrackM specific tables
        self._addTable("ids",                               # housekeeping information about highest ids
                       {
                        "sqid" : "INT",                     # highest sqid in the seqs table (-1 for None)
                        "hid" : "INT",                      # highest hid in the hits table (-1 for None)
                        "cid" : "INT",                      # highest cid in the contigs table (-1 for None)
                        "oid" : "INT",                      # highest oid in the orfs table (-1 for None)
                       },
                       force=True)

        self._addTable("paths",                             # information about the location of genome files
                       {
                        "gid" : "TEXT",                     # genomeTree compatible ID for genome 1 (primary key)
                        "path" : "TEXT"                     # absolute path to fasta file
                       },
                       force=True)

        self._addTable("pairs",                             # information about comparisons to make
                       {
                        "pid" : "INT",                      # unique ID describing this pair (primary key)
                        "gid_1" : "TEXT",                   # genomeTree compatible ID for genome 1 (foreign key)
                        "gid_2" : "TEXT",                   # genomeTree compatible ID for genome 2 (gid_2 > gid_1) (foreign key)
                        "ani_1" : "REAL",                   # ANI from gid_1 -> gid_2
                        "ani_2" : "REAL",                   # ANI from gid_2 -> gid_1
                        "ident_tree" : "REAL",              # GTDB (83 marker) identity
                        "ident_16S" : "REAL",               # 16S pairwise identity
                        "batch" : "INT",                    # batch processed in
                        "ani_comp" : "REAL"                 # which ani has this pair been processed as?
                       },
                       force=True)

        self._addTable("hits",                              # information about hits between genome pairs
                       {
                        "hid" : "INT",                      # unique ID for this hit (primary key)
                        "pid" : "INT",                      # pair ID from pairs table (foreign key)
                        "cid_1" : "INT",                    # ID of the contig for this hit
                        "start_1" : "INT",                  # start position for this hit in gid_1
                        "len_1" : "INT",                    # length of this hit in gid_1
                        "strand_1" : "INT",                 # strand of the hit in gid_1
                        "sqid_1" : "INT",                   # sequence ID for this hit from seqs table (foreign key)
                        "cid_2" : "TEXT",
                        "start_2" : "INT",
                        "len_2" : "INT",
                        "strand_2" : "INT",
                        "sqid_2" : "INT",
                        "ident" : "REAL",                   # percent identity for this hit
                       },
                       force=True)

        self._addTable("contigs",                           # store contig headers once
                       {
                        "cid" : "INT",                      # unique ID for this contig (primary key)
                        "header" : "TEXT",                  # contig header
                       },
                       force=True)

        self._addTable("seqs",                              # information about orfs found on hits
                       {
                        "sqid" : "INT",                     # unique ID for this sequence (primary key)
                        "seq" : "TEXT",                     # nucelotide sequence
                       },
                       force=True)

        self._addTable("orfs",                              # information about orfs found on hits
                       {
                        "oid" : "INT",                      # unique ID for this orf (primary key)
                        "sqid" : "INT",                     # hit ID for this orf from in hits table (foreign key)
                        "start" : "INT",                    # start of the ORF
                        "len" : "INT",                      # length of the orf
                        "strand" : "INT",                   # strand of the orf
                        "seq_nucl" : "TEXT",                # nucelotide sequence
                        "seq_prot" : "TEXT"                 # protein sequence
                       },
                       force=True)

        self._addTable("ann_kegg",                          # kegg annotatations for orfs
                       {
                        "oid" : "INT",                      # orf id for this hit (foreign key)
                        "hit_1" : "TEXT",                   # best hit
                        "hit_1" : "TEXT",                   # second best hit
                        "hit_1" : "TEXT"                    # third best hit
                       },
                       force=True)

        self._addTable("ann_nog",                           # eggnog annotations for orfs
                       {
                        "oid" : "INT",                      # orf id for this hit (foreign key)
                        "hit_1" : "TEXT",                   # best hit
                        "hit_1" : "TEXT",                   # second best hit
                        "hit_1" : "TEXT"                    # third best hit
                       },
                       force=True)

        self._addTable("ann_nr",                            # nr annotations for orfs
                       {
                        "oid" : "INT",                      # orf id for this hit (foreign key)
                        "hit_1" : "TEXT",                   # best hit
                        "hit_1" : "TEXT",                   # second best hit
                        "hit_1" : "TEXT"                    # third best hit
                       },
                       force=True)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
