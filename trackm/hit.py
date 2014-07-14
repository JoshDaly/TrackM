#!/usr/bin/env python
###############################################################################
#                                                                             #
#    hit.py                                                                   #
#                                                                             #
#    Class for storing information about a single pairwise hit                #
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

import random
import string

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Hit(object):
    """Simple storage class"""
    def __init__(self,
                 contig1=None,
                 start1=None,
                 len1=None,
                 strand1=None,      # == 0 forward, == 1 reverse
                 seq1=None,
                 contig2=None,
                 start2=None,
                 len2=None,
                 strand2=None,
                 seq2=None,
                 identity=None
                 ):
         self.contig1 = contig1
         self.start1 = start1
         self.len1 = len1
         self.strand1 = strand1
         self.seq1 = seq1
         self.contig2 = contig2
         self.start2 = start2
         self.len2 = len2
         self.strand2 = strand2
         self.seq2 = seq2
         self.identity = identity

    def createRandom(self):
        """fill this hit with random values"""
        self.start1 = random.randint(0,5000)
        self.len1 = random.randint(0,200)
        self.strand1 = random.randint(0,1)
        self.seq1 = "".join([['A','C','G','T'][random.randint(0,3)] for _ in range(self.len1)])

        self.contig2 = "".join([string.ascii_letters[random.randint(0,51)] for _ in range(100)])
        self.start2 = random.randint(0,5000)
        self.len2 = random.randint(0,200)
        self.strand2 = random.randint(0,1)
        self.seq2 = "".join([['A','C','G','T'][random.randint(0,3)] for _ in range(self.len2)])

        self.identity = float(random.randint(0,100))/100.

        # just to test contig header storage
        if random.randint(0,99) > 97:
            self.contig1 = "repeated"
        else:
            self.contig1 = "".join([string.ascii_letters[random.randint(0,51)] for _ in range(100)])


    def __str__(self):
        return "\t".join([str(i) for i in [self.contig1,
                                           self.start1,
                                           self.len1,
                                           self.strand1,
                                           self.seq1,
                                           self.contig2,
                                           self.start2,
                                           self.len2,
                                           self.strand2,
                                           self.seq2,
                                           self.identity
                                           ]
                          ]
                         )

###############################################################################
###############################################################################
###############################################################################
###############################################################################
