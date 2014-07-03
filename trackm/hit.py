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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Hit(object):
    """Simple storage class"""
    def __init__(self,
                 start1,
                 len1,
                 strand1,
                 seq1,
                 start2,
                 len2,
                 strand2,
                 seq2,
                 identity
                 ):
         self.start1 = start1
         self.len1 = start2
         self.strand1 = strand1
         self.seq1 = seq1
         self.start2 = start2
         self.len2 = len2
         self.strand2 = strand2
         self.seq2 = seq2
         self.identity = identity

    def __str__(self):
        return "\t".join([str(i) for i in [self.start1,
                                           self.len1,
                                           self.strand1,
                                           self.seq1,
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
