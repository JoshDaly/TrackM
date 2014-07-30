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
    """Class for encoding an alignment found between two contigs

    The hit is recorded in terms of its start position, length and strand with
    to how the contig was written in the file passed to the alignment algorithm.

    strand == 0 implies that the hit is on the contig as written
    strand == 1 implies it is on the reverse strand

    start and len always refers to the contig as written

    seq is the sequence that hit, if strand == 1 then it is the reverse complement
    the substring of length len starting at start from the contig as written.
    """
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

class HitData(object):
    """class for capturing hit data"""
    def __init__(self):
        self.hits               = {} # dict to store hit data
        self.length             = {} # dict to store length data
        self.distance16S        = {} # dict to store 16S distance
        self.roundedDistance    = {} # dict to store rounded 16S distance and total hits per 100 comparisons
        self.standardDeviation  = {} # dict to store standard deviation at each percentage
    
    def addLen(self,
               _ID_1,
               _ID_2,
               _LGT_LEN,
               ):   
        """create cumulative contig length object"""
        try: 
            self.length[_ID_1][_ID_2] += _LGT_LEN
        except KeyError:
            try:
                self.length[_ID_1][_ID_2]  = _LGT_LEN
            except KeyError:
                self.length[_ID_1] = {_ID_2 : _LGT_LEN}
        
        try: 
            self.length[_ID_2][_ID_1] += _LGT_LEN
        except KeyError:
            try:
                self.length[_ID_2][_ID_1] = _LGT_LEN
            except KeyError:
                self.length[_ID_2] = {_ID_1 : _LGT_LEN}
                    
    def add16SDist(self,
                   _ID_1,
                   _ID_2,
                   _PERC_ID
                   ): 
        """create 16S distance object, round up %"""
        try: 
            self.distance16S[_ID_1][_ID_2] = math.ceil(_PERC_ID)
        except KeyError:
            self.distance16S[_ID_1] = {_ID_2 : math.ceil(_PERC_ID)}
        try: 
            self.distance16S[_ID_2][_ID_1] = math.ceil(_PERC_ID)
        except KeyError:
            self.distance16S[_ID_2] = {_ID_1 : math.ceil(_PERC_ID)}
                
    def addHit(self,
               _ID_1,
               _ID_2,
               ):
        """add an LGT instance to the data store"""
        try: 
            self.hits[_ID_1][_ID_2] += 1
        except KeyError:
            try:
                self.hits[_ID_1][_ID_2] = 1
            except KeyError:
                self.hits[_ID_1] = {_ID_2 : 1}  
            
    def getIDS(self):
        """return array of IDs"""
        return self.hits.keys()

    def groupBy16S(self):
        """calculate the number of transfers per rounded 16S score"""
        for id_1 in self.hits.keys():
            for id_2 in self.hits[id_1]:
                try:
                    #self.roundedDistance[self.distance[id_1][id_2]] += [self.hits[id_1][id_2]]
                    self.roundedDistance[self.distance16S[id_1][id_2]] += self.hits[id_1][id_2]     # total hits per percentage
                    #self.roundedDistance[self.distance[id_1][id_2]] += 1     # total hits per percentage
                    self.standardDeviation[self.distance16S[id_1][id_2]] += [self.hits[id_1][id_2]] # hit array per percentage
                except KeyError:
                    #self.roundedDistance[self.distance[id_1][id_2]] = [self.hits[id_1][id_2]]
                    self.roundedDistance[self.distance16S[id_1][id_2]] = self.hits[id_1][id_2]     # total hits per percentage
                    #self.roundedDistance[self.distance[id_1][id_2]] = 1     # total hits per percentage
                    self.standardDeviation[self.distance16S[id_1][id_2]] = [self.hits[id_1][id_2]] # hit array per percentage
                  
    def calculateStD(self,perc):
        """return the standard deviation for each percentage"""
        hitList = []
        hitList = self.standardDeviation[perc]
        print perc,hitList
        stdev = np.array(hitList)
        #print  str(np.std(stdev, dtype=np.float64))
        return  np.std(stdev, dtype=np.float64) 
                    
    def numHits16S(self):
        """print stats about the number of hits in the roundedDistance dict"""
        for perc in self.roundedDistance.keys():
            totalPerPercent = 0
            length = len(self.roundedDistance[perc])
            for i in self.roundedDistance[perc]:
                totalPerPercent = totalPerPercent + i
            print "********************************************************************"
            print "%i" % (perc)
            print "Total hits: %i" % (totalPerPercent)
            print "Number of comparisons: %i" % (length)
            print "Average hits: %f" % (totalPerPercent/float(length))
            print "Min number of hits: %i" % (min(self.roundedDistance[perc]))
            print "Max number of hits: %i" % (max(self.roundedDistance[perc]))
            print "********************************************************************"
            
    def return16SPercArray(self):
        """return array of rounded 16S percentages """
        return self.roundedDistance.keys()
        
    def returnNormValues(self):
        """return array of hits normalised by 100 comparisons per rounded 16S score"""
        return 

###############################################################################
###############################################################################
###############################################################################
###############################################################################
