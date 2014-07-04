#!/usr/bin/env python
###############################################################################
#                                                                             #
#    worker.py                                                                #
#                                                                             #
#    Worker class designed to be run in Queue environment and communicate     #
#    with a TrackM server.                                                    #
#                                                                             #
#    The worker class does the heavy lifting of comparing pairwise genomes,   #
#    filtering results and collating information before returning the data    #
#    to the main server.                                                      #
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
from sys import exc_info
import os
import shutil

# local imports
from trackm.hit import Hit
from trackm.exceptions import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NucMerParser:
    """Wrapper class for parsing nucmer output"""
    # constants to make the code more readable
    _START_1  = 0
    _END_1    = 1
    _START_2  = 2
    _END_2    = 3
    _LEN_1    = 4
    _LEN_2    = 5
    _IDENTITY = 6
    _ID_1     = 7
    _ID_2     = 8

    def __init__(self):
        self.prepped = False

    def reset(self):
        self.prepped = False

    def readNuc(self, fp):
        """Read through a nucmer coords file

        this is a generator function
        """
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fp: # search for the first record
                        
                        if l[0] == '=': # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fp:
                fields = l.split('|')
                yield ([int(i) for i in fields[0].split()] +
                       [int(i) for i in fields[1].split()] +
                       [int(i) for i in fields[2].split()] +
                       [float(i) for i in fields[3].split()] +
                       fields[4].split())
            break # done!

class Worker(object):
    def __init__(self,
                 gPath1,            # absolute path to the first genome
                 gPath2,            # absolute path to the second genome
                 workID,            # workId for this task
                 serverURL,         # URL of the commanding TrackM server
                 serverPort         # port of the commanding TrackM server
                 ):
        self.gPath1 = gPath1
        self.gPath2 = gPath2
        self.workID = workID
        self.serverURL = serverURL
        self.serverPort = serverPort

        # this dictionary will store all the results of the analysis
        # essentially it will be a list of Hit instances
        self.results = {}

    def compareGenomes(self):
        """Do the work of pairwise comparison of genomes"""
        try:
            # this is the jumping off point for the comparison pipeline
            # already established in other TrackM scripts. Code inserted here should
            # resemble a pipeline that calls other functions defined in this class
            
            # make the tmp directory using workID
            os.mkdir("/tmp/trackm/%s" % self.workID)
            
            # move into tmp space 
            os.chdir("/tmp/trackm/%s" % self.workID)
            
            # run pairwise nucmer instance
            os.system("nucmer %s %s --mum --coords -p %s" % (self.gPath1,self.gPath2,self.workID))
            
            # parse nucmer .coords file
            self.getHitData()
            
            # Nuke tmp output directory
            #shutil.rmtree("/tmp/trackm/%s" % self.workID)
            
            # print out results
            for hit in self.results:
                print hit

            # once all the comparisons are done invoke phoneHome to send results back to the server
            self.phoneHome()
        except:
            # catch all exception, if anything goes wrong in the
            # comparison stage and we fail to catch it then we
            # can catch it here and return the exception to the
            # controlling server, that way it will know that the
            # job was aborted.
            print exc_info()[0]
            self.phoneHome(exception=exc_info()[0])
            
    def getHitData(self):
        # read in nucmer coords file
        NP = NucMerParser()
        with open('/tmp/trackm/%s' % self.workID, 'r') as fh:
            for hit in NP.readNuc(fh):
                # apply filter >500bp and >99% 
                #if hit[NP._IDENTITY] > 99 and hit[NP._LEN_1] > 500 and hit[NP._LEN_2] > 500:
                try:
                    self.results[self.workID] += [hit[NP._START_1],
                                                  hit[NP._END_1],
                                                  hit[NP._START_2].
                                                  hit[NP._END_2],
                                                  hit[NP._LEN_1],
                                                  hit[NP._LEN_2],
                                                  hit[NP._IDENTITY]]
                except KeyError:
                    self.results[self.workID] = [hit[NP._START_1],
                                                  hit[NP._END_1],
                                                  hit[NP._START_2].
                                                  hit[NP._END_2],
                                                  hit[NP._LEN_1],
                                                  hit[NP._LEN_2],
                                                  hit[NP._IDENTITY]]

    def phoneHome(self,
                  exception=None):
        """Send the results to the TrackM server

        If exception is None then we assume that it was a success!
        """
        pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################
