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
import zlib
from subprocess import Popen, PIPE

# local imports
from trackm.hit import Hit
from trackm.exceptions import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NucMerParser(object):
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
        
class ContigParser(object):
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    header = l.rstrip()[1:].partition(" ")[0] # save the header we just saw
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break

class genome_contigs(object):
    def __init__(self):
        self.genome_1_contigs = {} # dictionary to store first genome data
        self.genome_2_contigs = {} # dictionary to store second genome data

    def addContig(self,which_genome,generator):
        if which_genome == "genome1":
            contigs = generator
            for contig in contigs:
                contig_name = contig[0] 
                contig_seq  = contig[1]
                self.genome_1_contigs[contig_name] = contig_seq 
                
        if which_genome == "genome2":
            contigs = generator
            for contig in contigs:
                contig_name = contig[0] 
                contig_seq  = contig[1]
                self.genome_2_contigs[contig_name] =  contig_seq
            
    def returnContig(self,contig_name,which_genome,start,length):
        start = start - 1 
        stop = (start - 1) + length
        if which_genome == "genome1":
            return self.genome_1_contigs[contig_name][start:stop]
        if which_genome == "genome2":
            return self.genome_2_contigs[contig_name][start:stop]

class Worker(object):
    def __init__(self,
                 workID,            # workId for this task
                 gPath1,            # absolute path to the first genome
                 gid1,              # genome tree id for first genome
                 gPath2,            # absolute path to the second genome
                 gid2,              # genome tree id for second genome
                 ani,
                 serverURL,         # URL of the commanding TrackM server
                 ):
        self.gPath1 = gPath1
        self.gPath2 = gPath2
        self.gid1   = gid1 
        self.gid2   = gid2
        self.workID = workID
        self.serverURL = serverURL
        GC = genome_contigs()
        CP = ContigParser() 
        
        # this dictionary will store all the results of the analysis
        # essentially it will be a list of Hit instances
        self.results = [self.workID, int(ani*1000.)]
        
        # Capture genome1's contigs in dictionary
        with open(self.gPath1,'r') as fh:
            GC.addContig("genome1", CP.readFasta(fh))
            
        # Capture genome2's contigs in dictionary
        with open(self.gPath2,'r') as fh:
            GC.addContig("genome2", CP.readFasta(fh))
        
    def runCommand(self, cmd):
        """Run a command and take care of stdout

        expects 'cmd' to be a string like "foo -b ar"

        returns (stdout, stderr)
        """
        p = Popen(cmd.split(' '), stdout=PIPE)
        return p.communicate()

    def compareGenomes(self,
                       minLength=500,       # minimum length to be called a hit
                       minIdentity=99       # minimum identity to be called a hit
                       ):
        """Do the work of pairwise comparison of genomes"""
        try:
            # run pairwise nucmer instance
            # we expect this to be run from within a directory (handled by trackm.sge)
            # so there's no need to worry about file names
            output = self.runCommand("nucmer %s %s --mum --coords" % (self.gPath1, self.gPath2))
            self.results.append(output)

            # check to see that the 'out.coords' file exists and is non-empty
            ran_OK = True

            if ran_OK:
                # parse nucmer .coords file
                self.getHitData(minLength, minIdentity)
            else:
                self.results.append("ERROR")
                self.results.append(output)

        except Exception:
            # catch all exception, if anything goes wrong in the
            # comparison stage and we fail to catch it then we
            # can catch it here and return the exception to the
            # controlling server, that way it will know that the
            # job was aborted.
            print exc_info()[0]
            self.results.append("ERROR")
            self.results.append(exc_info()[0])
            raise
        
        # once all the comparisons are done (or have failed....)
        # invoke phoneHome to send results back to the calling server
        #self.phoneHome() ## REMOVE THE # TO RUN #### 


    def getHitData(self, minLength, minIdentity):
        """Filter Nucmer hits and add them to the result list"""
        NP = NucMerParser()
        GC = genome_contigs()
        with open('out.coords' % self.workID, 'r') as fh:
            for hit in NP.readNuc(fh):
                # apply filter >500bp and >99%
                if hit[NP._IDENTITY] >= minIdentity and hit[NP._LEN_1] >= minLength and hit[NP._LEN_2] > minLength:
                    # work out strandedness
                    if hit[NP._END_1] > hit[NP._START_1]:
                        # forward strand
                        strand1 = 0
                        start1 = hit[NP._START_1]
                    else:
                        strand1 = 1
                        start1 = hit[NP._END_1]

                    if hit[NP._END_2] > hit[NP._START_2]:
                        # forward strand
                        strand2 = 0
                        start2 = hit[NP._START_2]
                    else:
                        strand2 = 1
                        start2 = hit[NP._END_2]

                    # Get the seqs!
                    seq1 = GC.returnContig(hit[NP._ID_1], 'genome1', start1, hit[NP._LEN_1])
                    seq2 = GC.returnContig(hit[NP._ID_2], 'genome2', start2, hit[NP._LEN_2])

                    H = Hit(hit[NP._ID_1]+"_"+self.gid1,
                            start1,
                            hit[NP._LEN_1],
                            strand1,
                            seq1,
                            hit[NP._ID_2]+"_"+self.gid2,
                            start2,
                            hit[NP._LEN_2],
                            strand2,
                            seq2,
                            hit[NP._IDENTITY])

                    self.results.append(H)
        for hit in self.results:
            print "aasda",hit

    def phoneHome(self,
                  exception=None            # if there was some problem then this will not be None
                  ):
        # set up the socket
        context = zmq.Context()
        #print "W [%d] : Connecting to listener with port %s" % (self.id, self.port)
        socket = context.socket(zmq.REQ)
        socket.connect(self.serverURL)
        #print "W [%d] : Sending result" % self.id
        socket.send(jp.encode(self.results).encode("zlib"))
        message = socket.recv()
        #print "W [%d] : Received reply [ %s ]" % (self.id, message)
        if message == "DIE":
            #print "W [%d] : Ordered to die! Exiting now..." % self.id
            return


###############################################################################
###############################################################################
###############################################################################
###############################################################################
