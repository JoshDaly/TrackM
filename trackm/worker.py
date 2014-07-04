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

# local imports
from trackm.hit import Hit
from trackm.exceptions import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
        self.results = []

    def compareGenomes(self):
        """Do the work of pairwise comparison of genomes"""
        try:
            # this is the jumping off point for the comparison pipeline
            # already established in other TrackM scripts. Code inserted here should
            # resemble a pipeline that calls other functions defined in this class
            
            subgrp = 1000

            # make sure the out dir exists
            if not os.path.exists(args.out_dir):
                os.makedirs(args.out_dir)    
            
            """Read through tab-delimited file containing pairwise comparisons to be run"""
            pool = Pool(args.num_threads)
            count = 0 # counter for checkpoint of runtime
            cmds = [[]]
            stdouts = []
        
            with open(args.genome_comparisons) as fh: 
                """ ID_1         ID_2         %IDENTITY
                    2500069000    2501799900    93.58
                """
                header = fh.readline()
                counter = 0 
                for l in fh:
                    fields = l.rstrip().split("\t")
                    #cmds[-1].append("nucmer %s.fna %s.fna --mum --coords -p %s.gen > /%s" % (fields[0], fields[1], "%sv%s" %(fields[0],fields[1]), args.out_dir))
                    fasta1 = os.path.join(args.fasta_dir, "%s.%s" % (fields[0], args.fasta_suffix))
                    fasta2 = os.path.join(args.fasta_dir, "%s.%s" % (fields[1], args.fasta_suffix))
                    nucmer_prefix = os.path.join(args.out_dir, "%sv%s.gen" %(fields[0],fields[1]))
                    
                    cmds[-1].append("nucmer %s %s --mum --coords -p %s" % (fasta1,
                                                                           fasta2,
                                                                           nucmer_prefix)
                                    )
                    #cmds[-1].append("echo %s" %("%s %s %s" %(fasta1, fasta2, nucmer_prefix)))
                    counter += 1
                    if counter >= subgrp:
                        cmds.append([])
                        counter = 0
            
            print "start", datetime.datetime.now()
            for sub_cmds in cmds:
                stdouts.append(pool.map(runCommand, sub_cmds))            # list of tuples [(stdout, stderr)]
                print "%d done" % subgrp, datetime.datetime.now()
            
            print "finish", datetime.datetime.now()
            
            print "writing stdouts"
            for (out, err) in stdouts[0]:
                err_file = "%s.txt" % err.split('/')[2].split('.')[0]
                with open(os.path.join(args.out_dir, err_file), 'w') as err_fh:
                    for line in err:
                        err_fh.write(line)

            # >>>>>>>>>> REMOVE THIS WHEN YOU HAVE CODE HERE <<<<<<<<<<<<<<<<<<
            from inspect import currentframe, getframeinfo
            frameinfo = getframeinfo(currentframe())
            print "Insert comparison code in here!! I live at File: %s Line: %s" % (frameinfo.filename, frameinfo.lineno)
            # >>>>>>>>>> END <<<<<<<<<<<<<<<<<<

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
