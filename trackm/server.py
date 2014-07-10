#!/usr/bin/env python
###############################################################################
#                                                                             #
#    server.py                                                                #
#                                                                             #
#    Implements the TrackM server                                             #
#                                                                             #
#   The server wraps the DB and responds to communications over network       #
#   ports. The server directs workers to do pairwise calculations on genomes  #
#   as well an being the base interface to the DB for TrackM view operations. #
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# system imports
import multiprocessing
import zmq
import json
import time
import random
import jsonpickle as jp
import zlib
import socket

# local imports
from trackm.importInterface import ImportInterface
from trackm.viewInterface import ViewInterface
from trackm.hit import Hit
from trackm.sge import SGE

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TestProcessWorker(object):
    def __init__(self,
                 port,
                 id,
                 ani):
        # set up the worker
        self.port = port
        self.id = id
        self.result = [self.id, int(ani*1000.)]      # HACK: jsonpickle will convert some floats to None
        #print "W [%d] : worker made" % self.id

    def start(self):
        """start the worker"""
        # do the work we came here to do
        for i in range(random.randint(2,6)):
            H = Hit()
            H.createRandom()
            self.result.append(H)
        self.phoneHome()

    def phoneHome(self):
        # set up the socket
        context = zmq.Context()
        #print "W [%d] : Connecting to listener with port %s" % (self.id, self.port)
        socket = context.socket(zmq.REQ)
        socket.connect ("tcp://localhost:%s" % self.port)
        #print "W [%d] : Sending result" % self.id
        socket.send(jp.encode(self.result).encode("zlib"))
        message = socket.recv()
        #print "W [%d] : Received reply [ %s ]" % (self.id, message)
        if message == "DIE":
            #print "W [%d] : Ordered to die! Exiting now..." % self.id
            return

class ProcessListener(object):
    def __init__(self,
                 ip,            # ip address of the machine this listener is on
                 port,          # port to communicate on with worker
                 resultQueue,   # put results on this queue to add hits to DB
                 queueManager,  # SGE queue to place jobs on
                 scriptsDir,     # where to write SGE scripts to
                 workingDir,     # where tmp files will be stored
                 (id, gPath1, gPath2, batch, ani)):
        # set up the listener
        self.ip = ip
        self.port = port
        self.resultQueue = resultQueue
        self.queueManager = queueManager
        self.scriptsDir = scriptsDir
        self.workingDir = workingDir
        self.id = id
        self.gPath1 = gPath1
        self.gPath2 = gPath2

        self.worker = TestProcessWorker(self.port, self.id, ani)

    def start(self):
        """start the listener"""
        # set up the socket
        context = zmq.Context()
        socket = context.socket(zmq.REP)
        socket.bind("tcp://*:%s" % self.port)
        #print "L [%d] : ProcessListener on port: %s" % (self.id, self.port)

        # set the worker going
        ret_str, sge_script_fn = self.queueManager.makeSgeScript(self.scriptsDir,
                                                                 self.workingDir,
                                                                 self.id,
                                                                 self.gPath1,
                                                                 self.gPath2,
                                                                 "tcp://%s:%d" % (self.ip, self.port)
                                                                 )
        self.queueManager.lodgeJob(ret_str, sge_script_fn)

        multiprocessing.Process(target=self.worker.start).start()

        # wait for result from worker and decode (blocking)
        result = jp.decode(socket.recv().decode("zlib"))
        result[1] = float(result[1])/1000.

        # place the result on the queue
        self.resultQueue.put(result)

        # tell the worker to exit
        socket.send("DIE")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Server(object):
    def __init__(self,
                 dbFileName,        # path to db file to work with
                 port = None,       # port to listen on, None if not required
                 ):
        self.port = port
        self.dbFileName = dbFileName
        self.lock = multiprocessing.Lock()
        self.highestHitId = -1
        self.queueManager = None
        self.ip = self.getIpAddress()

    def getIpAddress(self):
       s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
       s.connect(("gmail.com",80))
       return s.getsockname()[0]

    def importNewPairs(self,
                      pairs,     # csv containing information about new pairs to process
                      paths     # absolute paths to new genome fasta files referenced in pairs
                      ):
        print "Importing new pairs from: %s" % pairs
        print "Import full contig file paths from: %s" % paths

        # get an interface to the file
        II = ImportInterface(self.dbFileName)

        # import paths first
        with open(paths, "r") as paths_fh:
            II.importGenomes(paths_fh)

        # import genomes next
        with open(pairs, "r") as pairs_fh:
            II.importPairs(pairs_fh)

    def makeHits(self,
                 queueURL,       # the queue we'll be sending work to
                 commsPort,      # port we'll be asked for progress etc on
                 portRange,      # the total number of concurrent pairs we can handle
                 scriptsDir,     # where to write SGE scripts to
                 workingDir,     # where tmp files will be stored
                 batches=[]      # the order to do pairs in
                 ):
        """process any specified outstanding pairs"""
        print "Processing outstanding pairs"

        # set tmp dirs
        self.scriptsDir = scriptsDir
        self.workingDir = workingDir

        # get hold of the server queue
        self.queueManager = SGE(queueURL)

        VI = ViewInterface(self.dbFileName)
        contig_headers = {}
        if batches != []:
            for batch in batches:
                pairs = VI.getOutstandingPairs(batch=batch)
                if len(pairs) > 0:
                    contig_headers = VI.getContigHeaders()
                    self.processOutstandingPairs(pairs, portRange, contig_headers, batch=batch)
        else:
            pairs = VI.getOutstandingPairs()
            if len(pairs) > 0:
                contig_headers = VI.getContigHeaders()
                self.processOutstandingPairs(pairs, portRange, contig_headers)

    def processOutstandingPairs(self, pairs, portRange, contigHeaders, batch=None):
        """Set up a worker to process pairs"""
        port_range = [int(i) for i in portRange.split(":")]
        port_range = range(port_range[0], port_range[1]+1)
        num_threads = len(port_range)
        if batch == None:
            print "Processing %d pairs on %d threads" % (len(pairs), num_threads)
        else:
            print "Processing %d pairs on %d threads (batch: %d)" % (len(pairs), num_threads, batch)

        # get the highest hit id in the db
        result_queue = multiprocessing.Queue()
        # the None is needed!
        updater_process = multiprocessing.Process(target=self._updateHits, args=(contigHeaders, result_queue)).start()
        # give the updater time to fire up
        time.sleep(1)

        # implement a dodgy thread pool
        all_procs = []
        for pair in range(len(pairs)):
            # get the port for this thread
            port = port_range[pair%num_threads]

            # create a Listener and start it on it's own thread
            L = ProcessListener(self.ip,
                                port,
                                result_queue,
                                self.queueManager,
                                self.scriptsDir,
                                self.workingDir,
                                pairs[pair])
            all_procs.append(multiprocessing.Process(target=L.start).start())

            # kill some time while we're waiting for some of the threads to clear
            while len(multiprocessing.active_children()) >= (num_threads + 1):
                time.sleep(1)

        # make sure all the processes are done before we go on
        while len(multiprocessing.active_children()) > 1:
            time.sleep(1)

        # wait for the workers to finish up
        for proc in all_procs:
            if proc is not None:
                proc.join()

        # kill the updater process
        result_queue.put(None)

        # wait for it to complete
        if updater_process is not None:
            updater_process.join()

    def _updateHits(self, contigHeaders, resultQueue):
        while(True):
            # block waiting for the next item in the queue to appear
            hits = resultQueue.get(block=True, timeout=None)
            if hits == None:
                break

            # get an import interface
            II = ImportInterface(self.dbFileName)
            II.importHits(contigHeaders, hits)

    def runAsDaemon(self):
        """Run in the background and serve requests for TrackM view etc"""
        pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################
