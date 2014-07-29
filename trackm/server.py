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
__version__ = "0.2.1"
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
import datetime
import random
import jsonpickle as jp
import zlib
import socket
import sys

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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HitTask(object):
    """Encapsulate information needed to carry out a comparison between two genomes"""
    def __init__(self, (id, gPath1, gPath2, gid1, gid2, batch, ani) ):
        self.id = id            # unique ID for the comparison (pid)
        self.gPath1 = gPath1    # path to the first genome
        self.gPath2 = gPath2    # path to the second genome
        self.gid1 = gid1        # uid for the first genome
        self.gid2 = gid2        # uid for the second genome
        self.ani = ani          # highest ani between genome1 and genome2

    def __str__(self):
        return "\t".join([str(i) for i in [self.id,
                                           self.gPath1,
                                           self.gPath2,
                                           self.gid1,
                                           self.gid2,
                                           self.ani]
                          ]
                         )

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ProcessListener(object):
    def __init__(self,
                 ip,             # ip address of the machine this listener is on
                 port,           # port to communicate on with worker
                 taskQueue,      # consume tasks from this queue
                 resultQueue,    # put results on this queue to add hits to DB
                 queueManager,   # SGE queue to place jobs on
                 sgeBaseDir,     # where to write SGE scripts to
                 workingDir,     # where tmp files will be stored
                 resQLowerLimit, # halt adding new jobs to the queue till the result queu reaches this size
                 resQUpperLimit  # when the result queue reaches this size
                 ):
        # set up the listener
        self.ip = ip
        self.port = port
        self.taskQueue = taskQueue
        self.resultQueue = resultQueue
        self.queueManager = queueManager
        self.sgeBaseDir = sgeBaseDir
        self.workingDir = workingDir
        self.resQUpperLimit = resQUpperLimit
        self.resQLowerLimit = resQLowerLimit

    def start(self):
        """start the listener"""
        # set up the socket
        while True:
            # get the next task
            #print ">> DBG about to get task"

            # don't let the result queue get too big!
            if int(self.resultQueue.qsize()) > self.resQUpperLimit:
                while int(self.resultQueue.qsize()) > self.resQLowerLimit:
                    time.sleep(10)

            task = self.taskQueue.get(block=True, timeout=None)
            if task == None:
                # poison pill
                #print ">> DBG Last task"
                break

            #print ">> DBG got task: %d" % task.id
            
            # configure zmq to listen for the result
            socket_OK = False
            context = zmq.Context()
            socket = context.socket(zmq.REP)
            while not socket_OK:
                try:
                    socket.bind("tcp://*:%s" % self.port)
                    socket_OK  = True
                except zmq.ZMQError:
                    #print sys.exc_info()[0], "Job: %d" % task.id
                    time.sleep(1)

            # set the worker going
            ret_str, sge_script_fn = self.queueManager.makeSgeScript(self.sgeBaseDir,
                                                                     self.workingDir,
                                                                     task.id,
                                                                     task.gPath1,
                                                                     task.gPath2,
                                                                     task.gid1,
                                                                     task.gid2,
                                                                     task.ani,
                                                                     "tcp://%s:%d" % (self.ip, self.port)
                                                                     )
            exit_status = self.queueManager.lodgeJob(ret_str, sge_script_fn)

            if exit_status == 0:

                # TO DO: send the deets of this job off to an external management thread which
                # monitors the queue to make sure the job isn't just dropped. If it is then is can send
                # a "DIE" signal to this listener

                # wait for result from worker and decode (blocking)
                #print ">> DBG waiting for worker: %d" % task.id
                result = jp.decode(socket.recv().decode("zlib"))

                # check to see we've not been told to die
                if result == "DIE":
                    # we abandon the worker and simply exit
                    return

                # check to see that there was no issue running the worker
                # basically, check to see that the last item in the result array is
                # actually a hit
                if len(result) > 3:
                    # there is something on the end of this array
                    if result[-2] == "ERROR":
                        # something went wrong. Print it out!
                        # TODO use logging module
                        print self.id, self.gPath1, self.gPath2
                        print result[-1]

                result[1] = float(result[1])/1000.

                # place the result on the result queue
                self.resultQueue.put(result)
                #print ">> DBG %d The length of result queue is %d" % (task.id,int(self.resultQueue.qsize())) 
                # tell the worker to exit
                socket.send("DIE")
            else:
                print " OH NO! EXIT STATUS == %d : %d" % (exit_status, task.id)
                time.sleep(8)
                #print ">> DBG No time for sleeping: %d" % task.id

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Server(object):
    def __init__(self,
                 dbFileName,        # path to db file to work with
                 port=None,         # port to listen on, None if not required
                 ):
        self.port = port
        self.dbFileName = dbFileName
        self.highestHitId = -1
        self.queueManager = None
        self.ip = self.getIpAddress()
        #self.DBL = DbLogger(self.dbFileName)

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
                 sgeBaseDir,     # where to write SGE scripts to
                 workingDir,     # where tmp files will be stored
                 batches=[],     # the order to do pairs in
                 hitCache=1000   # number of hits to cache before writing to database
                 ):
        """process any specified outstanding pairs"""
        print ">> DBG " + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        print "Processing outstanding pairs"

        # set tmp dirs
        self.sgeBaseDir = sgeBaseDir
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
                    self.processOutstandingPairs(pairs, portRange, contig_headers, hitCache, batch=batch)
        else:
            pairs = VI.getOutstandingPairs()
            if len(pairs) > 0:
                contig_headers = VI.getContigHeaders()
                self.processOutstandingPairs(pairs, portRange, contig_headers, hitCache)

    def processOutstandingPairs(self, pairs, portRange, contigHeaders, hitCache, batch=None):
        """Set up a worker to process pairs"""
        port_range = [int(i) for i in portRange.split(":")]
        port_range = range(port_range[0], port_range[1]+1)
        num_threads = int(len(port_range))
        if batch == None:
            print "Processing %d pairs on %d threads" % (len(pairs), num_threads)
        else:
            print "Processing %d pairs on %d threads (batch: %d)" % (len(pairs), num_threads, batch)

        # use these queues to handle threading
        task_queue = multiprocessing.Queue()
        result_queue = multiprocessing.Queue()

        # make the list of tasks
        for pair in range(len(pairs)):
            task_queue.put( HitTask(pairs[pair]) )

        # place poison pills and set up the listeners
        listeners = []
        for t in range(num_threads):
            task_queue.put(None)
            port = port_range[t%num_threads]        # 1 port per process
            L = ProcessListener(self.ip,
                                port,
                                task_queue,
                                result_queue,
                                self.queueManager,
                                self.sgeBaseDir,
                                self.workingDir,
                                hitCache*2,
                                hitCache*4)

            listeners.append(multiprocessing.Process(target=L.start))

        # start the thread that will handle the placing of results in the DB
        result_handling_process = multiprocessing.Process(target=self._updateHits, args=(contigHeaders, result_queue, hitCache)).start()
        # give the updater time to fire up
        time.sleep(1)
        
        for l in listeners:
            l.start()
            time.sleep(1)       # not allowed to attack the ssh daemon

        for l in listeners:
            if l is not None:
                l.join()

        # kill the updater process
        result_queue.put(None)

        # wait for it to complete
        if result_handling_process is not None:
            result_handling_process.join()

    def _updateHits(self, contigHeaders, resultQueue, hitCache):
        II = ImportInterface(self.dbFileName,verbosity=-1)
        tmp = []
        while(True):
            current = resultQueue.get(block=True, timeout=None)
            #print "current capacity of resultQueue %d" % int(resultQueue.qsize())
            if current == None: 
                #print ">> DBG Current is None"
                II.importHits(contigHeaders, tmp)
                break
            else:
                tmp.append(current)
            
            #print ">> DBG Length of tmp is %d" % len(tmp)
             
            if len(tmp) >= hitCache: # set size limit of queue
                print ">> DBG " + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + " Length of tmp is %d :: current capacity of resultQueue %d" % (len(tmp), int(resultQueue.qsize()))
                II.importHits(contigHeaders, tmp)
                tmp = []

    def runAsDaemon(self):
        """Run in the background and serve requests for TrackM view etc"""
        pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################
