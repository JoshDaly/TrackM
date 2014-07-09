#!/usr/bin/env python
###############################################################################
#                                                                             #
#    sge.py                                                                   #
#                                                                             #
#    Place jobs on an SGE queue                                               #
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
__credits__ = ["Josh Daly", "Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Dev"

###############################################################################

import argparse
import sys
import os
import shutil
import errno
from subprocess import Popen, PIPE

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class SGE(object):
    """Class that handles putting SGE jobs on a queue"""
    def __init__(self): pass

    def runCommand(self, cmd):
        """Run a command and take care of stdout

        expects 'cmd' to be a string like "foo -b ar"

        returns (stdout, stderr)
        """
        p = Popen(cmd.split(' '), stdout=PIPE)
        return p.communicate()

    def makeSurePathExists(self, path):
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def removePath(self, path):
        try:
            shutil.rmtree(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def sgeIt(self,
              sgeWorkingDir,        # full path to the working directory
              sgeId,                # a unique id for this job
              commands,             # a list of commands that need to be carried out
              numProcessors=1,      # how many processors do we want?
              ):
        """ Main wrapper"""
        # make the sge script
        self.makeSurePathExists(sgeWorkingDir)
        sge_fn = os.path.join(sgeWorkingDir, "sge_%s" % sgeId)
        sge_out_fn = "%s.out" % sge_fn
        sge_err_fn = "%s.err" % sge_fn
        sge_fn += ".sh"
        with open(sge_fn, "w") as sge_fh:
            sge_fh.write("#!/bin/bash\n")
            sge_fh.write("#$ -r y\n")
            sge_fh.write("#$ -o %s\n" % sge_out_fn)
            sge_fh.write("#$ -e %s\n" % sge_err_fn)
            for command in commands:
                sge_fh.write("%s\n" % command)

        # lodge the sge job on the queue
        self.runCommand("qsub %s" % sge_fn)

        # remove all the evidence NOTE: perhaps you want to parse some stuff first?
        self.removePath(sgeWorkingDir)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
