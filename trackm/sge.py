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
    def __init__(self,
                 queueURL,          # where we will be submitting jobs
                 ):
        self.queueURL = queueURL

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

    def makeSgeScript(self,
                      sgeScriptsDir,        # full path to the location of sge scripts
                      baseWorkingDir,       # full path to the working directory
                      jobId,                # the unique id for this job (pid)
                      gPath1,               # full path to the first genome
                      gPath2,               # full path to the second genome
                      gid1,                 # the unique ID for the first genome
                      gid2,                 # the unique ID for the second genome
                      ani,                  # highest ani between the two servers
                      url,                  # url of the server to pass results to
                      ):
        """Create an SGE script that will run a worker task"""
        # make the sge script
        working_dir = os.path.join(baseWorkingDir, "job_%d" % jobId)

        sge_base_fn = os.path.join(working_dir, "sge_%s" % jobId)
        sge_out_fn = "%s.out" % sge_base_fn
        sge_err_fn = "%s.err" % sge_base_fn

        sge_script_fn = os.path.join(sgeScriptsDir, "job_%d.sh" % jobId)

        ret_str = "#!/bin/bash\n"
        ret_str += "#$ -r y\n"
        ret_str += "#$ -o %s\n" % sge_out_fn
        ret_str += "#$ -e %s\n" % sge_err_fn
        ret_str += "mkdir -p %s\n" % working_dir
        ret_str += "cd %s\n" % working_dir
        ret_str += "module load trackm\n"
        ret_str += "module load mummer\n"
        ret_str += "trackm worker %d %s %s %s %s %0.3f %s\n" % (jobId, gPath1, gid1, gPath2, gid2, ani, url)

        # clean up the filesystem
        ret_str += "cd ..\n"
        ret_str += "rm -rf %s\n" % working_dir
        ret_str += "rm %s\n" % sge_script_fn
        return (ret_str, sge_script_fn)

    def lodgeJob(self,
                 sgeStr,
                 scriptPath):
        """lodge an SGE job on the queue"""
        with open(scriptPath, 'w') as fh:
            fh.write(sgeStr)
        return
        self.runCommand("ssh -f %s 'qsub %s'" % (self.queueURL, scriptPath))

###############################################################################
###############################################################################
###############################################################################
###############################################################################
