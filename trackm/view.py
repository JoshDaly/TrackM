#!/usr/bin/env python
###############################################################################
#                                                                             #
#    view.py                                                                  #
#                                                                             #
#    The View class handles all TrackM based visualisations.                  #
#    NOTE: This class expects a TrackM server to be running in 'daemon' mode. #
#    when it is started.                                                      #
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

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    def __init__(self,
                 serverURL,         # URL of the commanding TrackM server
                 serverPort         # port of the commanding TrackM server
                 ):
        self.serverURL = serverURL
        self.serverPort = serverPort

    def connect(self):
        """Try connect to the TrackM server"""
        pass

    def plotSomething(self):
        """Place holder"""
        # an example entry point from the main TrackM entry point
        # Ideally, you should have a separate plot function
        # defined for each type of plot.

        # >>>>>>>>>> REMOVE THIS WHEN YOU HAVE CODE HERE <<<<<<<<<<<<<<<<<<
        from inspect import currentframe, getframeinfo
        frameinfo = getframeinfo(currentframe())
        print "Make me plot something! I live at File: %s Line: %s" % (frameinfo.filename, frameinfo.lineno)
        # >>>>>>>>>> END <<<<<<<<<<<<<<<<<<


###############################################################################
###############################################################################
###############################################################################
###############################################################################
