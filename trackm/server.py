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
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################


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

###############################################################################
###############################################################################
###############################################################################
###############################################################################
