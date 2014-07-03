#!/usr/bin/env python
###############################################################################
#                                                                             #
#    exceptions.py                                                            #
#                                                                             #
#    All the exceptions we'd like to raise.                                   #
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TM_Exception(BaseException): pass

#------------------------------------------------------------------------------
# SHELL  -- deal with failure of calling external applications like Nucmer etc.
class TM_ShellException(TM_Exception): pass
class TM_ExternalProgramFailedException(TM_ShellException): pass

#------------------------------------------------------------------------------
# COMMS  -- deal with failure to communicate with the server
class TM_CommsException(TM_Exception): pass
class TM_ConnectionException(TM_CommsException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

