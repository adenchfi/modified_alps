#!/usr/bin/python
#
#  copyright (C) 2008  troy d. straszheim  <troy@resophonic.com>
#  
#  Distributed under the Boost Software License, Version 1.0.
#  See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt
#

#
#  "Passthru" driver, only responsible for 'flipping' exit status of
#  tests that are expected to fail.  See driver.py for the version
#  that is run when BOOST_BUILD_SLAVE is on, which does xmlizaton
#  and the like
#
import sys, os, os.path
from subprocess import Popen, PIPE

def verbose(what):
        print what

# ignored
# log = os.path.join(sys.argv[1], "Log.xml")
op = sys.argv[2]
# target = sys.argv[3]
argv = sys.argv[4:]
expect_fail = op.endswith("fail")

#
#  execute subprocess
#
subproc = None
returncode = None
ex = None
stdout = None
stderr = None
try:
    print argv
    subproc = Popen(argv, stdout=PIPE, stderr=PIPE)
    (stdout, stderr) = subproc.communicate()
except EnvironmentError, e:
    ex = e

returncode = subproc.returncode

if stdout:
    print stdout
if stderr:
    print stderr

if not ex: 
    # possibly flip the return code
    if not expect_fail:
        if not returncode:
            verbose("ok.")
        else:
            verbose("error.")
        sys.exit(returncode)
    else:
        if returncode != 0:
            verbose("ok.")
            sys.exit(0)
        else: 
            verbose("*** UNEXPECTED SUCCESS ***")
            sys.exit(1)  # we need an exit status for 'unexpected success'
else:
    # if there is an os error 'above' the actual exit status of the subprocess,
    # use the errno
    print "Error in build system: " + str(ex.strerror)
    sys.exit(ex.errno)

    

