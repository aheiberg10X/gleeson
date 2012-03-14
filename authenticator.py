#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')

from web_utils import *

printHeader()

sys.stdout = open("debug/authenticate_out.txt",'w')
sys.stderr = open("debug/authenticate_err.txt",'w')


LOGIN_FAIL = "incorrect password"

fields = cgi.FieldStorage()
try :
    passwd = fields["passwd"].value

    if passIsValid( passwd ) :
        content_file = fields["content_file"].value
        doExec = bool(int(fields["doExec"].value))
        if doExec :
            exec( "import %s" % content_file )
            printToServer( exec( "%s.main()" % content_file ) )
        else :
            fin = open( "/home/Gleeson/database/src/html/forms/%s" % \
                        content_file )
            content = fin.read()
            fin.close()
        printToServer( content )
    else :
        printToServer( LOGIN_FAIL )

except KeyError :
    printToServer( "You got here from somewhere strange.  Try accessing from gleesonentry.ucsd.edu" )


