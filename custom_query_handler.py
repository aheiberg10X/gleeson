#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import *
from queries import makeReport

printHeader()

fields = cgi.FieldStorage()

try :
    names = ["variant_sql","call_where","filename"]
    values = [fields[name].value for name in names]

    sys.stdout = open("debug/custom_out.txt",'w')
    sys.stderr = open("debug/custom_err.txt",'w')

    try :
        makeReport( values[0], values[1], values[2] )
        printToServer( '<a href="../reports/%s.tsv">report generated</a>' % values[2] )
    except Exception, (e) :
        printToServer( str(e) )

except KeyError, (e) :
    printToServer( "You most likely left %s blank." % str(e) )

except Exception, (e) :
    printToServer( str(e) )
