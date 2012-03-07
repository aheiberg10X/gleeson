#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import printHeader, printToServer

printHeader()

fields = cgi.FieldStorage()


names = ["family_id"]
values = [fields[name].value for name in names]

sys.stdout = open("debug/parent_child_out.txt",'w')
sys.stderr = open("debug/parent_child_err.txt",'w')

try :
    #TODO
    #adapt erics parents.py script to work
    printToServer( "pass" ) 
except Exception, (e) :
    printToServer( str(e) )
