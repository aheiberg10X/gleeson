#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import *
from plate_macros import addPlateMacro, deletePlateMacro
from queries import deletePlate
from db import Conn

printHeader()

fields = cgi.FieldStorage()

sys.stdout = open("debug/delete_plate_out.txt",'w')
sys.stderr = open("debug/delete_plate_err.txt",'w')

try :
    conn = Conn("localhost",dry_run=False)
    plate_id = int(fields["plate_id"].value)
    deletePlateMacro( plate_id )
    printToServer( "Plate metadata deleted from plate_ids.txt and plates.py" )
    printToServer( listify( deletePlate( conn, plate_id ) ) )

except Exception, (e) :
    printToServer( str(e) )
