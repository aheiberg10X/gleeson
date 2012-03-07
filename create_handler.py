#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import printHeader, printToServer
from plate_macros import addPlateMacro

printHeader()

fields = cgi.FieldStorage()

sys.stdout = open("debug/create_out.txt",'w')
sys.stderr = open("debug/create_err.txt",'w')

names = ["plate_name","snp_vcf_file","indel_vcf_file","snp_seattle_file","indel_seattle_file"]

values = [fields[name].value for name in names]
try :
    addPlateMacro( *values )
    printToServer( 'Plate created successfully. Go <a href="../upload_plate_form.html">here</a> to add it to the database' )
except Exception, (e) :
    printToServer( str(e) )
