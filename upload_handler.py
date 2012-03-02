#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import printHeader
from plate_macros import addPlateMacro

printHeader()

fields = cgi.FieldStorage()

sys.stdout = open("debug/upload_out.txt",'w')
sys.stderr = open("debug/upload_err.txt",'w')

from importer import main
main(fields["plate_name"].value)

