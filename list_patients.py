#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')

import queries
import db
from web_utils import *


LOGIN_FAIL = "incorrect password"

fields = cgi.FieldStorage()
passwd = fields['passwd'].value

def main() :
    printToServer( "Content-Type: text/html" )
    printToServer( "" )
    printToServer("hello")
