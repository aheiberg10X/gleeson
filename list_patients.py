#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')

import queries
import db
from web_utils import *


def main() :
    conn = db.Conn("localhost")
    headers = ['id','name','family','disease','generation','affected']
    query = ''' Select %s
        FROM TempPatients
        where valid = 1
        ''' % ', '.join(headers)


    to_print = []
    to_print.append('''
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Patient List</title>
<link rel="stylesheet" type="text/css" href="../style/general.css" media="all">
''')

    to_print.append( tablifiy( conn.query( query ), headers[1:] ) )

    return ''.join(to_print)

