#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')

import queries
import db
from web_utils import *


def main() :
    to_print = [] 
    fin = open( "plate_ids.txt" )
    plate_list = []
    for line in fin.readlines() :
        print line
        plate_list.append( line.strip().split(':') )

    fin.close()

    to_print.append('''
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>Plate List</title>
<link rel="stylesheet" type="text/css" href="../style/general.css" media="all">
''')

    to_print.append( tablifiy( plate_list, ['Plate Name','id'] ) )
    return ''.join( to_print )

if __name__ == '__main__' :
    main()

