#!/usr/bin/python2.4

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')

import queries
import db
from web_utils import *

printHeader()

LOGIN_FAIL = "incorrect password"

fields = cgi.FieldStorage()
passwd = fields['passwd'].value

if 'passwd' in fields :
    if passIsValid( passwd ) :
        links = []
        links.append( authenticLink ( passwd, \
                                     "create_plate_form.html", \
                                     "Create a New Plate" ) )
        links.append( authenticLink( passwd, \
                                     "upload_plate_form.html", \
                                     "Upload a New Plate" ) )
        links.append( authenticLink( passwd, \
                                     "custom_query_form.html", \
                                     "Custom Report" ) )
        links.append( authenticLink( passwd, \
                                     "parent_child_form.html", \
                                     "Parent/Child Report" ) )

        print listify( links )

    else :
        print LOGIN_FAIL

else :
    print "You got here somewhere strange.  Try accessing from gleesonentry.ucsd.edu"

#if 'interval' not in fields and 'passwd' in fields :
    #if passIsValid( fields['passwd'].value ) :
        #print '<form name="query" action="hello.py">'
        #print '    <p>Report variants in  <select name="interval">'
        #print '        <option value="Genes">Genes</option>'
        ##print '        <option value="Exons">Exons</option>'
        #print '        <option value="Variants">Variants</option>'
        #print '    </select>'
        #print '    With >= <input type="text" name="include_atleast" size="3"/>(int)'
        #print '    <select name="include_genotype" />'
        #print '        <option value="Hom">Hom</option>'
        #print '        <option value="Het">Het</option>'
        ##print '        <option value="Both">Both</option>'
        #print '    </select> calls that all belong to someone in'
        #print '    <input type="text" name="group" value="ex. 200,202,203"/>'
        #print '    excluding those Genes/Exons/Varaints with <select name="exclude_genotype" />'
        #print '        <option value="Hom">Hom</option>'
        #print '        <option value="Het">Het</option>'
        #print '        <option value="Both">Both</option>'
        #print '    </select>'
        #print '    call that belong to someone in<input type="text" name="exclude_group" value="ex. 201,204" /></p>'
        #print '    <p>Calls must have depth >= <input type="text" name="call_depth" value="0"size="3"/> and qual >= <input type="text" name="call_qual" value="0" size="3"/></p>'
        #print "    <p>Variants considered must meet these requirements: <textarea rows=\"10\" cols=\"60\" name=\"variant_restrictions\"> v.dbSNP = '.' and (i.ss_functionGVS = 'missense' or i.ss_functionGVS = 'nonsense' or i.ss_functionGVS = 'splice-3' or ss_functionGVS = 'splice-5' or (ss_functionGVS is null and effect = 'NON_SYNONYMOUS_CODING')) </textarea> (edit this at your own risk)</p>"
        #print '    <p>Report Name:<input type="text" name="filename" value="myreport" /></p>'
        #print '    <p><input type="submit" value="Run Query" /></p>'
        #print '    <input type="hidden" name="passwd" value="%s" />' % fields['passwd'].value
        #print '    <p> Problem? andrew.heiberg@gmail.com </p>'
        #print '</form>'
#
        #print '<br /><br /> Patient IDs<br />'
        #conn = db.Conn()
        #r = conn.iterateQuery("select * from Patients")
        #for t in r :
            #print "<br />%d ----- %s" % (int(t[0]),t[1])
#
#
    #else :
        #print '<p>get lost</p>'
#
#elif 'passwd' in fields and 'interval' in fields :
    #if passIsValid( fields['passwd'].value ) : 
        #inputs = ["interval","include_atleast","include_genotype","group","exclude_genotype","exclude_group","call_depth","call_qual","variant_restrictions","filename"]
        #params = {}
        #not_valid = False
        #for i in inputs :
            #if i not in fields : 
                #print "You left [%s] blank" % i
                #not_valid = True
                #break
            #else :
                #params[i] = fields[i].value
#
        #if not not_valid :
            #print "making report"
            #url = query.makeReport( params )
            #print "<br /><br /><a href=\"%s\"> Download report </a>" % url
            ##print "<br /><a href=\"http://132.239.160.134/cgi-bin/andrew/hello.py\">Make another query</a>"
#
    #else :
        #print '<p>get lost</p>'
#
#else :
    #print "nothing to see here"
#
#
#

