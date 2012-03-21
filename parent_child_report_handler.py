#!/usr/bin/python2.4

# A query that will perform a compound het analysis on a 
# family with an affected and one or two parents

import os.path
import csv
import db
from time import clock, time

import cgi
import sys
sys.path.append('/home/Gleeson/database/src')
from web_utils import printHeader, printToServer
import queries

printHeader()

sys.stdout = open("debug/parent_child_out.txt",'w')
sys.stderr = open("debug/parent_child_err.txt",'w')

fields = cgi.FieldStorage()

######################################################################
# QUERIES
######################################################################
#
# This SQL query find all hom mutations in the affected that are shared
# with two parents, but not that are shared between both parents
#
#TODO
#modified this query to be reportable, need to tweak rest of script
TWO_PARENT_SHARED_HOM_QUERY = '''
    SELECT %s
    FROM Variants AS v
    INNER JOIN Isoforms AS i ON i.var_id = v.id
    INNER JOIN Genes AS g ON g.id = i.gene_id,
    (
        SELECT var_id
        FROM Calls
        WHERE pat_id IN (%%s,%%s)
        AND GT = 1
        GROUP BY var_id
        HAVING COUNT( pat_id ) =2
    ) AS t1

    INNER JOIN (

        SELECT var_id
        FROM Calls AS c
        INNER JOIN Variants AS v ON c.var_id = v.id
        WHERE AF < .03
        AND (
            pat_id = %%s
            AND GT = 2
            )
    ) AS t2 ON t1.var_id = t2.var_id
    WHERE v.id = t1.var_id''' % (queries.select_string)

SINGLE_PARENT_SHARED_HOM_QUERY = '''
    SELECT t1.var_id, v.chrom, v.pos
    FROM Variants AS v, (
        SELECT var_id
        FROM Calls
        WHERE pat_id = $$P1$$
        AND GT = 1
        GROUP BY var_id
        HAVING COUNT( pat_id ) =2
    ) AS t1

    INNER JOIN (

    SELECT var_id
    FROM Calls AS c
    INNER JOIN Variants AS v ON c.var_id = v.id
    WHERE AF < .03
    AND (
        pat_id = $$A$$
        AND GT = 2
        )
    ) AS t2 ON t1.var_id = t2.var_id
    WHERE v.id = t1.var_id
     '''

#
# Given a family ID, find all members of the family, and determine
# the affected by sorting. In this situation, we expect each family
# we are interogating to have parents that have a higher generation
# code than our affected patient
#
PATIENT_QUERY = ''' Select id,name,family,disease,generation,affected
        FROM TempPatients
        where family =%s 
        and valid = 1
        '''


######################################################################
# getPatientData
######################################################################

#TODO
#we can't always assume the affected is always the last
#can't always assume there is only one affected
#can't always assume there is only one child
def getPatientData(conn, familyid) :
    myquery = PATIENT_QUERY % familyid
    generation_gap = '''
    select max(generation), min(generation)
    from TempPatients
    where family = %s''' % familyid

    (maxgen, mingen) = conn.query( generation_gap )[0]
    if int(maxgen) - int(mingen) > 1 :
        print "whowhowhwow multigeneration thing going on"

    #parents assumed to be unaffected
    parents = []
    affecteds = []
    unaffecteds = []
    #allpatients = []
    for row in conn.query(myquery) :
        eyed,name,family,disease,generation,affected = row
        affected = int(affected)
        if affected == 1 and generation == maxgen :
            affecteds.append( eyed )
        elif affected == 0 and generation == mingen :
            parents.append( eyed )

    #sanity check
    if len(affecteds) != 1 :
        raise Exception("There are %d != 1 affecteds for family %s." % \
                        (len(affecteds), familyid))
    elif len(parents) != 2 :
        raise Exception("There are %d != 2 parents for family %s" % \
                        (len(parents), familyid))
    return affecteds, parents
# END getPatientData

try :

    conn = db.Conn("localhost")

    familyid = fields.getvalue("family_id")
    filename = fields.getvalue("filename")
    (affecteds, parents) = getPatientData(conn, familyid)
    q = TWO_PARENT_SHARED_HOM_QUERY % (parents[0], parents[1], affecteds[0])
    queries.makeReport( q, "1=1", filename )
    printToServer( '<a href="../reports/%s.tsv">Report Ready</a>' % filename )

except KeyError, (e) :
    printToServer( "You most likely left %s blank." % str(e) )

except Exception, (e) :
    printToServer( str(e) )
