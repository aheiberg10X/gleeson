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

printHeader()

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
        WHERE pat_id IN ($$P1$$,$$P2$$)
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
        FROM Patients
        where family = $$famid$$"
        and valid = 1
        '''


######################################################################
# getHomVarInfo
######################################################################
def getHomVarInfo(conn, affected, parents) :

    myquery = ""
    if len(parents) == 2 : 
        print "Performing Two Parent Query"

        myquery = TWO_PARENT_SHARED_HOM_QUERY
        myquery = myquery.replace("$$A$$", str(affected))
        myquery = myquery.replace("$$P1$$", str(parents[0]))
        myquery = myquery.replace("$$P2$$", str(parents[1]))
        #print myquery
    elif len(parents) == 1 :
        print "Performing Single Parent Query"

        myquery = SINGLE_PARENT_SHARED_HOM_QUERY
        myquery = myquery.replace("$$A$$", str(affected))
        myquery = myquery.replace("$$P1$$", str(parents[0]))
    #elif len(parents) == 0 :
    #    print "Performing Affected Het Query"

    #    myquery = AFFECTED_HET_QUERY
    #    myquery = myquery.replace("$$A$$", str(affected))
    else :
        sys.err("What kind of family is this?" )
        sys.exit(1)

    # Time the Query
    start = time()

    variantdata = []
    index =0
    for row in conn.query( myquery ):
        #print index,":",row
        variantdata.append( "%s:%s" % (row[1], row[2]) )
        index += 1

    # Print Run Time
    elapsed = (time() - start)
    print "Elapsed:",elapsed
    return variantdata
# END getHetVarInfo

######################################################################
# getPatientData
######################################################################

#TODO
#we can't always assume the affected is always the last
#can't always assume there is only one affected
#can't always assume there is only one child
def getPatientData(conn, familyid) :
    myquery = PATIENT_QUERY
    myquery = myquery.replace("$$famid$$", familyid)
    #parents assumed to be unaffected
    parents = []
    affecteds = []
    unaffecteds = []
    #allpatients = []
    for row in conn.query(myquery) :
        allpatients.append(row)
    allpatients = sorted(allpatients)
    affected = allpatients[-1]
    parents = allpatients[:-1]
    return affected, parents
# END getPatientData

######################################################################
# getPatientFile
######################################################################
#TODO
#want to generate this on the fly, not go digging around the reports
#NOw that we have 
def getPatientFile( patient, type ):
    patientfile = "/home/Gleeson/database/output/%s_%s.tsv" % (patient, type)
    assert( os.path.exists( patientfile ) )
    return patientfile
# END getPatientFile

######################################################################
# twoParentAnalysis
######################################################################
def twoParentAnalysis( variantdata, patientfile, outfile ) :

    # Open target file
    FILE = open( patientfile, "r" )
    reader = csv.reader( FILE, delimiter="\t" )

    header = []
    positions = {}
    uniquify = []
    for row in reader :
        # Get the Header
        if len(header) == 0 :
            header = row
            continue
        if len(row) < 3 :
            continue
        mykey = "%s:%s" % (row[0], row[1])

        #If the key has been found before, skip it
        if mykey in uniquify :
            continue

        # Find matches from p1
        if mykey in variantdata :
            print mykey
            if not positions.has_key(row[6]) :
                positions[row[6]] = []
            positions[row[6]].append(row)
            uniquify.append(mykey)
    FILE.close()

    # Print out compound hets where both p1 and p2 have a het
    FOUT = open(outfile, "w" )

    FOUT.write( "\t".join(header) +"\n" )
    for gene in sorted(positions) :
        for row in positions[gene] :
            FOUT.write( "\t".join(row)+"\n" )
    return
# END twoParentAnalysis

######################################################################
# Main
######################################################################
def main(familyid) :
    #file = open("./output.csv",'w')

    conn = db.Conn("localhost")

    # Get Affected and Parents
    (affected, parents) = getPatientData(conn, familyid)
    #print affected, parents
    print "Patients:",affected[1], parents[0][1], parents[1][1]

    # Find the PatientFile
    patientfile = getPatientFile( affected[0], "homs" )

    # Make the OutputFile
    outfile = "./reports/%s_homanalysis.csv" % familyid

    # Get Het Information
    if len(parents) == 2 :
        variantdata = getHomVarInfo(conn, affected[1], [parents[0][1], parents[1][1]])
        print variantdata[1:10]

        print "Matches VariantData:",len(variantdata)
        twoParentAnalysis( variantdata, patientfile, outfile )
    elif len(parents) == 1 :
        variantdata = getHomVarInfo(conn, affected[1], [parents[0][1]] )

        print "Matches Variantdata:",len(variantdata)
        #twoParentAnalysis( hetdata_p1, hetdata_a, patientfile, outfile )
    else :
        sys.err("This family doesnt work for this query")

# END MAIN

# Run Stuff
if __name__ == "__main__" :
    #familyid = "109"
    familyid = "1268"
    if len(sys.argv) > 1 :
        familyid = sys.argv[1]
    print "Family ID:",familyid
    main(familyid)
