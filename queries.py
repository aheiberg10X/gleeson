import db
import re
import csv
import sys
import os.path
import broad
import globes
#configname = sys.argv[1]
#(head,tail) = os.path.split(configname)
#(filename,ext) = os.path.splitext(tail)


COVERAGE = 8

#all reports generally need the same columns
#these are they, broken up by table
vcols = ["id","chrom","pos","dbSNP","ref","mut","type","qual",
         "filter","AF","granthamScore","scorePhastCons",
         "consScoreGERP","distanceToSplice","AfricanHapMapFreq",
         "EuropeanHapMapFreq", "AsianHapMapFreq","clinicalAssociation"]

icols = ["functionGVS","polyPhen","codon_pos","codon_total","ref_aa","mut_aa","accession"]

gcols = ["geneSymbol","omim_disease"]

#going in the output, to format a query row to match this ordering
#modify formatQueryRow
column_headers = ["chrom", "pos", "dbSNP", "ref", "mut", "filter", \
                  "gene", "accession", "AF", \
                  "functionGVS", "polyPhen", "AA_Change", "AA_Pos", \
                  "granthamScore", "scorePhastCons", "consScoreGERP", \
                  "distanceToSplice", "omim", "clinicalAssociation", \
                  "GT:DP:GQ", "#HomShares", "Hom Shares", \
                  "#HetShares", "Het Shares"]

vcols_string = ', '.join(["v.%s" % c for c in vcols])
icols_string = ', '.join(["i.%s" % c for c in icols])
gcols_string = ', '.join(["g.%s" % c for c in gcols])
select_string = "%s, %s, %s" % (vcols_string, icols_string, gcols_string)
dont_want = ["intron","near-gene-5","intergenic","near-gene-3","coding-synonymous","coding-notMod3"]
gvs = ["functionGVS <> '%s'" % dw for dw in dont_want]
gvs = ' and '.join(gvs)

# input: a query row, where SELECT was vcols+icols+gcols
# return: the values corresponding to those of column_headers
#    We queried for vcols, icols, gcols and want to print out the appropriate
#    values given the column_headers we've chosen.  
#    THis functions takes a row returned by the query and selects only the
#    stuff we care to print
#    Note this doesn't get us all the way.  This list will stil be missing
# GT:DP:GQ and all the share information

#offset lets us SELECT $offset arbitrary columns before the usual
#vcols+gcols+icols
def formatQueryRow( row, offset=0 ) :
    output = []
    #basic var stuff
    output.extend( row[offset+1:offset+6] )
    #filter
    output.append( row[offset+8] )
    #gene
    output.append( row[-2] )
    #accession
    output.append( row[-3] )
    #AF
    output.append( row[offset+9] )
    #functionGVS
    output.append( row[-9] )
    #polyPhen
    output.append( row[-8] )
    #ref/mut aa
    output.append( "%s/%s" % (row[-5],row[-4]) )
    #pos/tot
    output.append( "%s/%s" % (row[-7],row[-6]) )
    #grantham,phast,gerp,splice
    output.extend( row[offset+10:offset+14] )
    #omim 
    output.append( row[-1] )
    #clinicalAssociation
    output.append( row[offset+17] )
    return output

#take some information about a row (coming from either Calls or JointCalls)
#and turn it into a string
def callString( calls_row, table ) :
    if table == 'Calls' :
        GT = broad.encodeGT( int(calls_row[3]) )
        dp_gq = [str(t) for t in calls_row[4:6]]
        return ':'.join( [GT] + dp_gq )
    elif table == 'JointCalls' :
        return ", ".join( [str(t) for t in calls_row] )

#returns 3 lists, noinfs, hets, and homs
#each list is comprised of tuples (pat_id, pat_name, callInfo)
def getPatients( conn, var_id, where_clause="", exclude=[], table="Calls" ) :

    q = '''
    select c.*, p.name, p.id
    from %s as c inner join Patients as p on c.pat_id = p.id
    where p.valid = 1 and c.var_id = %d %s''' % (table, var_id, where_clause)

    if table == 'Calls' :
        noinfs, homs, hets = [],[],[]
        lookup = {0 : noinfs, \
                  1 : hets, \
                  2 : homs}

    #quick extension for a project with Jeong-ho, Sangwoo, and Vineet
    elif table == 'JointCalls' :
        lookup = {}
        #genotype codes corresponding to:
        #AA_AA,AA_AB,AA_BB,AB_AA,AB_AB,AB_BB,BB_AA,BB_AB,BB_BB
        for gt in [0,1,2,10,11,12,20,21,22] :
            lookup[gt] = []

    r = conn.iterateQuery( q )
    for row in r :
        call_row = row[:-2]
        call_string = callString(call_row, table)
        GT = row[3]
        if int(row[-1]) not in exclude :
            lookup[int(GT)].append( [call_row[1],row[-2],call_string] )

    return lookup

#find variants where child is hom and both parents hets
def parentFind( parent1, parent2, child) :
    conn = db.Conn("localhost")
    query = '''
    select t1.*
    from (
    select v.id, chrom, pos
    from Calls as c inner join Variants as v on v.id = c.var_id
    where (pat_id = %d or pat_id = %d) and GT = 1 and AF < .1
    group by var_id
    having count(pat_id) = 2) as t1
    inner join
    (select v.id, chrom, pos
     from Calls as c inner join Variants as v on v.id = c.var_id
     where pat_id = %d and GT = 2 and AF < .1) as t2 on t1.id = t2.id''' \
    % (parent1, parent2, child)

#an old report for Eric, find shared Calls between two pats
def intersect() :
    conn = db.Conn("localhost")
    conn2 = db.Conn("localhost")
    cols = "g.id, v.id, i.id, v.chrom, v.pos, v.dbSNP, v.ref, v.mut, v.filter, g.geneSymbol, v.AF, i.functionGVS, i.polyPhen, i.ref_aa, i.mut_aa, i.codon_pos, i.codon_total, v.granthamScore, v.scorePhastCons, v.consScoreGERP, v.distanceToSplice, g.omim_disease"
 
    query = '''
    SELECT %s FROM Isoforms AS i
    INNER JOIN (
        SELECT var_id
        FROM Calls
        WHERE GT =1 AND ( pat_id =899 OR pat_id =1077  )
        GROUP BY var_id
        HAVING count( pat_id ) =2
    ) AS t ON t.var_id = i.var_id
    INNER JOIN Variants AS v ON v.id = i.var_id
    INNER JOIN Genes AS g ON i.gene_id = g.id
    WHERE TYPE = 1  AND ( functionGVS LIKE 'missense%%'
                       OR functionGVS = 'nonsense'
                       OR functionGVS LIKE 'stop%%'
                       OR functionGVS LIKE 'frameshift%%' )

    ORDER BY i.gene_id, i.id''' % cols

    fh = open("intersect.txt",'wb')
    fout = csv.writer( fh, \
                    delimiter='\t', \
                    quoting=csv.QUOTE_MINIMAL )

    def writeOut( fout, buffer, var_count ) :
        if var_count > 1 :
            fout.writerows(buffer)

    #print headers
    fout.writerow( [t.split('.')[1] for t in cols.split(', ')][3:] + \
                   ["Hom Count", "Hom Shares", "Het Count", "Het Shares"] )

    prev_gid = -42
    prev_vid = -42
    var_count = 0
    buffer = []
    for i,row in enumerate(conn.query(query)) :
        if i % 1000 == 0 :
            print i
        (gid,vid,iid) = row[0],row[1],row[2]

        (noinfs, hets, homs) = getPatients( conn, vid, exclude=[899,1077] )
        het_string = '; '.join([ht[1] for ht in hets])
        hom_string = '; '.join([ht[1] for hm in homs])
        output = row[3:]+(len(homs),hom_string,len(hets),het_string)

        if gid != prev_gid :
            writeOut( fout, buffer, var_count )
            buffer = [output]
            var_count = 1
            prev_gid = gid
            prev_vid = -42
        else :
            buffer.append( output )
            if vid != prev_vid :
                var_count += 1
                prev_vid = vid
    writeOut( fout, buffer, var_count )
    fh.close()

#strip off the <table_alias> from <table_alias>.<column_name>
def makeColsReadable( cols ) :
    return [c.split('.')[1] for c in cols]


#var_query : a query on the variants,isoforms, and genes table
#            the columns selected out must be vcols+icols+gcols
#Each variant returned by the query is formatted and output,
#along with a rollup of the het and hom pats
#The returned patients can be filtered by call_where, which is supplied
#to getPatients().  eg where DP >= 8 and ...
#the report will be named outfile_name
#call_detail = True will report the call info alongside the patient name
def makeReport( var_query, call_where, outfile_name, call_detail=False ) :
    fout = open( '/home/Gleeson/database/src/html/reports/%s.tsv' % outfile_name, 'w')
    csvout = csv.writer( fout, delimiter='\t' )
    csvout.writerow( column_headers )

    conn = db.Conn("localhost")
    conn2 = db.Conn("localhost")

    #check the select statement of the incoming query for correct columns
    no_white = re.compile(r'\s+')
    input_query = no_white.sub( "", var_query ).lower()
    target = no_white.sub( "", "select %s" % select_string ).lower()
    if not input_query.startswith( target ) :
        raise Exception("The entered query does not start with 'SELECT %s'" % select_string)

    for varix,r in enumerate(conn.query( var_query )) :
        output_row = formatQueryRow( r )
        #Report is not per patient, call info column no meaning here
        output_row.append("-")
        lookup = getPatients(conn2, r[0], " and (%s)" % call_where)
        hets = lookup[1]
        homs = lookup[2]
        if call_detail :
            hom_pats = [p[1]+p[2] for p in homs]
            het_pats = [p[1]+p[2] for p in hets]
        else :
            hom_pats = [p[1] for p in homs]
            het_pats = [p[1] for p in hets]

        num_homs = len(hom_pats)
        hom_string = '; '.join(hom_pats)

        het_string = '; '.join(het_pats)
        num_hets = len(het_pats)

        output_row.extend( [num_homs,hom_string,num_hets,het_string] )
        csvout.writerow( output_row )

    fout.close()

#make a hom and het report for each patient
#run this after a new plate gets loaded (and you have updatedAF)
#files goes to output/
def familyReports() :
    #print "changed getPatients(), things will break"
    outdir = globes.OUT_DIR
    conn = db.Conn("localhost")

    conn2 = db.Conn("localhost")
    print "connetions made"

    num_vcols = len(vcols)

    #handles to the file names
    fouts = {}
    #buffers to accumulate writes to the fouts
    fbuffers = {}

    query = '''select name from Patients'''

    #open files for each patient, one for hets, one for homs
    #print the column headers for each file
    #also initialize the buffer space
    for r in conn.iterateQuery( query ) :
        patient = broad.sanitizePatientName( r[0] )
        fouts[patient] = {}
        fbuffers[patient] = {}
        for gt in ["hets","homs"] :
            filename = '%s/%s_%s.tsv' % (outdir,patient,gt)
            fouts[patient][gt] = filename

            f = open(filename, 'wb')
            fout = csv.writer( f,\
                               delimiter='\t', \
                               quoting=csv.QUOTE_MINIMAL )
            #fouts[patient][gt].writerow( column_headers )
            fout.writerow( column_headers )
            f.close()

            fbuffers[patient][gt] = []

    #what are the interesting variants
    query = '''select %s, %s, %s
           from Variants as v inner join Isoforms as i on v.id = i.var_id
                              inner join Genes as g on g.id = i.gene_id
           where (%s)  and v.AF < 0.1
           order by AF''' % (vcols_string, icols_string, gcols_string, gvs)

    print query

    #write the buffers out to the respective files
    def flush() :
        for pat in fbuffers :
            for gt in ["homs","hets"] :
                f = open( fouts[pat][gt], 'a')
                fout = csv.writer( f,\
                                   delimiter='\t', \
                                   quoting=csv.QUOTE_MINIMAL )
                fout.writerows( fbuffers[pat][gt] )
                f.close()
                #fouts[pat][gt].writerows( fbuffers[pat][gt] )
                fbuffers[pat][gt] = []


    for varix,r in enumerate(conn.query( query )) :
        #do a buffer flush
        if varix % 10000 == 0 :
            flush()
            print varix

        var_id = r[0]
        isIndel = int(r[6]) == 2
        if isIndel :
            where = ""
        else :
   #        only look at patients meeting these call reqs
            where = " and c.DP >= %d" % COVERAGE
        lookup = getPatients( conn2, var_id, where_clause=where )
        (noinfs,hets,homs) = [lookup[gt] for gt in [0,1,2]]

        if len(hets) == len(homs) == 0 : continue

        hom_pats = [p[1] for p in homs]
        num_homs = len(hom_pats)
        hom_string = '; '.join(hom_pats)

        het_pats = [p[1] for p in hets]
        het_string = '; '.join(het_pats)
        num_hets = len(het_pats)

        output_row = formatQueryRow( r )

        for ix,(pat_id,pat,call) in enumerate(homs) :
            pat = broad.sanitizePatientName( pat )
            hom_shares = hom_pats[:ix] + hom_pats[ix+1:]
            new_hom_string = '; '.join(hom_shares)
            fbuffers[pat]["homs"].append( output_row + \
                                         [call, num_homs-1, new_hom_string, \
                                          num_hets, het_string] )

        for ix,(pat_id,pat,call) in enumerate(hets) :
            pat = broad.sanitizePatientName( pat )
            het_shares = het_pats[:ix] + het_pats[ix+1:]
            new_het_string = '; '.join(het_shares)
            fbuffers[pat]["hets"].append( output_row + \
                                         [call, num_homs, hom_string, \
                                          num_hets-1, new_het_string] )
            #fouts[pat]["hets"].writerow( output_row + \
                                         #[call, num_homs, hom_string, \
                                          #num_hets-1, new_het_string] )

    flush()

#once a new VCF has been loaded, AF will have to be updated
def updateAF(conn) :
    query = "select count(*) from Patients where valid = 1"
    num_pats = conn.queryScalar( query, int )
    query = '''update Variants, (select var_id, sum(GT)/%d as newAF
                                 from Calls as c inner join Patients as p
                                      on c.pat_id = p.id
                                 where p.valid = 1 and DP >= %d
                                 group by var_id) as t
               set AF = t.newAF
               where id = t.var_id''' % (2*num_pats, COVERAGE)
    conn.put( query )

if __name__ == '__main__' :
    #print vcols_string,icols_string,gcols_string
    #intersect()

    #conn = db.Conn("localhost", dry_run=False)
    #(noinfs, hets, homs) = getPatients( conn, 123 )
    #print hets
    print "hwerew"
    familyReports()
    #updateAF(conn)
    
    #print genQ1(params)
    #printNewColsDict()
    
    #conn = db.Conn()
    #(noinfs,hets,homs) = getPatients( conn, 1 )
    #print homs

    #makeReport()


    #delete from Isoforms where var_id in (select id from Variants where type = 2);# 53994 row(s) affected.


    #delete from Calls where var_id in (select id from Variants where type = 2);# 991508 row(s) affected.


############################################################################
###        Very Old Stuff
############################################################################
#f = open( configname )
#params = json.loads( f.read() )
#f.close()
#params = {"interval" : "Genes", \
          #"include_atleast" : 2, \
          #"include_genotype" : "Hom", \
          #"call_depth" : 8, \
          #"call_qual" : False, \
          #"group" : "200,202,203", \
          #"exclude_genotype" : "Hom", \
          #"exclude_group" : "201,204", \
          #"variant_restrictions" : ''' v.dbSNP = '.'
                  #and (i.ss_functionGVS = 'missense' 
                    #or i.ss_functionGVS = 'nonsense' 
                    #or i.ss_functionGVS = 'splice-3' 
                    #or ss_functionGVS = 'splice-5'
                    #or (ss_functionGVS is null and effect = 'NON_SYNONYMOUS_CODING'))'''}


#delete from Variants where type = 2# 35200 row(s) affected.
q0 = "delete from TempIntVars"

def printNewColsDict() :
    cols = ["tiv1.int_id", "tiv1.var_id", "v.id", "v.chrom", \
            "v.pos", "v.dbSNP", "v.ref", "v.mut", "i.ref_aa", \
            "i.mut_aa", "v.qual", "v.filter", "v.AF", \
            "v.ss_granthamScore", "v.ss_scorePhastCons", \
            "v.ss_consScoreGERP", "v.ss_distanceToSplice", \
            "v.ss_cDNAPosition", "v.ss_AfricanHapMapFreq", \
            "v.ss_EuropeanHapMapFreq", "v.ss_AsianHapMapFreq", \
            "v.sift_region", "v.sift_type", "v.sift_prediction", \
            "v.sift_score", "v.sift_omim", "i.effect", "i.ss_functionGVS", \
            "i.ss_polyPhen", "i.transcript_id", "i.exon_rank", \
            "i.codon_pos", "i.codon_total"]
    for i,c in enumerate(cols) :
        print "%d : '%s', \\" % (i,c)


#does it make more sense to put variant filtering here?
def genQ1( params ) :
    v = params["interval"]
    template = "inner join %s as gi on gi.var_id = c.var_id"
    if v == "Genes" : interval_table = template % "GeneIntervals"
    elif v == "Exons" : interval_table = template % "ExonIntervals"
    else : interval_table = "join VariantIntervalPlaceholder as gi"

    v = params["include_atleast"]
    include_atleast = int(v)

    v = params["include_genotype"]
    if v == "Hom" : include_weigher = "IndifferentHomWeigher"
    elif v == "Het" : include_weigher = "IndifferentHetWeigher"
    else : include_weigher = "BothWeigher"

    v1,v2 = params["call_depth"], params["call_qual"]
    call_restrictions = []
    if v1 : call_restrictions.append( "c.DP >= %s" % v1 )
    if v2 : call_restrictions.append( "c.GQ >= %s" % v2 )
    call_restrictions = " and ".join( call_restrictions )

    ors = ["%%s = %d" % int(pat.strip()) for pat in params["group"].split(',')]
    calls_group = "(%s)" % ' or '.join([t % "c.pat_id" for t in ors])
    patient_group = "(%s)" % ' or '.join([t % "id" for t in ors])

    ors = ["%%s = %d" % int(pat.strip()) for pat in params["exclude_group"].split(',')]
    calls_exclude_group = "(%s)" % ' or '.join([t % "c.pat_id" for t in ors])
    patient_exclude_group = "(%s)" % ' or '.join([t % "id" for t in ors])

    v = params["exclude_genotype"]
    if v == "Hom" : penalty_weigher = "PenaltyHomWeigher"
    elif v == "Het" : penalty_weigher = "PenaltyHetWeigher"
    else : penalty_weigher = "PenaltyBothWeigher"

    variant_restrictions = params["variant_restrictions"]
    null_weight = 1000

    return \
    '''insert into TempIntVars
    select gi.int_id, gi.name, c.var_id, sum(weight)
    from Calls as c
          left join
         (
           (select id,GT,weight
            from Patients join %s
            where %s)
                   union
           (select id,GT,weight
            from Patients join %s 
            where %s)
         ) as w on
                   w.id = c.pat_id
               and w.GT = c.GT
        %s
        inner join (
                select v.id
                from Variants as v inner join Isoforms as i on i.var_id = v.id
                where %s
                group by v.id) as vi on vi.id = c.var_id
    where ((%s or %s))
       and %s
    group by c.var_id
    having sum( w.weight ) < 0;''' \
    % (include_weigher, patient_group, penalty_weigher, patient_exclude_group, interval_table, variant_restrictions, calls_group, calls_exclude_group, call_restrictions) #, null_weight)

dcols ={0 : 'tiv1.int_id', \
        1 : 'i.gene', \
        2 : 'v.id', \
        3 : 'v.chrom', \
        4 : 'v.pos', \
        5 : 'v.dbSNP', \
        6 : 'v.ref', \
        7 : 'v.mut', \
        8 : 'i.ref_aa', \
        9 : 'i.mut_aa', \
        10 : 'v.qual', \
        11 : 'v.filter', \
        12 : 'v.AF', \
        13 : 'v.ss_granthamScore', \
        14 : 'v.ss_scorePhastCons', \
        15 : 'v.ss_consScoreGERP', \
        16 : 'v.ss_distanceToSplice', \
        17 : 'v.ss_cDNAPosition', \
        18 : 'v.ss_AfricanHapMapFreq', \
        19 : 'v.ss_EuropeanHapMapFreq', \
        20 : 'v.ss_AsianHapMapFreq', \
        21 : 'v.sift_region', \
        22 : 'v.sift_type', \
        23 : 'v.sift_prediction', \
        24 : 'v.sift_score', \
        25 : 'v.sift_omim', \
        26 : 'i.effect', \
        27 : 'i.ss_functionGVS', \
        28 : 'i.ss_polyPhen', \
        29 : 'i.transcript_id', \
        30 : 'i.exon_rank', \
        31 : 'i.codon_pos', \
        32 : 'i.codon_total'}

def genQ2( params ) :
    if params['interval'] != 'Variants' :
        group_level = "int_id"
    else :
        group_level = "var_id"
    print "gl: %s" % group_level
    return \
    '''select %s from TempIntVars tiv1 inner join
         (select %s
          from TempIntVars tiv2
          group by %s
          having sum(weight) <= -%s) as t on tiv1.%s = t.%s
         inner join Variants as v on v.id = tiv1.var_id
         inner join Isoforms as i on i.var_id = v.id
         left join Genes as g on g.ucsc_id = tiv1.name
    where %s''' \
    % (', '.join( dcols.values() ), group_level, group_level, params["include_atleast"], group_level, group_level, params["variant_restrictions"]) 

int_cols = ["Interval Name","int_id"]
var_cols = ["Var ID", \
            "Chrom", \
            "Pos", \
            "dbSNP", \
            "Ref", \
            "Mut", \
            "Ref AA", \
            "Mut AA"]
iso_cols = []
#def makeReport(params) :
    #try :
        #conn = db.Conn()
        #r = conn.cur.execute( q0 )
        #q1 = genQ1(params)
        #q2 = genQ2(params)
        #print "<br /><br />Q1: %s" % q1
        #print "<br /><br />Q2: %s" % q2
       # 
        #r = conn.cur.execute( q1 )
        #r = conn.iterateQuery( q2)
        #reportname = '../../html/reports/%s.tsv' % params["filename"]
        #returnname = '../../reports/%s.tsv' % params["filename"]
        #freport = csv.writer( open(reportname,'wb'), \
                    #delimiter='\t', \
                    #quoting=csv.QUOTE_MINIMAL )
        #prev_int_id = -1
        #prev_var_id = -1
#
        ##this is duplicating call)_restrictions in genQ1... can do better
        #where_clause = ' and c.DP >= %s and c.GQ >= %s' % (params['call_depth'], params['call_qual'])
        #freport.writerow( makeColsReadable(dcols.values()) + ['Hom Shares','Het Shares'])
        #for row in r :
            #var_id = row[2]
            #(noinfs,hets,homs) = ['; '.join([str(t[1]) for t in p]) for p in getPatients( conn, var_id, where_clause )]
            #freport.writerow( list(row) + [homs, hets] )
#
#
            ####OLD PRETTY(IER) PRINTING CODE
            ##int_id, int_name = row[:2]
            ##var_id = row[2]
            ###print int_id, var_id
            ##if int_id != prev_int_id :
                ##freport.writerow( [] )
                ##freport.writerow( [int_id,int_name] )
                ##freport.writerow( ['',''] + makeColsReadable(dcols.values())[2:] + ['Hom Shares','Het Shares'])
                ##prev_int_id = int_id
##
            ###if var_id != prev_var_id :
            ##freport.writerow( ['',''] + list(row[2:]) + [homs, hets]  )
            ##prev_var_id = var_id
        #return returnname
    #except Exception, (e) :
        #print "<br /><br />%s" % str(e)
############
