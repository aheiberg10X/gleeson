import db
import csv
import simplejson as json
import sys
import os.path
import broad
import globes
#configname = sys.argv[1]
#(head,tail) = os.path.split(configname)
#(filename,ext) = os.path.splitext(tail)


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

q0 = "delete from TempIntVars"

def printNewColsDict() :
    cols = ["tiv1.int_id", "tiv1.var_id", "v.id", "v.chrom", \
            "v.pos", "v.dbSNP", "v.ref", "v.mut", "v.ref_aa", \
            "v.mut_aa", "v.qual", "v.filter", "v.AF", \
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
        1 : 'g.geneSymbol', \
        2 : 'v.id', \
        3 : 'v.chrom', \
        4 : 'v.pos', \
        5 : 'v.dbSNP', \
        6 : 'v.ref', \
        7 : 'v.mut', \
        8 : 'v.ref_aa', \
        9 : 'v.mut_aa', \
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

def callString( calls_row ) :
    GT = broad.encodeGT( int(calls_row[2]) )
    dp_gq = [str(t) for t in calls_row[3:5]]

    return ':'.join( [GT] + dp_gq )

def getPatients( conn, var_id, where_clause="" ) :
    q = '''
    select c.*, p.name
    from Calls as c inner join Patients as p on c.pat_id = p.id
    where c.var_id = %d %s''' % (var_id, where_clause)
    noinfs, homs, hets = [],[],[]
    lookup = {0 : noinfs, \
              1 : hets, \
              2 : homs}
    r = conn.iterateQuery( q )
    for row in r :
        call_row = row[:-1]
        call_string = callString(call_row)
        GT = row[2]
        lookup[int(GT)].append( (call_row[1],row[-1],call_string) )

    return (noinfs,hets,homs)

#TODO
#need some kind of general rollup function here


def makeColsReadable( cols ) :
    return [c.split('.')[1] for c in cols]

def makeReport(params) :
    try :
        conn = db.Conn()
        r = conn.cur.execute( q0 )
        q1 = genQ1(params)
        q2 = genQ2(params)
        print "<br /><br />Q1: %s" % q1
        print "<br /><br />Q2: %s" % q2
        
        r = conn.cur.execute( q1 )
        r = conn.iterateQuery( q2)
        reportname = '../../html/reports/%s.tsv' % params["filename"]
        returnname = '../../reports/%s.tsv' % params["filename"]
        freport = csv.writer( open(reportname,'wb'), \
                    delimiter='\t', \
                    quoting=csv.QUOTE_MINIMAL )
        prev_int_id = -1
        prev_var_id = -1

        #this is duplicating call)_restrictions in genQ1... can do better
        where_clause = ' and c.DP >= %s and c.GQ >= %s' % (params['call_depth'], params['call_qual'])
        freport.writerow( makeColsReadable(dcols.values()) + ['Hom Shares','Het Shares'])
        for row in r :
            var_id = row[2]
            (noinfs,hets,homs) = ['; '.join([str(t[1]) for t in p]) for p in getPatients( conn, var_id, where_clause )]
            freport.writerow( list(row) + [homs, hets] )


            ###OLD PRETTY(IER) PRINTING CODE
            #int_id, int_name = row[:2]
            #var_id = row[2]
            ##print int_id, var_id
            #if int_id != prev_int_id :
                #freport.writerow( [] )
                #freport.writerow( [int_id,int_name] )
                #freport.writerow( ['',''] + makeColsReadable(dcols.values())[2:] + ['Hom Shares','Het Shares'])
                #prev_int_id = int_id
#
            ##if var_id != prev_var_id :
            #freport.writerow( ['',''] + list(row[2:]) + [homs, hets]  )
            #prev_var_id = var_id
        return returnname
    except Exception, (e) :
        print "<br /><br />%s" % str(e)

def indelReport() :
    outdir = globes.OUT_DIR
    conn = db.Conn("localhost")
    conn2 = db.Conn("localhost")
    print "connetions made"
    vcols = conn2.getColumns('Variants')
    icols = conn2.getColumns('Isoforms')

    #the general report
    fout = open("%s/indelReport.tsv" % (outdir),'w')
    freport = csv.writer( fout, \
                          delimiter='\t', \
                          quoting=csv.QUOTE_MINIMAL )
    freport.writerow( vcols + icols + ["Homs","Hets"] )

    #the per family reports
    fouts = {}
    query = '''select name from Patients'''
    for r in conn.iterateQuery( query ) :
        patient = broad.sanitizeFamilyName( r[0] )
        filename = '%s/%s_indels.tsv' % (outdir,patient)
        fouts[patient]=csv.writer( open(filename, 'wb'),\
                                   delimiter='\t', \
                                   quoting=csv.QUOTE_MINIMAL )
        fouts[patient].writerow( vcols + icols + ["GT:DP:GQ", "Hom Shares"] ) 

    query = '''select v.*, i.*
               from Variants as v inner join Isoforms as i on v.id = i.var_id
               where v.id < 1000'''

    for r in conn.iterateQuery( query ) :
        var_id = r[0]
        (noinfs,hets,homs) = getPatients( conn, var_id )
        for ix,(pat_id,pat,call) in enumerate(homs) :
            pat = broad.sanitizeFamilyName( pat )
            other_shares = '; '.join([p[1] for p in homs if p[1] != pat])
            fouts[pat].writerow( list(r) + [call,other_shares] )

        hom_names = [t[1] for t in homs]
        het_names = [t[1] for t in hets]
        (hom_string, het_string) = ['; '.join(t) for t in (het_names,hom_names)]
        freport.writerow( list(r)+[hom_string,het_string] )

    fout.close()

if __name__ == '__main__' :
    indelReport()
    #print genQ1(params)
    #printNewColsDict()
    
    #conn = db.Conn()
    #(noinfs,hets,homs) = getPatients( conn, 1 )
    #print homs

    #makeReport()
