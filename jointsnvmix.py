import db
import queries
from collimator import CollimatorSource
from globes import toBaseN, splitIterator, chromNum
#TODO
#basically everything in variant has a more natural place to be
from variant import JointCall, Variant
import csv
from plates import Plate_JSM_HCD_1577_2_1, Plate_JSM_HME_1563_2_1, Plate_JSM_HME_1565_2_4, Plate_JSM_HME_1573_2_1, Plate_JSM_HME_1574_2_2, Plate_JSM_HME_1620_2_2


#"chrom", "position", "ref_base", "var_base"
headers = ["chrom", "pos", "ref", "mut", "normal_counts_a", "normal_counts_b", "tumour_counts_a", "tumour_counts_b", "p_AA_AA", "p_AA_AB", "p_AA_BB", "p_AB_AA", "p_AB_AB", "p_AB_BB", "p_BB_AA", "p_BB_AB", "p_BB_BB"]

indexOf = {}
for i,h in enumerate(headers) :
    indexOf[h] = i

#take the list of 9 probabilities and covert it into the correct int code
#AA_AA -> 0
#AA_AB -> 1
#AA_BB -> 2
#AB_AA -> 10
# etc
def convertGT( GT_probs ) :
    return toBaseN( GT_probs.index(max(GT_probs)), 3, [] )

def sortHelper( it ) :
    return (chromNum(it[0]),int(it[1]))

def sortFile( filename ) :
    fin = open(filename)
    fout = open("%s.sorted" % filename,'w')

    sorted_splts = sorted( splitIterator( fin ), key = sortHelper )
    fout.write( "\n".join( ['\t'.join(splt) for splt in sorted_splts] ) )

    fin.close()
    fout.close()

def skipper( splt ) :
    return "GL" in splt[0]

class JointSNVMixSource( CollimatorSource ) :
    def __init__(self, filename, pat_name) :
        self.pat_name = pat_name
        self.fin = open(filename)
        self.allow_absent = False
        self.group_repeats = False
        self.iterator = splitIterator( self.fin, \
                                       skipLine = skipper )

    def eqkey(self, it) :
        return (it[0],it[1],it[2],it[3])

    def integrator(self, target, splts) :
        #this is a one to many match with SeattleSeq
        assert len(splts) == 1

        splt = splts[0]

        #transition from variant fields to JointCall fields
        ix = 4
        splt[0] = chromNum( splt[0] )
        fields = dict( zip( headers[:ix], splt[:ix] ) )

        jc = JointCall( self.pat_name )
        jc.fields = dict( zip( headers[ix:], splt[ix:] ) )
        jc.fields["GT"] = convertGT( splt[-9:] )

        return Variant( fields, [jc] )


        #TODO
        #create a variant with the chrom,pos,ref,mut data
        #give variants a new type (3) to distinguish them as JointSNV vs GATK
        #create a table called jointCall
        #create a call with the AD and use PL fields for 9 probs
        #make calls with pat_name

def makeReport() :
#SELECT jc.var_id, COUNT( DISTINCT jc.pat_id ) 
#FROM  `JointCalls` AS jc
#INNER JOIN Variants AS v ON v.id = jc.var_id
#INNER JOIN Isoforms AS i ON i.var_id = v.id
#WHERE AF =0
#AND (
#functionGVS <>  'intron'
#AND functionGVS <>  'near-gene-5'
#AND functionGVS <>  'intergenic'
#AND functionGVS <>  'near-gene-3'
#AND functionGVS <>  'coding-synonymous'
#AND functionGVS <>  'coding-notMod3'
#)
#AND (
#p_AA_AB > .9
#OR p_AA_AB > .9
#OR p_AB_BB > .9
#)
#GROUP BY jc.var_id
#HAVING COUNT( DISTINCT jc.pat_id ) >1
#ORDER BY COUNT( DISTINCT jc.pat_id ) DESC /*inner join Isoforms as i on i.var_id = jc.var_id WHERE p_AA_BB > .90 or p_AA_AB > .9*/
#LIMIT 0 , 30

    fout = open("../output/joint.tsv",'w')
    csvout = csv.writer( fout, \
                    delimiter='\t', \
                    quoting=csv.QUOTE_MINIMAL )

    extend = ["patient"] + headers[4:]
    reuse = queries.column_headers[:-5] 
    csvout.writerow( reuse + extend )

    conn = db.Conn("localhost")
    conn2 = db.Conn("localhost")

    call_where = "and (p_AA_AB > .9 or p_AA_AB > .9 or p_AB_BB > .9)"
    call_where = ''' and (normal_counts_a / (normal_counts_a + normal_counts_b)
                          > .95 ) and normal_counts_a + normal_counts_b > 10'''
    query = '''
SELECT i.id, %s,%s,%s
FROM `JointCalls` as jc inner join Variants as v on v.id = jc.var_id
                        inner join Isoforms as i on i.var_id = v.id
                        inner join Genes as g on i.gene_id = g.id
where AF <= .0005 and
  (%s)
  %s
  and dbSNP is null
  and scorePhastCons >= .5
ORDER BY AF, jc.var_id, i.id''' % \
    (queries.vcols_string, queries.icols_string, queries.gcols_string, \
     queries.gvs, call_where)

    rows = conn.query( query )
    last_iid = -1
    for i, row in enumerate(rows) :
        print row
        outrow = queries.formatQueryRow(row,offset=1)

        iso_id = row[0]
        var_id = row[1]
        gt_lists = queries.getPatients( conn2, var_id, table='JointCalls', where_clause = call_where )
        #print " ----------- "
        #print gt_lists
        if last_iid != iso_id :
            for gt,calls in gt_lists :
                for call in calls :
                    pat_name = call[1]
                    rest = call[2].split(',')[6:]
                    together = outrow+[pat_name]+rest
                    together_string = "".join([str(t) for t in together])
                    csvout.writerow( together )
                    last = together_string

            last_iid = iso_id
        else :
            print "dupe!"
        #if i > 20 : break
        #hets_string = ", ".join([t[1] for t in hets])
        #homs_string = ", ".join([t[1] for t in homs])
        #if patient in hets_string or patient in homs_string :
            #already_called += 1
        #outrow.extend( ["-",len(homs),homs_string,len(hets),hets_string] )

    fout.close()

if __name__ == '__main__' :
    makeReport()

    #plates = [Plate_JSM_HME_1563_2_1(), Plate_JSM_HME_1565_2_4(), Plate_JSM_HME_1573_2_1(), Plate_JSM_HME_1574_2_2(), Plate_JSM_HME_1620_2_2()]
    #for plate in plates :
       #sortFile( plate.varFile('snp') )
