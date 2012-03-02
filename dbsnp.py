from collimator import Collimator, Source
import globes
import db
import csv

class DBSNPSource(Source) :
    def __init__(self, freq_file, allele_file ) :
        dbsnp_folder = "%s/dbsnp_132" % globes.ROOT_DIR
        self.freq_file = open( "%s/%s" % (dbsnp_folder,freq_file) )
        self.allele_file = open( "%s/%s" % (dbsnp_folder,allele_file) )

        #allele lookup
        self.dallele = {}
        for line in self.allele_file.readlines() :
            (eyed, allele) = line.split('\t')[:2]
            self.dallele[int(eyed)] = allele

        self.iterator = self.iterate()
        self.allow_absent = False
        self.group_repeats = False

    def iterate(self) :
        prev_rsid = -1
        #[rsid, allele, count, allele2, count]
        rsid_info = []
        for splt in globes.splitIterator( self.freq_file ) :
            #print splt
            try :
                rsid,alleleid,count = int(splt[0]), int(splt[1]), \
                                      int(float(splt[2]))
            except ValueError :
                print splt
                assert "Value" == "Error"

            if alleleid > 13 or alleleid == 5 or alleleid == 8 :
                # where the indels start
                # + (used for generic indels?) 
                # - (used for generic indels?)
                continue

            if rsid != prev_rsid :
                if prev_rsid > 0 :
                    yield rsid_info
                    prev_rsid = rsid_info[0]
                    rsid_info = []
                rsid_info.append(rsid)
                rsid_info.append( self.dallele[alleleid] )
                rsid_info.append( count )
                #rsid_info.append( {self.dallele[alleleid] : (count,freq)} )
                prev_rsid = rsid
            else :
                rsid_info.append( self.dallele[alleleid] )
                rsid_info.append( count )
                #rsid_info[1][self.dallele[alleleid]] = (count,freq)

    def eqkey(self, rsid_info) :
        return rsid_info[0]

    #target is (var_id, rsid, ref_freq, count)
    #first two will be provided by the variant table iterator
    def integrator(self, target, rsid_infos) :
        if len(rsid_infos) != 1 : 
            print "rsid_infos length != 1"
            assert False
        rsid_info = rsid_infos[0]
        count = 0
        for k in rsid_info[1].keys() :
            count += rsid_info[1][k][0]
        target["count"] = count
        target["alleles"] = rsid_info[1]
        return target

class VariantSource :
    def __init__(self,conn) :
        self.conn = conn
        self.allow_absent = True 
        self.group_repeats = False
        self.iterator = self.iterate()

    def iterate(self) :
        query = '''select id,dbSNP,ref,mut
                   from Variants
                   where dbSNP is not NULL
                   order by dbSNP'''
        for row in self.conn.query(query) :
            yield row

    def eqkey(self,it) :
        return int(it[1])  #rsid = it[1]

    def integrator(self, target, its) :
        if len(its) != 1 :
            print "should only see one variant at a time"
            assert False

        it = its[0]
        target["var_id"] = it[0]
        target["rsid"] = it[1]
        target["ref"] = it[2]
        target["mut"] = it[3]
        return target

def newTarget() : return {}

#assuming sources is VariantSource, DBSNPSource
#def handleAbsent( ae ) {
    #if ae.ix == 0 : pass

def flipGenotype( base ) :
    if   base == 'A' : return 'T'
    elif base == 'T' : return 'A'
    elif base == 'C' : return 'G'
    elif base == 'G' : return 'C'
    else :
        print base, "is not a valid base"
        assert False

if __name__ == '__main__' :
    conn = db.Conn("gleeson-closet")
    conn2 = db.Conn("gleeson-closet")
    query = "select id,dbSNP,ref,mut from Variants"
    for ix,row in enumerate( conn.query( query ) ) :
        var_id, rsid, ref, mut = row[0],row[1],row[2],row[3]
        if rsid :
            query = "select * from dbSNP where id = %d" % rsid
            rows = conn2.query(query)
            if len(rows) != 1 :
                print "no match for rsid: %d" % rsid
                continue

            (rsid,allele1,count1,allele2,count2) = rows[0]
            if allele2 :
                total = float(count1 + count2)
                if   ref == allele1 and mut == allele2 :
                    freq = count2 / total
                elif ref == allele2 and mut == allele1 :
                    freq = count1 / total
                else :
                    ref = flipGenotype(ref)
                    mut = flipGenotype(mut)
                    if   ref == allele1 and mut == allele2 :
                        freq = count2 / total
                    elif ref == allele2 and mut == allele1 :
                        freq = count1 / total
                    else :
                        print row, rows[0]
                        print "no match, even after flip\n"
                        #assert False
            else :
                if mut == allele1 :
                    freq = 1
                else :
                    mut = flipGenotype(mut)
                    if mut == allele1 :
                        freq = 1
                    else:
                        print row, rows[0]
                        print "no match\n"
                        #assert False

                #print "ref,mut", ref, mut, "1,2:", allele1, allele2 


        if ix > 500 : break


    #dbss = DBSNPSource("SNPAlleleFreq.bcp","Allele.bcp")
#
    #name = "dbSNP.txt"
    #fout = open(name,'wb')
    #fbulk = csv.writer( fout, \
                        #delimiter=',', \
                        #quoting=csv.QUOTE_MINIMAL )
    #for d in dbss.iterate() :
        #fbulk.writerow( d )
#
    #fout.close()

#LOAD DATA LOCAL INFILE 'dbSNP.txt' INTO TABLE dbSNP FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n' (id,allele1,count1,allele2,count2);

