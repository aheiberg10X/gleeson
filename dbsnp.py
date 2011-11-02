from collimator import Collimator, Source
import globes
import db

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
        #[rsid, {allele:(count,freq),allele2:(count,freq)]
        rsid_info = []
        for splt in globes.splitIterator( self.freq_file ) :
            #print splt
            rsid,alleleid,count,freq = int(splt[0]), int(splt[1]), \
                                       int(float(splt[2])), float(splt[3])

            if len(self.dallele[alleleid]) > 1 \
               or alleleid == 5 \
               or alleleid == 8 :
                continue
            if rsid != prev_rsid :
                if prev_rsid > 0 :
                    freq = 0
                    for k in rsid_info[1] :
                        freq += rsid_info[1][k][1]
                    if freq < .999 :
                        print "the way we skip indels isnt getting it right"
                        assert False

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

if __name__ == '__main__' :
    dbss = DBSNPSource("SNPAlleleFreq.bcp","Allele.bcp")
    conn = db.Conn("localhost")

    print 'wuuuuuuuuuuuut'
    fbulk = open("dbsnp_bulk_inserts.txt",'wb')
    for d in dbss.iterate() :
        fbulk.write( "%s\n" % ",".join(d) )

    fbulk.close()


    #vs = VariantSource(conn)
#
    #sources = [vs,dbss]
    #coll = Collimator( sources, globes.compareHelper, newTarget )
#
    #for t in coll :
        #if 'rsid' not in t :
            #print "no info in db about this rsid"
        #else :
            #print t
