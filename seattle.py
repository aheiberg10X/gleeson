import globes
import indel
import broad
from annotator import Annotator,doNothing
import db
import importer
from variant import Isoform

#######################################################################
##############    SeattleSeq   ########################################
#######################################################################
cols = ["# inDBSNPOrNot","chromosome","position","referenceBase","sampleGenotype","sampleAlleles","allelesDBSNP","accession","functionGVS","functionDBSNP","rsID","aminoAcids","proteinPosition","cDNAPosition","polyPhen","granthamScore","scorePhastCons","consScoreGERP","chimpAllele","CNV","geneList","AfricanHapMapFreq","EuropeanHapMapFreq","AsianHapMapFreq","hasGenotypes","dbSNPValidation","repeatMasker","tandemRepeat","clinicalAssociation","distanceToSplice","microRNAs","proteinSequence"]

indexOf = {}
for i,col in enumerate(cols) :
    indexOf[col] = i

class SeattleAnnotator(Annotator) :
    def __init__(self) :
        Annotator.__init__(self,'seattle')
        self.indexOf = indexOf
        self.cols = ["accession", "functionGVS", "polyPhen", \
                     "granthamScore", "scorePhastCons", "consScoreGERP", \
                     "distanceToSplice", "cDNAPosition", "AfricanHapMapFreq", \
                     "EuropeanHapMapFreq", "AsianHapMapFreq"]

        self.db_cols = ["accession","ss_functionGVS","ss_polyPhen",\
                        "ss_granthamScore","ss_scorePhastCons",\
                        "ss_consScoreGERP","ss_distanceToSplice",\
                        "ss_cDNAPosition", \
                        "ss_AfricanHapMapFreq","ss_EuropeanHapMapFreq",\
                        "ss_AsianHapMapFreq"]

    def parse( self, varList ) : pass

    def run( self ) : pass

    def getPosition( self, out_splt ) :
        if self.switch == 'snp' :
            keys = ["chromosome","position","referenceBase","sampleAlleles"]
            (chrom2,pos2,ref2,als) = [out_splt[self.indexOf[k]] for k in keys]
            sp = als.split('/')
            if len(sp) == 2 :
                (a1,a2) = sp
                if a1 == ref2 : mut2 = a2
                elif a2 == ref2 : mut2 = a1
                else : assert "Have a problem" == "with figuring out mut2"
            elif len(sp) == 1 :
                mut2 = sp[0]
            else :
                assert "length of sampleAllelels" == "not == 2 or 1"


        elif self.switch == 'indel' :
            keys = ['chromosome','position','referenceBase','sampleGenotype']
            (chrom2,pos2,ref2,sg) = [out_splt[self.indexOf[k]] for k in keys]
            #print "readiung: ", chrom2, pos2, ref2, sg
            samples = sg.split(',')
            mut2 = ""

            isInsertion = '-' in ref2
            if isInsertion :
                ref2 = ref2[0]
                for sample in samples :
                    (one,two) = sample.split('/')
                    if ref2 != one :
                        mut2 = one
                        break
                    elif ref2 != two :
                        mut2 = two
                        break
                    else : continue
                if mut2 == 'N' or mut2 == '' : mut2 = '*'
            else :
                for sample in samples :
                    (one,two) = sample.split('/')
                    if one[1:] == ref2 :
                        ref2 = one
                        mut2 = one[0]
                        break
                    elif two[1:] == ref2 :
                        ref2 = two
                        mut2 = two[0]
                        break
                    else : continue
                if mut2 == "" :
                    ref2 = '*'
                    mut2 = '*'

        else : assert 'switch must be' == 'snp or indel'

        #not ideal, but if SeattleSeq is going to be a bitch
        #it's the best we can do
        if mut2 == 'N' : mut2 = '*'
        #print "returning: ", chrom2, pos2, ref2, mut2
        return (chrom2,pos2,ref2,mut2)


    def sqlComparator( self, sqlrow, out_splt ) :
        (chrom1,pos1,ref1,mut1) = [sqlrow[i] for i in range(1,5)]
        (chrom2,pos2,ref2,mut2) = self.getPosition( out_splt )
        return globes.compareVariants( chrom1,pos1,ref1,mut1,\
                                       chrom2,pos2,ref2,mut2 )

    def nullify( self, value ) :
        if value == 'NA' or value == 'unknown' : 
            return ''
        else : return value

    @db.catch
    def sqlIntegrator( self, sqlrow, out_splt ) :
        (eyeD,chrom,pos,ref,mut,gene_id) = sqlrow
        values = [ self.nullify( out_splt[self.indexOf[c]] ) for c in self.cols]

        #self.conn.update("Genes", [values[0]], [self.db_cols[0]], gene_id)
        self.conn.update("Variants", values[1:], self.db_cols[1:], eyeD)


    def varListComparator( self, variant, out_splt ) :
        #print "is seattle's varListComparator even being called?"
        (chrom1,pos1,ref1,mut1) = variant.getPosition()
        (chrom2,pos2,ref2,mut2) = self.getPosition( out_splt )
        #print chrom1,pos1,ref1,mut1,chrom2,pos2,ref2,mut2
        return globes.compareVariants( chrom1,pos1,ref1,mut1, \
                                       chrom2,pos2,ref2,mut2 )

    def varListIntegrator( self, variant, out_splt ) :
        pp = out_splt[indexOf["proteinPosition"]]
        pos,tot = -1,-1
        if not pp == '' and not pp == 'NA':
            splt = pp.split('/')
            pos,tot = int(splt[0]), int(splt[1])

        accession = out_splt[indexOf["accession"]]
        unmatched_gene_ids = []
        message = []
        if accession and not accession == 'none' :
            query = "select id from Genes where refseq = '%s'" % accession
            gene_ids = self.conn.query( query )
            if not gene_ids :
                gid = importer.makeEmptyGene( self.conn, 'refseq', accession )
                unmatched_gene_ids.append( gid )
                #create a gene, append to unmatched_gene_ids
                message.append("Accession: %s is missing from Genes table, created a new Gene entry for it, id: %d" % (accession,gid))

            for gene_id in [int(row[0]) for row in gene_ids] :
                gene_id = int(gene_id)
                match = False
                for iso in variant.isoforms :
                    k = "gene_id"
                    if k in iso.fields :
                        message.append("iso_gene_id: %d" % iso.fields[k] )
                    if k in iso.fields and iso.fields[k] == gene_id :
                        iso.fields["ss_functionGVS"] = out_splt[indexOf["functionGVS"]]
                        iso.fields["ss_polyPhen"] = out_splt[indexOf["polyPhen"] ]
                        match = True
                if not match :
                    for iso in variant.isoforms :

                        if iso.getFields(["codon_pos"])[0] == pos and \
                           iso.getFields(["codon_total"])[0] == tot : 
                            message.append("Accession: %s has gene_ids, this one: %s doesn't match %s but does match the codon_pos: %d and codon_total: %d. This situtation current treated as unmatched." % (accession,gene_id,str(iso),pos,tot) )

                    unmatched_gene_ids.append( gene_id )
        else :
            unmatched_gene_ids.append(-1)

        for gene_id in unmatched_gene_ids :
            iso = Isoform()
            iso.fields["gene_id"] = gene_id
            dbcols = ["ss_functionGVS","ss_polyPhen", \
                      "codon_pos","codon_total"]

            iso.fields["ss_functionGVS"] = out_splt[indexOf["functionGVS"]]
            iso.fields["ss_polyPhen"] = out_splt[indexOf["polyPhen"]]
            iso.fields["codon_pos"] = pos
            iso.fields["codon_total"] = tot
            variant.isoforms.append(iso)


        for (c,dbc) in zip(self.cols,self.db_cols) :
            variant.fields[dbc] = self.nullify( out_splt[indexOf[c]] )

        return "| \t |".join(message)

    def stopper( self, splt ) :
        return len(splt) == 1

    def register( self, dargs ) :
        self.switch = dargs['switch'].lower()
        self.iterator = globes.splitIterator( dargs["file"], \
                                              burn=1, \
                                              stopper=self.stopper )

        #the default interaction will be with a varList
        self.comp = self.varListComparator
        self.eqfunc = self.varListIntegrator
        self.ltfunc = doNothing
        self.gtfunc = doNothing
        self.allow_unmatched = True
        self.conn = dargs["dbconn"]

        k = 'target'
        if k in dargs and dargs[k] == "sql" :
            self.comp = self.sqlComparator
            self.eqfunc = self.sqlIntegrator

#TODO: finish
def filterGERP( SNPList, seattle_seq_annotations_filename, opts=['n',-10] ) :
    fin = open( seattle_seq_annotations_filename, 'rb' )
    header = fin.readline()
    newSNPList = []

    #options
    exclude_slice_utr_snps = (opts[0].lower() == 'y')
    GERP_lower_thresh = opts[1]

    indexOf = {"inDBSNP" :     0,
               "chrom" :       1,
               "pos" :         2,
               "ref" :         3,
               "sampGT" :      4,
               "sampAlleles" :  5,
               "dbSNPAlleles" : 6,
               "transcript" :   7,

              }
    printColumnWarning( seattle_seq_annotations_filename, indexOf )

    #are we asserting: len(SNPList) == (num lines in file) ?
    for j in range( len(SNPList) ) :
        if entrySameAsLast( SNPList, j ) :
            pass#fill some shit in
        else :
            line = fin.readline().strip()
            splt = line.split("\t")

def uniqueWrapper( splt ) :
    calls = splt[len(broad.COLUMN_MAP)-1:]
    patients = range( len(calls) )
    (unique,message) = indel.indelUniqueToDisease( patients, calls )
    return unique

def checkCoverage( call ) :
    call_splt = broad.splitCall( call )
    if broad.isCovered( call_splt, thresh = 5 ) :
        return call
    else :
        return "./."

#take the file genreated here and put through SeattleSeq
def parseIndels( indelList ) :
    groups = {"indel_input" : globes.FAMILIES}
    out = "%s/seattle/input" % (globes.INT_DIR, )
    broad.pickOutFamilies( globes.INDEL_FILE, out, \
                           groups,\
                           lineFilter = uniqueWrapper, \
                           callToString = checkCoverage)

#take the data back from seattle seq and separate 
def separateOutputToFamilies() :
    fin = open( "%s/seattle/input/indel_input.vcf" % (globes.INT_DIR) )
    patients = broad.getPatients( fin )
    fouts = [ open("%s/indels_by_fam/%s.tsv" % (globes.OUT_DIR, \
                                                pat.replace('/','-')), 'wb' ) \
              for pat in patients ]
    fin.close()

    #errrgg so I can re-get out the original read data
    finin = open( "%s/seattle/input/indel_input.vcf" % (globes.INT_DIR) )
    patients = broad.getPatients( finin )
    finin_splt = finin.readline().strip().split('\t')

    fin = open( "%s/seattle/output/indel_output.tsv" % (globes.INT_DIR) )
    column_splt = fin.readline().strip().split('\t')
    bp = indexOf["sampleAlleles"]
    column_splt = column_splt[:bp] + ["originalBroadCall"] + column_splt[bp:]
    new_columns = "\t".join( column_splt )
    for fout in fouts :
        fout.write( "%s\n" % new_columns )

    for line in fin :
        #ignore the comment lines at the end
        if '#' in line : continue

        #get the necessary column values
        splt = line.strip().split('\t')
        cols = ["chromosome","position","refBase","sampleGenotype"]
        chrom,pos,refBase,sampleGTs = [ splt[ indexOf[c] ] for c in cols ]
        sampleGTs = sampleGTs.split(',')

        #find line in input file that corresponds to the output line
        num_incs = 0
        while True :
            cols = ["chrom","pos"]
            values = [ finin_splt[ broad.COLUMN_MAP[c] ] for c in cols ]
            finin_chrom, finin_pos = values
            finin_calls = finin_splt[ broad.COLUMN_MAP["calls"]: ]
            if pos == finin_pos and chrom == finin_chrom :
                break
            else :
                num_incs += 1
                finin_splt = finin.readline().strip().split('\t')

        #The output may have multiple lines for each input line, corresponding
        #to the different transcripts. This means that if 'line' no longer
        #matches the finin_line, we should only have to jump next once
        assert num_incs <= 1

        #isMutated is a function that takes a GT from the output file
        #and determines if it is a mutation
        isInsertion = '-' in refBase
        if isInsertion :
            l,r = refBase.split('-')
            isMutated = lambda gt : \
                            not gt == '%s/%s' % (l,l) and not gt == "N/N"
        else :
            isMutated = lambda gt : \
                            not refBase in gt.split('/')[1] and not gt == 'N/N'

        num_mutations = 0
        for i,(fout,gt) in enumerate( zip(fouts,sampleGTs) ) :
            if isMutated(gt) :
                num_mutations += 1
                #splt_copy = list(splt)
                #print splt_copy
                splt[ indexOf["sampleGenotype"] ] = "%s" % (gt)
                newline = "\t".join( splt[:bp] + [ finin_calls[i] ] + splt[bp:] )
                fout.write( "%s\n" % newline )

        #because of indel.indelUniqueToDisease
        assert num_mutations == 1

    [f.close() for f in fouts]
    fin.close()

if __name__ == '__main__' :
    #indelList = indel.inputINDELS( globes.INDEL_FILE )
    #parseIndels( 5 )
    separateOutputToFamilies()
