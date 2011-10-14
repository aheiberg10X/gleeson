import globes
import indel
import broad
from annotator import Annotator,doNothing
import db
import importer
from variant import Isoform
from collimator import Source

#######################################################################
##############    SeattleSeq   ########################################
#######################################################################
cols = ["# inDBSNPOrNot","chromosome","position","referenceBase","sampleGenotype","sampleAlleles","allelesDBSNP","accession","functionGVS","functionDBSNP","rsID","aminoAcids","proteinPosition","cDNAPosition","polyPhen","granthamScore","scorePhastCons","consScoreGERP","chimpAllele","CNV","geneList","AfricanHapMapFreq","EuropeanHapMapFreq","AsianHapMapFreq","hasGenotypes","dbSNPValidation","repeatMasker","tandemRepeat","clinicalAssociation","distanceToSplice","microRNAs","proteinSequence"]
#no cDNA
cols = ["# inDBSNPOrNot","chromosome","position","referenceBase","sampleGenotype","sampleAlleles","allelesDBSNP","accession","functionGVS","functionDBSNP","rsID","aminoAcids","proteinPosition","polyPhen","granthamScore","scorePhastCons","consScoreGERP","chimpAllele","CNV","geneList","AfricanHapMapFreq","EuropeanHapMapFreq","AsianHapMapFreq","hasGenotypes","dbSNPValidation","repeatMasker","tandemRepeat","clinicalAssociation","distanceToSplice","microRNAs","proteinSequence"]

indexOf = {}
for i,col in enumerate(cols) :
    indexOf[col] = i

class SeattleAnnotator(Source) :
    def __init__(self, file, switch) :
        self.switch = switch
        self.indexOf = indexOf
        self.iterator = globes.splitIterator( file, \
                                              burn=1, \
                                              stopper=self.stopper )

        self.allow_unmatched = False
        self.group_repeats = True

        #want this? cDNAPosition
        self.cols = ["accession", "functionGVS", "polyPhen", \
                     "granthamScore", "scorePhastCons", "consScoreGERP", \
                     "distanceToSplice", "AfricanHapMapFreq", \
                     "EuropeanHapMapFreq", "AsianHapMapFreq","clinicalAssociation"]

        self.db_cols = ["accession","ss_functionGVS","ss_polyPhen",\
                        "ss_granthamScore","ss_scorePhastCons",\
                        "ss_consScoreGERP","ss_distanceToSplice",\
                        "ss_cDNAPosition", \
                        "ss_AfricanHapMapFreq","ss_EuropeanHapMapFreq",\
                        "ss_AsianHapMapFreq"]

    def getPosition( self, out_splt ) :
        if self.switch == 'snp' :
            keys = ["chromosome","position","referenceBase","sampleGenotype"]
            (chrom2,pos2,ref2,mut2) = [out_splt[self.indexOf[k]] for k in keys]
            #keys = ["chromosome","position","referenceBase","sampleAlleles"]
            #(chrom2,pos2,ref2,als) = [out_splt[self.indexOf[k]] for k in keys]
            #sp = als.split('/')
            #print "sp:",sp
            #print "ref2: ", ref2
            #if len(sp) == 2 :
                #(a1,a2) = sp
                #if a1 == ref2 : mut2 = a2
                #elif a2 == ref2 : mut2 = a1
                #else : assert "Have a problem" == "with figuring out mut2"
            #elif len(sp) == 1 :
                #mut2 = sp[0]
            #else :
                #assert "length of sampleAllelels" == "not == 2 or 1"


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
    
    def nullify( self, value ) :
        if value == 'NA' or value == 'unknown' : 
            return ''
        else : return value

    def eqkey( self, out_splt ) :
        return self.getPosition( out_splt )

    def integrator( self, variant, out_splts ) :
        #get the variant field information from the first out_splt
        #remember there may be multiple corresponding to different isoforms
        for (c,dbc) in zip(self.cols,self.db_cols) :
            variant.fields[c] = self.nullify( out_splts[0][indexOf[c]] )

        #make variant.isoforms out of each splt
        for out_splt in out_splts :
            pp = out_splt[indexOf["proteinPosition"]]
            pos,tot = -1,-1
            if not pp == '' and not pp == 'NA':
                splt = pp.split('/')
                pos,tot = int(splt[0]), int(splt[1])

            iso = Isoform()
            iso.fields["accession"] = out_splt[indexOf["accession"]]
            iso.fields["ss_functionGVS"] = out_splt[indexOf["functionGVS"]]
            iso.fields["ss_polyPhen"] = out_splt[indexOf["polyPhen"]]
            iso.fields["codon_pos"] = pos
            iso.fields["codon_total"] = tot
            variant.isoforms.append(iso)

        return variant

    def stopper( self, splt ) :
        return len(splt) == 1


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
