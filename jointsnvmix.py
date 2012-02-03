from collimator import Source
from globes import toBaseN, splitIterator, chromNum
#TODO
#basically everything in variant has a more natural place to be
from variant import JointCall, Variant
from plates import Plate_JSM_HCD_1577_2_1

#"chrom", "position", "ref_base", "var_base"
headers = ["chrom", "pos", "ref", "mut", "normal_counts_a", "normal_counts_b", "tumour_counts_a", "tumour_counts_b", "p_AA_AA", "p_AA_AB", "p_AA_BB", "p_AB_AA", "p_AB_AB", "p_AB_BB", "p_BB_AA", "p_BB_AB", "p_BB_BB"]

indexOf = {}
for i,h in enumerate(headers) :
    indexOf[h] = i

#take the list of 9 probabilities and covert it into the correct int code
def convertGT( GT_probs ) :
    return toBaseN( GT_probs.index(max(GT_probs)), 3 )

def sortHelper( it ) :
    return (chromNum(it[0]),int(it[1]))

def sortFile( filename ) :
    #TODO complete
    fin = open(filename)
    fout = open("%s.sorted" % filename,'w')

    sorted_splts = sorted( splitIterator( fin ), key = sortHelper )
    fout.write( "\n".join( ['\t'.join(splt) for splt in sorted_splts] ) )

    fin.close()
    fout.close()

class JointSNVMixSource( Source ) :
    def __init__(self, filename, pat_name) :
        self.pat_name = pat_name
        self.fin = open(filename)
        self.allow_absent = False
        self.group_repeats = False
        self.iterator = splitIterator( self.fin )

    def eqkey(self, it) :
        return (it[0],it[1],it[2],it[3])

    def integrator(self, target, splts) :
        #this is a one to many match with SeattleSeq
        assert len(splts) == 1

        splt = splts[0]

        #transition from variant fields to JointCall fields
        ix = 4
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

if __name__ == '__main__' :
    sortFile( Plate_JSM_HCD_1577_2_1().varFile('snp') )
