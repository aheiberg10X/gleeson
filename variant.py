import broad
import globes
from math import pow
#from sys import getsizeof

#for more info use broad .vcf file
#AD,"Allelic depths 
#DP,"Read Depth"
#GQ,"Genotype Quality"
#GT,"Genotype"
#PL,"Phred-scaled likelihoods"
class SQLizable :
    def getFields( self, columns ) :
        fields = []
        for c in columns :
            if c in self.fields :
                fields.append( self.fields[c] )
            else :
                fields.append("NULL")
        return fields

class BaseCall(SQLizable) :
    def __init__(self, call_splt, pat_ix, variant_ix=1) :
        self.pat_ix = pat_ix
        self.fields = {}
        self.fields["GT"] = broad.convertGT( call_splt, variant_ix )
        if len(call_splt) > 1 :
            (GT,AD,DP,GQ,PL) = call_splt
            self.fields["DP"] = int(DP)
            self.fields["GQ"] = float(GQ)

            (r,m) = [int(t) for t in AD.split(',')]
            self.fields["AD_ref"] = r
            self.fields["AD_mut"] = m

            (aa,ab,bb) = [pow(10, -1*float(t)) for t in PL.split(',')]
            self.fields["PL_AA"] = aa
            self.fields["PL_AB"] = ab
            self.fields["PL_BB"] = bb

    #def __sizeof__(self) :
        #s = getsizeof(self.fields)
        #print "base call is: ", s
        #return s

    def prettyPrint(self) :
        skeys = globes.sortKeysByValues( broad.CALL_MAP )
        t = [str(getattr(self,k)) for k in skeys]
        gt = int(t[0])
        if gt == 0 : return broad.encodeGT( gt )
        else :
            t[0] = broad.encodeGT( int(t[0]) )
            return ":".join( t )

    def isMutated(self) :
        return broad.isMutated( self.GT )

class Isoform(SQLizable) :
    def __init__(self) :
        self.fields = {}

    def clone(self) :
        new = Isoform()
        new.fields = self.fields.copy()
        return new

    #def __sizeof__(self) :
        #return getsizeof(self.fields)

    def __str__(self) :
        s = []
        for f in self.fields :
            s.append("||%s : %s||" % (f,self.fields[f]) )
        return '\n'.join(s)

class Variant(SQLizable) :
    #def __sizeof__(self) : 
        #s = getsizeof(self.fields)
        #for iso in self.isoforms :
            #s += iso.__sizeof__()
        #for c in self.base_calls :
            #s += c.__sizeof__() 
        #print "has: ", len(self.base_calls), " base calls"
        #return s

    def __init__(self) :
        self.fields = {}
        self.base_calls = []
        self.isoforms = []

    def __repr__(self) : 
        return str(self.getPosition())

    def getPosition(self):
        fields = ["chrom","pos","ref","mut"]
        return [self.fields[f] for f in fields]

#Generalizing SNPList and IndelList
class VariantList :
    def __init__(self, vcf_file, working_set_size=1000) :
        self.file = vcf_file
        fin = open( vcf_file, "rb" )
        self.patients = broad.getPatients( fin )
        self.input_iter = globes.splitIterator( fin, burn=0 )
        self.working_set_size = working_set_size
        self.working_set = []
        indexOf = broad.COLUMN_MAP
        globes.printColumnWarning( vcf_file, indexOf )
        self.annotators = []
        self.allow_unmatched = False
        self.iterator = self.iterate()
        #header lines already burnt by getPatients
        #self.variants = map( makerFunc, it )

  #def __sizeof__(self) :
        #s = 0
        #for v in self.variants :
            #s += v.__sizeof__()
        #return s

    def registerAnnotator( self, annotator, args ) :
        annotator.register(args)
        self.annotators.append( annotator )

    def refillWorkingSet(self) :
        #fill up the working set from the input file iterator 
        self.working_set = []
        goon = True
        try :
            print "Building next: %d variants" % self.working_set_size
            for i in range(self.working_set_size) :
                self.working_set.append( self.maker( self.input_iter.next() ) )
        except StopIteration :
            goon = False

        ##give each annotator a go at the working set
        #for anno in self.annotators :
            #print "%s annotator passing over working set" % anno.name
            #try :
                #anno.annotate( iter(self.working_set) )
            #except StopIteration :
                #continue
        return goon

    def iterate(self,fast_forward=0) :
        while fast_forward > 0 :
            self.input_iter.next()
            fast_forward -= 1
        goon = True
        while goon :
            goon = self.refillWorkingSet()
            print "working set refilled"
            for var in self.working_set : yield var

            #print "==========================================="
            #print "STOPPING EARLY TO DEBUG: variant.iterate()"
            #print "==========================================="
            #goon = False

    #take a line splt and turn it into a variant, overridden by SNP nad INDEL
    def maker(self) : pass


