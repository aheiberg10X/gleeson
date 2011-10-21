import broad
import globes
from math import pow
from collimator import Source
#from sys import getsizeof

#########SHOULD REALLY NAME THIS MODULE VCF########################3


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
        print self.fields
        #skeys = globes.sortKeysByValues( broad.CALL_MAP )
        #t = [str(getattr(self,k)) for k in skeys]
        #gt = int(t[0])
        #if gt == 0 : return broad.encodeGT( gt )
        #else :
            #t[0] = broad.encodeGT( int(t[0]) )
            #return ":".join( t )

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
    def __init__(self, fields, calls) :
        self.fields = fields
        self.base_calls = calls
        self.isoforms = []

    def __repr__(self) : 
        return str(self.getPosition()) + " and calls..."

    def getPosition(self):
        fields = ["chrom","pos","ref","mut"]
        return [self.fields[f] for f in fields]

#Generalizing SNPList and IndelList
class VariantList(Source) :
    def __init__(self, vcf_file, fast_forward=0) :
        self.indexOf = broad.COLUMN_MAP
        globes.printColumnWarning( vcf_file, self.indexOf )
        self.fin = open( vcf_file, "rb" )
        self.patients = broad.getPatients( self.fin )

        self.allow_absent = False
        self.group_repeats = False
        self.iterator = self.iterate(fast_forward)

    def iterate(self, fast_forward = 0) :
        count = 0
        for row in globes.splitIterator( self.fin, burn=0 ) :
            if count < fast_forward :
                count += 1
                continue
            else : yield row

    def eqkey(self, it) :
        fields = ["chrom","pos","ref","mut"]
        return [it[self.indexOf[f]] for f in fields]

    def integrator( self, target, splts ) :
        if len(splts) != 1 : assert "len isn't right"
        for splt in splts :
            ##TODO generalize this to make it vendor independent, call start column is feature of VCF, not broad???
            calls = splt[ broad.CALL_START: ]
            base_calls = []
            for pat_ix,c in enumerate(calls) :
                sc = broad.splitCall(c)
                gt = broad.convertGT( sc )
                if broad.isMutated( gt ) or broad.noInf( gt ) :
                    base_calls.append( BaseCall(sc,pat_ix) )

            fields = {}
            keys = broad.COLUMN_MAP.keys()
            for k in keys :
                if k == "chrom" :
                    fields[k] = globes.chromNum( splt[self.indexOf[k]] )
                elif k == "info" :
                    ##TODO generalize this to make it vendor independent
                    dinfo = broad.makeInfoDict( splt[ self.indexOf[k] ] )
                    fields["AF"] = dinfo["AF"]
                else :
                    fields[k] = splt[ self.indexOf[k] ]


            #according to: http://www.broadinstitute.org/gsa/wiki/index.php/Understanding_the_Unified_Genotyper's_VCF_files
            #ref and alt are always given for the forward strand
            fields['strand'] = True

        return Variant( fields, base_calls )




