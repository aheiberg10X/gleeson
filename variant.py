import broad
import globes
from math import pow

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
    def __init__(self, call_splt, pat_name, variant_ix=1) :
        self.pat_name = pat_name
        self.fields = {}
        self.fields["GT"] = broad.convertGT( call_splt, variant_ix )
        if len(call_splt) > 1 :
            (GT,AD,DP,GQ,PL) = call_splt
            self.fields["DP"] = int(DP)
            self.fields["GQ"] = float(GQ)

            (r,m) = [int(t) for t in AD.split(',')]
            self.fields["AD_ref"] = r
            self.fields["AD_mut"] = m

            (aa,ab,bb) = [pow(10, -1*(float(t)/10)) for t in PL.split(',')]
            self.fields["PL_AA"] = aa
            self.fields["PL_AB"] = ab
            self.fields["PL_BB"] = bb

    def prettyPrint(self) :
        print self.fields

    def isMutated(self) :
        return broad.isMutated( self.GT )


class Isoform(SQLizable) :
    def __init__(self) :
        self.fields = {}

    # not used?
    def clone(self) :
        new = Isoform()
        new.fields = self.fields.copy()
        return new

    def __str__(self) :
        s = []
        for f in self.fields :
            s.append("||%s : %s||" % (f,self.fields[f]) )
        return '\n'.join(s)

class SNPEffIsoform(SQLizable) :
    def __init__(self) :
        self.fields = {}

# A variant object is something that can be inserted into the db 
# see see insertVariant() in importer.py
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




