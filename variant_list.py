import broad
import globes

#Generalizing SNPList and IndelList
class VariantList :
    def build(self, vcf_file, makerFunc) :
        self.file = vcf_file

        print "Inputting SNP data\n"

        #maps 3 letter codons to 1 letter code, unused in Gleeson100.cpp
        #codonTable = getCodonTable()

        fin = open( vcf_file, "rb" )
        self.patients = broad.getPatients( fin )
        indexOf = broad.COLUMN_MAP

        globes.printColumnWarning( vcf_file, indexOf )

        self.variants = [ makerFunc(line,self.patients) \
                              for line in fin.readlines() ]

        fin.close()

    #take a file and reforge into VariantList
    def slurpVCF( file ) :
        pass
