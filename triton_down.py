#add indexOf vcf file to globes rather than retyping it for each of 

class Gene :
    def __init__(self) :
        self.exonList = []

class Exon :
    def __init__( self, start, end ) :
        self.start = start
        self.end = end

def inputGeneList( gene_list_filename ) :
    geneList = []

    fin = open( gene_list_filename )
    indexOf = {"#bin" :       0,\
               "name" :       1,\
               "chrom" :      2,\
               "strand" :     3,\
               "txStart" :    4,\
               "txEnd" :      5,\
               "cdsStart" :   6,\
               "cdsEnd" :     7,\
               "exonCount" :  8,\
               "exonStarts" : 9,\
               "exonEnds" :   10,\
               "name2" :      12}

    header = fin.readline()
    for line in fin :
        g = Gene()
        splt = line.strip().split('\t')

        #name attribute is fancy, handle it explicitly
        g.name = "%s(%s)" % (splt[indexOf["name2"]], splt[indexOf["name"]])

        #otherwise set attributes of g en masse
        attributes = ["chrom","strand","txStart","txEnd",\
                      "cdsStart","cdsEnd"]
        for attr in attributes :
            setattr( g, attr, splt[ indexOf[attr] ] )

        #turn strand attr into a boolean
        g.strand = (g.strand == '+')

        #for some reason there is a trailing ',' on starts and ends
        exonCount = int(splt[ indexOf["exonCount"] ])
        exonStarts = splt[ indexOf["exonStarts"] ].strip(',').split(',')
        exonEnds = splt[ indexOf["exonEnds"] ].strip(',').split(',')
        assert len(exonStarts) == exonCount == len(exonEnds)

        exonList = [ Exon( int(se[0]),int(se[1]) ) \
                     for se in zip(exonStarts,exonEnds)]
        g.exonList.append( exonList )
        geneList.append( g )

    fin.close()
    return geneList


#TODO: make sure samLocAsPrev returns false if j==0
def parseForSIFT( SNPList, output_filename ) :
    fout = open( input_filename, 'rb' )
    
    for j,snp in SNPList :
        if not sameLocationAsPrev( j, SNPList ) :
            string = "%s,%s,1,%s/%s\n" % \
                     (snp.chrom[3:], snp.genomicLoc, snp.refNT, snp.mutNT)
            fout.print( string )


def filterUnique( indelList, freq ) :
    newIndelList = []
    for patient in indelList :
        for indel in patient :
            


if __name__ == "__main__" :
    inputGeneList( "/home/Gleeson/Desktop/Gleeson-dir/nitin/NCBIFiles/RefSeqGeneList_hg19.txt" )

#def addData( SNPList, old_vcf_file ) :
    #fin = open( old_vcf_file )
#
    ##skip all the header files
    #while True :
        #line = fin.readline().strip()
        #lineIsHeader = line.find( "CHROM" ) == -1
        #if not lineIsHeader : break
#
    #splt = line.split('\t')
    #patients = splt[ indexOf["family_data"]: ]
#
    ##snpix
    #i = 0
    #for dataline in fin.readlines() :
        #pass
#
#def inputIndels( indels_filename ) :
    #fin = open( indels_filename, 'rb' )
    ##this is where things differ from previous
    ##nitin got/generated an individual indel.vcf file for each patient,
    ##I have one big list
#
##1st dim: patients
##2nd dim: indels
##???Frequency at which we can classify indel as homozygous???
## Allele Count / Depth
#def filterIndels( indelList, freq) :
    #print "Filtering indel frequencies"
    #newIndelsList = []
    #for i,patient in enumerate(indelList) :
        #for snp in patient :
            #if snp.ac / snp.dp >= freq :
                #newIndelsList[i].append( snp )
#
    #return newIndelList
#
#
