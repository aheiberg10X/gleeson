import sys
import os

def readInChunks(file_object, MB=64):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 32MB."""
    chunk_size = 1024 #1KB
    end_reached = False
    while not end_reached :
        chunk_list = []
        for i in range( MB*1024 ) :
            data = file_object.read(chunk_size)
            if not data :
                end_reached = True
                break
            else : chunk_list.append( data )
        yield "".join(chunk_list)


def splitVCF( vcf_filename ) :
    [path, file] = vcf_filename.rsplit("/",1)
    print path, " | ", file
    [name, ext] = file.split(".")

    fin = open( vcf_filename, 'rb' )
    subdir =  "%s/%s" % (path,name)
    if not os.path.exists( subdir ) :
        os.mkdir( subdir )

    header_lines = []
    while True :
        line = fin.readline()
        lineIsHeader = line.find( "CHROM" ) == -1
        if not lineIsHeader : 
            break
        else :
            header_lines.append( line )

    #append the line that caused the 'break'
    header_lines.append( line )
    header_string = "".join(header_lines)

    chunk_num = 1
    for chunk_string in readInChunks( fin ) :
        print "Chunk num: %d" % chunk_num
        chunk_filename = "%s/chunk_%d.%s" % (subdir,chunk_num,ext)
        fchunk = open( chunk_filename, "wb" )
        fchunk.write( header_string )
        fchunk.write( chunk_string )
        chunk_num += 1


    #i = 0
    #chunklines = 10000
    #for line in fin.readlines() :
        #if i % chunklines == i :
            #print "starting new chunk"
            #if fchunk : fchunk.close()
       # 
        #fchunk.write( line )
        #i += 1

    fin.close()

def writeChunk( chunk_filename, merged_file ) :
    fchunk = open( chunk_filename, 'rb' )
    merged_file.write( fchunk.read() )
    fchunk.close()

def combineResults( jobdir, chunk_count ) :

    ssdir = "../seattleseq"
    fmerged = open("%s/%s/merged.txt" % (ssdir,jobdir),'wb')
    chunk = 1
    files = [f for f in os.listdir( "%s/%s" % (ssdir,jobdir) ) if not f.endswith( ".gz" )]
    while chunk <= chunk_count :
        prefix = "SeattleSeqAnnotation.chunk_%d.vcf" % chunk
        [writeChunk( "%s/%s/%s" % (ssdir,jobdir,t),fmerged ) for t in files if t.startswith(prefix)]
        chunk += 1

    fmerged.close()

if __name__ == "__main__" :

    argc = len(sys.argv) 
    if argc == 1 :
        splitVCF( "../data/Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf" )
    elif argc == 2 :
        vcf_filename = sys.argv[1].strip()
        splitVCF( vcf_filename )
    else :
        print "argc > 2, stopping"

    #combineResults("0526",17)
