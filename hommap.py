import broad
import globes
import plateII

indexOf = broad.COLUMN_MAP
callIndexOf = broad.CALL_MAP

def formatGenotype( call, coverage_thresh = 8 ) :
    call_splt = broad.splitCall( call )
    #return str(broad.convertGT( call_splt ))

    gt = call_splt[ callIndexOf["GT"] ]
    if gt == './.' : return "--"
    else :
        dp = int( call_splt[ callIndexOf["DP"] ] )
        if dp >= coverage_thresh :
            if   gt == "0/0" : return "AA"
            elif gt == "0/1" : return "AB"
            elif gt == "1/1" : return "BB"
            else : raise Exception("Genotype not recognized")
        else : return "--"

def existsINdbSNP( line_splt ) :
    dbsnp = line_splt[ indexOf["dbSNP"] ]
    passes = line_splt[ indexOf["filter"] ]
    return not dbsnp == '.' and passes == 'PASS'

def makeInput( vcf_filename ) :
    header = "SNP*ID\t%s" % "\t".join( globes.FAMILIES )
    outdir = "%s/hommap" % globes.INT_DIR
    outname = "errybody"
    groups = {outname : globes.FAMILIES}
    broad.pickOutFamilies( vcf_filename, outdir, \
                           family_groups = groups, \
                           callToString = formatGenotype, \
                           shouldPrintLineSplt = existsINdbSNP, \
                           cols_to_use = [indexOf["dbSNP"]] )

    #This is extremely hacky; the whole purpose is to split
    #lines with multiple dbSNPs (sep by ';') into two lines
    fin = open("%s/%s.vcf" % (outdir,outname), 'rb')
    fout = open("%s/%s_done.vcf" % (outdir,outname), 'wb')
    (columns, headers) = broad.getColumnsAndHeaders( fin, "009-MTI" )
    columns[0] = "SNP*ID"
    fout.write( "%s\n" % "\t".join(columns) )
    for dataline in fin :
        splt = dataline.strip().split('\t',1)
        rss = splt[0].split(';')
        for rs in rss :
            fout.write( "%s\t%s\n" % (rs,splt[1]) )
    fin.close()
    fout.close()

if __name__ == '__main__' :
    makeInput( globes.SNP_FILE )
