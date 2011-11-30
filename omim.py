import db
import simplejson as json

def makeGeneIDs() :
    query = '''
    select id, gene
    from Isoforms
    where gene is not null
    order by gene,id'''

    conn = db.Conn("localhost")
    gmap = {}
    genes = []
    last = -1
    for row in conn.query( query ) :
        eyed = row[0]
        gene = row[1]
        if gene != last :
            genes.append( gene )
            #gmap[gene] = [eyed]
            last = gene
        #else :
            #gmap[gene].append(eyed)

    fout = open("genes.txt",'w')
    fout.write( json.dumps(genes) )
    fout.close()

def fillGenesTable() :
    conn = db.Conn("localhost",dry_run=False)
    fin = open("genes.txt",'r')
    genes = json.loads( fin.read() )
    fin.close()
    gene_id = 0
    for gene in genes :
        print gene
        conn.insert('Genes',[gene_id,gene],['id','geneSymbol'])
        up = "update Isoforms set gene_id = %d where gene = '%s'" % (gene_id,gene)
        conn.put( up )
        gene_id += 1

def updateGenesWithOMIM() :
    conn = db.Conn("localhost",dry_run=False)
    fin = open("../omim/refined.txt")
    max_len = 0
    for line in fin.readlines() :
        splt = line.strip().split('\t')
        disease_name = splt[1].split(';')[0]
        if len(disease_name) > max_len : max_len = len(disease_name)
        gene_name = splt[-1]
        query = "select id from Genes where geneSymbol = '%s'" % gene_name
        gid = conn.queryScalar( query, int )
        if gid :
            conn.update('Genes',[disease_name],['omim'],gid)
            pass
        else :
            print 'cant find gene_name %s' % gene_name
        #query = "update Genes set omim = '%s' where geneSymbol = '%s'" \
                #% (disease_name, gene_name)

    print max_len

    fin.close()

if __name__ == '__main__' :
    #makeGeneIDs()
    #fillGenesTable()
    updateGenesWithOMIM()

##Gene | Disease
#cat genemap | awk -F '|' '{ print $6 "||" $14 } '

##Below is all crap, the FIELD TI is not what we wanted


##just grab the 'FIELD TI' lines out of omim.txt, then strip of optional *|^|% char
##split mim_id and disease name with a tab
#cat omim.txt | awk '/\*FIELD\* TI/ {getline; print $0 }' | gawk 'match($0, /([0-9]+)(.*)/, ary) {print ary[1]"\t"ary[2]}' > omim_field_ti.txt

##join on the mim number of both files
#join -1 1 -2 1 omim_field_ti.txt mim2gene.txt > joined.txt

##only grab the lines that have a gene number
#then match mim_id, disease name, gene|pheno, gene_id, gene symbol
#we only care about gene symbol, as that is all we have to map on
#cat joined.txt | awk '/.*(gene|phenotype) [^-]/ { print $0 }' | gawk 'match( $0, /([0-9]+) (.*) (gene|phenotype) ([0-9]+) ([^-]+)/, ary) { print ary[1]"\t"ary[2]"\t"ary[3]"\t"ary[4]"\t"ary[5] }' > refined.txt

##Note that you can't use the gene name in the omim.txt FIELD TI field.  There are cases where it doesn't match the joined name in mim2gene (102645)
