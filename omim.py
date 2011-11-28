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

if __name__ == '__main__' :
    #makeGeneIDs()
    fillGenesTable()

#just grab the 'FIELD TI' lines out of omim.txt
cat omim.txt | awk '/\*FIELD\* TI/ {getline; print $0 }' | gawk 'match($0, /([0-9]+.*)/, ary) {print ary[1]}' > omim_field_ti.txt
