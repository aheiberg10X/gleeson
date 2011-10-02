import snpeff
import variant_list
import globes
import broad
import snp
import importer
from time import time
from sys import getsizeof
import cProfile
import db

def clean() :
    importer.reset()

def insert() :

    snpList = snp.SNPList( globes.SNP_FILE )
    #print "size(b): ", snpList.__sizeof__()
    effAnn = snpeff.SNPEffAnnotator( globes.SNP_FILE )
    snpList.registerAnnotator( effAnn, {} )
    effAnn.run()
    ##print snpList.variants[0].base_calls[0].fields["PL_AA"]
    ##importer.populatePatients( snpList.patients )
    #importer.insertVariantList( snpList )
    print "num homrefs not inserted:", snpList.homrefs

if __name__ == '__main__' :
    #clean()
    #cProfile.run('insert()')
    insert()

