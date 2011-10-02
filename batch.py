import globes
import re
from subprocess import Popen

class Batch :
    def __init__(self, cshfile) :
        self.file = cshfile
    def submit(self) :
        pop = Popen(['qsub',self.file])
        pop.wait()
        print "Job submitted"

def writeBatchFile( name, commands ) :
    header = \
'''#!/bin/csh
#PBS -q small 
#PBS -N sift
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -o %s/output/%s_out.txt
#PBS -e %s/output/%s_err.txt
#PBS -V
#PBS -M andrew.heiberg@gmail.com
#PBS -m abe
#PBS -A gleeson-lab''' % ( globes.BATCH_DIR, name, globes.BATCH_DIR, name )

    path = "%s/%s.csh" % (globes.BATCH_DIR, name)
    fout = open(path, 'wb' )
    fout.write( "%s\n\n" % header )
    commands = '\n\n'.join( commands )
    oneslash = re.compile(r'/+')
    commands = oneslash.sub( '/', commands )
    fout.write( commands )
    fout.write( "\n\n" )
    fout.close()
    return Batch(path)


if __name__ == '__main__' :
    #writeBatchFile('test', commands=['cd /home/aheiberg','ls'] )

    writeBatchFile('snpeff', commands=['cd /projects/gleeson-lab/bin/snpEff_v1_9_5' ,'java -Xmx20g -jar snpEff.jar -f hg37 /projects/gleeson-lab/BroadData/PlateII/raw_data/Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf'] )
