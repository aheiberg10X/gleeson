from annotator import Annotator, walker, doNothing
import batch
import globes
from collimator import Source
from variant import Isoform
from os.path import basename, splitext
from os import chmod
import importer
from subprocess import Popen

cols = ["chrom","position","reference","change","change type","homozygous","quality","coverage","warnings","gene_id","gene_name","bio_type","transcript_id","exon_id","exon_rank","effect","old_aa/new_aa","old_codon/new_codon","codon_num(cds)","cds_size","codons around","aas around","custom_interval_id"]

indexOf =  {}
for i,col in enumerate(cols) : indexOf[col] = i

#Chromo    Position    Reference   Change  Change type Homozygous  Quality Coverage    Warnings    Gene_ID Gene_name   Bio_type    Trancript_ID    Exon_ID Exon_Rank   Effect  old_AA/new_AA   Old_codon/New_codon Codon_Num(CDS)  Codon_Degeneracy    CDS_size    Codons around   AAs around  Custom_interval_ID

class SNPEffAnnotator(Source) :
    def __init__(self, file, switch, fast_forward=0 ) :
        self.file = file
        self.switch = switch
        self.indexOf = indexOf
        self.iterator = self.iterate(fast_forward)
        self.allow_absent = False
        self.group_repeats = True

    #iterate() helper
    def headerCheck( self, header_splt ) :
        if not cols == header_splt :
            return (False, "%s \n\ndo not match expected: \n\n%s" \
                            % (header_splt, cols))
        else :
            return (True, "good to go")

    #iterate() helper
    def stopper( self, out_splt ) :
        return True

    def skipper( self, out_splt ) :
        change_type = out_splt[ self.indexOf["change_type"] ]
        isSNP = change_type.lower() == 'snp'
        if self.switch == 'snp' :
            return not isSNP
        elif self.switch == 'indel' :
            return inSNP
        else :
            assert 'switch must be' == ' snp or indel'

    def iterate( self, fast_forward = 0 ) :
        count = 0
        for row in globes.splitIterator( \
                            self.file, \
                            burn = ?, \
                            header_line_num = ?, \
                            skipLine = self.skipper, \
                            headerSanityCheck = self.headerCheck, \
                            stopIter = self.stopper ) :
            if count < fast_forward :
                count += 1
                continue
            else : yield row

    def eqkey( self, out_splt ) :
        if self.switch == 'snp' :
            keys = ['chrom','position','reference','change']
            (chrom,pos,ref,mut) = [ out_splt[ self.indexOf[k] ] for k in keys ]
            return (globes.chromNum(chrom),pos,ref,mut)
        elif sel.switch == 'indel' :
            keys = ['chrom','position','reference','change', "change type"]
            (chrom2,pos2,ref2,mut2,change_type) = \
                              [ out_splt[ self.indexOf[k] ] for k in keys ]
            if change_type == 'INS' :
                #example:
                # 1   866511  rs60722469  C   CCCCT
                # 1   866512              *   +CCCT   INS
                return globes.compareVariants( chrom1, pos1, ref1, mut1[1:], \
                                               chrom2, int(pos2)-1, ref2, mut2[1:] )
            elif change_type == 'DEL' :
                #1   874864  .   CT  C
                #1   874865       * -T  DEL Het
                return globes.compareVariants( chrom1, pos1, ref1[1:], mut1, \
                                               chrom2, int(pos2)-1, mut2[1:], ref2 ) 
            else :
                assert "change_type %s is not" % change_type == " INS or DEL"

        else : assert 'switch must be' == ' snp or indel'

    def integrator( self, variant, out_splts ) :
        pass


    ########################################################################
    #######          OLD            ########################################
    
    
    def run(self, input_file, run_locally=True) :
        print input_file
        bname = splitext( basename(input_file) )[0]
        output_dir = "%s/%s/output" % (globes.INT_DIR, self.name)
        print output_dir
        output_file = "%s/%s.csv" \
                        % (output_dir, bname)
        if not run_locally :
            commands = ['cd /projects/gleeson-lab/bin/snpEff_v1_9_5',
                        'java -Xmx20g -jar snpEff.jar -f hg37 %s' \
                           % input_file,
                        'mv %s/output/%s_out.txt %s' \
                           % (globes.BATCH_DIR,self.name,output_file) ]
            batchJob = batch.writeBatchFile('snpeff', commands )
            batchJob.submit()
        else :
            commands = ['cd %s/bin/snpEff_v2_0_2' % (globes.ROOT_DIR), \
                        'java -Xmx5g -jar snpEff.jar hg37 %s' \
                           % input_file] #,

            filename = "run_snpeff.sh"
            shellfile = open(filename, 'wb')
            shellfile.write( '\n'.join(commands) )
            shellfile.close()
            chmod( filename, 0775 )

            pop = Popen( ['sh',  filename], \
                          stdout = open( "%s/stdout.txt" % output_dir, 'wb' ), \
                          stderr = open( "%s/std.err" % output_dir, 'wb' ) )
            pop.wait()
            print 'snpeff Popen finished'
#                        'mv %s/snpeff/output/%s_out.txt %s' \
 #                          % (globes.INT_DIR,self.name,self.output_file) ]


    def varListSNPComparator(self,variant,out_splt) :
        (chrom1,pos1,ref1,mut1) = variant.getPosition()
        keys = ['chrom','position','reference','change']
        (chrom2,pos2,ref2,mut2) = [ out_splt[ self.indexOf[k] ] for k in keys ]
        #print chrom1, pos1, ref1, mut1, chrom2, pos2, ref2, mut2
        return globes.compareVariants( chrom1,pos1,ref1,mut1,\
                                       chrom2,pos2,ref2,mut2 )

    def varListINDELComparator( self, variant, out_splt ) :
        #print 'comparing: %s to %s' % (variant.getPosition(), out_splt )
        (chrom1,pos1,ref1,mut1) = variant.getPosition()
        keys = ['chrom','position','reference','change', "change type"]
        (chrom2,pos2,ref2,mut2,change_type) = \
                          [ out_splt[ self.indexOf[k] ] for k in keys ]
        #if change_type == 'SNP' : return 1  #want to pass over
        if change_type == 'INS' :
            #example:
            # 1   866511  rs60722469  C   CCCCT
            # 1   866512              *   +CCCT   INS
            return globes.compareVariants( chrom1, pos1, ref1, mut1[1:], \
                                           chrom2, int(pos2)-1, ref2, mut2[1:] )
        elif change_type == 'DEL' :
            #1   874864  .   CT  C
            #1   874865       * -T  DEL Het
            return globes.compareVariants( chrom1, pos1, ref1[1:], mut1, \
                                           chrom2, int(pos2)-1, mut2[1:], ref2 ) 
    
    def varListIntegrator(self,variant,out_splt) :
        novel_variant = "ref_aa" not in variant.fields
        messages = []
        if novel_variant :
            aas = out_splt[ self.indexOf["old_aa/new_aa"] ].split('/')
            if len(aas) == 2 :
                variant.fields["ref_aa"] = aas[0]
                variant.fields["mut_aa"] = aas[1]


        iso_keys = ["transcript_id", \
                    "exon_rank", "effect"]
        iso = Isoform()
        for ik in iso_keys :
            iso.fields[ik] = out_splt[ self.indexOf[ik] ]

        codon_pos = out_splt[ self.indexOf["codon_num(cds)"] ]
        if codon_pos :
            iso.fields["codon_pos"] = int(codon_pos)
        #else : iso.fields["codon_pos"] = codon_pos

        cdssize = out_splt[ self.indexOf["cds_size"] ]
        dbname = "codon_total"
        if not cdssize == '' :
            iso.fields[dbname] = int(cdssize)/3
        #else : iso.fields[dbname] = cdssize

        k = "gene_id"
        ucsc_id = out_splt[self.indexOf[k]]
        if not ucsc_id == '' :
            query = "select id from Genes where ucsc_id = '%s'" \
                     % ucsc_id
            gene_ids = self.conn.query(query)
            if len(gene_ids) > 0 :
                for gene_id in [int(row[0]) for row in gene_ids] :
                    clone = iso.clone()
                    clone.fields[k] = gene_id
                    variant.isoforms.append(clone)
            else :
                gid = importer.makeEmptyGene( self.conn, 'ucsc_id', ucsc_id  )
                iso.fields[k] = gid
                messages.append( 'No gene id for ucsc_id=|%s|, create an empty gene, id: %d' % (ucsc_id,gid) )
                variant.isoforms.append(iso)
        else :
            iso.fields[k] = -1
            variant.isoforms.append(iso)


        #print "varListIntegrator, variant got: %d isoforms" % len(variant.isoforms)

        return '| \t |'.join(messages)

    def sqlComparator( self, sqlrow, out_splt ) :
        (chrom1,pos1,ref1,mut1) = sqlrow[1], sqlrow[2],sqlrow[3],sqlrow[4]
        keys = ['chrom','position','reference','change']
        (chrom2,pos2,ref2,mut2) = [ out_splt[ self.indexOf[k] ] for k in keys ]
        return globes.compareVariants(  chrom1,pos1,ref1,mut1,\
                                        chrom2,pos2,ref2,mut2 )
    ##OLD MAKE SURE TO SYNC WITH VARLIST_INT BEFOER USERING
    #def sqlIntegrator( self, sqlrow, out_splt ) :
        #(vid,chrom,pos,ref,mut,source) = sqlrow
        #variant_db_cols = ["ref_aa","mut_aa","gene_id"]
        #variant_values = []
        #aas = out_splt[ self.indexOf["old_aa/new_aa"] ].split('/')
        #if len(aas) == 2 :
            #variant_values.append(aas[0])
            #variant_values.append(aas[1])
#
        #ucsc_id = out_splt[ self.indexOf["gene_id"] ] 
        #if ucsc_id == "" :
            #variant_values.append(-1)
        #else :
            #q = "select id from Genes where ucsc_id = '%s'" % ucsc_id
            #gene_id = self.conn.queryScalar( q )
            #if gene_id :
                #variant_values.append(gene_id)
            #else :
                #variant_values.append(-2)
#
        #self.conn.update( "Variants", variant_values, variant_db_cols, vid )
#
        #query = "select count(*) from Isoforms where var_id = %d" % vid
        #count = self.conn.executeScalar( query )
        #if count == 0 :
            #isoform_db_cols = ["var_id","transcript_id", \
                        #"exon_rank","effect", "codon_num(cds)", "cds_size"]
            #isoform_values = [vid]
            #for ik in isoform_db_cols[1:-1] :
                #isoform_values.append( out_splt[ self.indexOf[ik] ] )
            #isoform_values.append( int(out_splt[self.indexOf["cds_size"]])/3 )
#
            #self.conn.insert( "Isoforms", isoform_values, isoform_db_cols )


    def skipSNPS( self, out_splt ) :
        return out_splt[indexOf['change type']] == 'SNP'

    def register( self, dargs ) :
        switch = dargs['switch'].lower()
        if switch == 'snp' : 
            skipper = globes.dontStop
        elif switch == 'indel' :
            skipper = self.skipSNPS
            self.allow_unmatched = True
        else : assert False
        self.iterator = globes.splitIterator(dargs['file'], \
                                             burn=3, \
                                             skipper=skipper, \
                                             stopper=globes.tritonStop)
        self.comp = self.varListINDELComparator
        self.eqfunc = self.varListIntegrator
        self.ltfunc = doNothing
        self.gtfunc = doNothing
        

        self.conn = dargs['dbconn']

        k = 'target'
        if k in dargs and dargs[k] == "sql" :
            self.comp = self.sqlComparator
            self.eqfunc = self.sqlIntegrator

if __name__=='__main__' :
    #effAnn = SNPEffAnnotator()
    ##effAnn.run( globes.INDEL_FILE )
    pass
