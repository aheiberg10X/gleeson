import globes
import broad
from annotator import Annotator, doNothing
import batch
import snp

import os
import subprocess

class SIFTAnnotator(Annotator) :
    #switch is redundant here,  can deduce from var_list type
    def __init__(self, switch) :
        self.switch = switch.lower()
        Annotator.__init__(self, "sift" )
        self.input_file = "%s/input/%s.csv" % (self.root_dir,switch)
        self.output_dir = "%s/output" % (self.root_dir)
        self.indexOf = {"coordinates" : 0, \
                        "codons" : 1, \
                        "transcript id" : 2, \
                        "protein id" : 3, \
                        "substitutions" : 4, \
                        "region" : 5, \
                        "dbsnp id" : 6, \
                        "snp type" : 7, \
                        "prediction" : 8, \
                        "score" : 9, \
                        "median info" : 10, \
                        "# seqs at position" : 11, \
                        "gene id" : 12, \
                        "gene name" : 13, \
                        "gene desc" : 14, \
                        "omim disease" : 15, \
                        "average allele freqs" : 16, \
                        "ceu allele freqs" : 17, \
                        "user comment" : 18}
        
        self.cols = ["region","snp type", "prediction","score","omim disease"]
        self.db_cols = ["sift_region","sift_type","sift_prediction","sift_score","sift_omim"]

    def parse( self, varList ) :
        if self.switch == 'snp' or self.switch == 'indel' :
            if os.path.isfile( self.input_file ) :
                print "There is already an input file for SIFT @ %s. ?Date?" + \
                      " This (probably) means SIFT has already be run" + \
                      " and these existing annotations can be merged." + \
                      " Would you like to skip to the merge step (s)" + \
                      " or re-run SIFT (r) ?"
                uin = raw_input()
                if uin.lower() == 's' :
                    return True

            fout = open( self.input_file, 'wb' )
            if self.switch == 'snp' :
                for snp in varList.iterate() :
                    keys = ["chrom","pos","ref","mut","strand"]
                    values = [snp.fields[k] for k in keys]
                    (chrom,pos,ref,mut,strand) = values
                    if strand : strand = 1
                    else : strand = -1
                    fout.write( "%s,%s,%d,%s/%s\n" % (chrom,pos,strand,ref,mut) )
            else :
                for var in varList.iterate() :
                    keys = ["chrom","pos","ref","mut","strand"]
                    values = [var.fields[k] for k in keys]
                    print values
                    (chrom,pos,ref,mut,strand) = values

                    if strand : strand = 1
                    else : strand = -1
                    isInsertion = len(ref) == 1 and len(mut) > 1
                    isDeletion = len(ref) > 1 and len(mut) == 1
                    if isInsertion :
                        start = int(pos)
                        end = start
                        allele = mut
                    elif isDeletion :
                        start = int(pos)
                        end = start + (len(ref)-len(mut))
                        allele = '-/'
                    else : assert False
                    fout.write("%s,%d,%d,%d,%s\n" % (chrom,start,end,strand,allele) )
            fout.close()

        else : raise Exception("switch should be 'snp' or 'indel'")

        return False


    ##SOMETHING IS EVER SO SLIGHTLY OFF IN TRITON BATCH SCRIPT
    def run(self, i="", o="", args={'run_locally':True}) :
        if not i : i = self.input_file
        if not o : o = self.output_dir
        run_locally = args['run_locally']
        #-i : input
        #-o : output
        #-d : variation db directory
        #-c : coding info directory
        snp_exec = "./SIFT_exome_nssnvs.pl" #% (globes.SIFT_BIN)
        snp_exec_and_args = [ snp_exec, \
                     '-i', i, \
                     '-o', o, \
                     '-d', \
                     '%s/db/Human_db_37' % (globes.SIFT_HOME), \
                     '-A', '1', \
                     '-B', '1', \
                     '-C', '1', \
                     '-J', '1', \
                     '-K', '1', \
                     '-L', '1']

        indel_exec = ['perl', "%s/SIFT_exome_indels.pl" % (globes.SIFT_BIN)]
        indel_exec_and_args = \
            indel_exec + \
             ['-i', i, \
             '-o', o, \
             '-c', \
             '%s/coding_info/Homo_sapien_37' % (globes.SIFT_HOME), \
             '-d', \
             '%s/db/Human_db_37' % (globes.SIFT_HOME) ]

        #choose which version to run
        if self.switch == 'snp' :
            exec_and_args = snp_exec_and_args
        else :
            exec_and_args = indel_exec_and_args

        if run_locally :
            pop = subprocess.Popen(
                     exec_and_args, \
                     stdout = open( "%s/stdout" % globes.SIFT_OUTPUT, 'wb' ), \
                     stderr = open( "%s/stderr"% globes.SIFT_OUTPUT, 'wb' ) )

            pop.wait()
        else :
            params = {'out' : 'sift_out.txt' , \
                   'err' : 'sift_err.txt' }
            commands = ['cd %s' % globes.SIFT_BIN, ' '.join( exec_and_args )]
            batchJob = batch.writeBatchFile( "sift_%s.csh" % self.switch, commands )
            batchJob.submit()

            #check output file by getting jobid, using that to check outpage
            #for status 'complete'.  Then can check every X seconds to move on once complete.  OR, just leave this up to user, and separate merge as standalone to be invoked once output is generated, will have to do this for SS anyway

    def clean( self, value ) :
        if value == 'N/A' : return ''
        elif value.startswith('DAMAGING') : return 'DAMAGING'

    def varListComparator( self, snp, out_splt ) :
        (chrom1,pos1,ref1,mut1) = snp.getPosition()
        coord = out_splt[ self.indexOf["coordinates"] ].split(',')
        (chrom2,pos2,strand,refmut) = coord
        (ref2,mut2) = refmut.split('/')
        return globes.comparePositions( chrom1,pos1,chrom2,pos2 )
        #return globes.compareVariants( chrom1,pos1,ref1,mut1,\
                                        #chrom2,pos2,ref2,mut2 )

    def varListIntegrator( self, snp, out_splt ) :
        for (c,dbc) in zip(self.cols,self.db_cols) :
            snp.fields[dbc] = self.clean( out_splt[self.indexOf[c]] )

    def sqlComparator( self, sqlrow, out_splt ) :
        (chrom1,pos1,ref1,mut1) = sqlrow[1], sqlrow[2],sqlrow[3],sqlrow[4]
        coord = out_splt[ self.indexOf["coordinates"] ].split(',')
        (chrom2,pos2,strand,refmut) = coord
        (ref2,mut2) = refmut.split('/')
        return globes.compareVariants(  chrom1,pos1,ref1,mut1,\
                                        chrom2,pos2,ref2,mut2 )

    def sqlIntegrator( self, sqlrow, out_splt ) :
        (eyeD,chrom,pos,ref,mut,gene_id) = sqlrow
        ##update the gene with any information
        #q = "select name from Genes where id = %d" % gene_id
        #name = self.conn.query( q )[0][0]
        #if not name :
            #values = (out_splt[ self.indexOf["Gene Name"] ])
            #columns = ("name")
            #self.conn.update("Genes", values, columns, gene_id)

        #update variant
        #don't really need 'snp type' as the ref/mut AA imply it 
        values = [out_splt[ self.indexOf[c] ] for c in self.cols]
        try : float(values[3])
        except ValueError : values[3] = 'NULL'
        self.conn.update("Variants", values, self.db_cols, eyeD)
        pass

    def register( self, dargs ) :
        prediction_num = dargs["prediction"]
        predictions_file = "%s/%d/%d_predictions.tsv" \
                                % (globes.SIFT_OUTPUT, \
                                   prediction_num, \
                                   prediction_num)
        predictions_file = sort( predictions_file )
        self.iterator = globes.splitIterator( predictions_file, stopper=globes.tritonStop )

        #the default interaction will be with a varList
        self.comp = self.varListComparator
        self.eqfunc = self.varListIntegrator
        self.ltfunc = doNothing
        self.gtfunc = doNothing

        #If merging into a SQL iterator (iSQL), then will need to be passed a db 
        #connection in the arg list
        k = "dbconn"
        if k in dargs :
            self.conn = dargs[k]
            self.comp = self.sqlComparator
            self.eqfunc = self.sqlIntegrator

##################################################################
#################    SIFT     ####################################
##################################################################


#takes in a line, returns the key that 'sortSIFT' uses to sort
def sortHelper( line ) :
    [chrom,loc] = line.strip().split('\t')[0].split(',')[0:2]
    return ( globes.chromNum(chrom), int(loc) )

def sort( sift_file ) :
    fin = open( sift_file, 'rb' )
    sorted_filename = "%s.sorted" % sift_file
    fout = open( sorted_filename, 'wb' )
    fout.write( "".join( sorted( fin.readlines()[1:], key = sortHelper ) ) )
    fout.close()
    return sorted_filename

#TODO: change this to accept SNPList, not file
def parseForSNP() :
    fin = open( globes.SNP_FILE )
    fout = open( "%s/intermediate_data/sift/input/sift_snp_input.csv" \
                    % (globes.DATA_DIR), 'wb' )

    fout.close()
    fin.close()

#TODO: change this to accept INDEL_LIST, not file
def parseForINDEL( indel_file ) :
    print indel_file
    fin = open( indel_file )
    (path,ext) = indel_file.split('.',1)
    fname = path.split('/')[-1]
    fout = open( "%s/intermediate_data/sift/input/%s_sift_input.csv" \
                    % (globes.DATA_DIR, fname), 'wb' )
    print fout

    #fast-forward through header lines
    patients = broad.getPatients( fin )

    indexOf = broad.COLUMN_MAP
    for dataline in fin :

        splt = dataline.strip().split('\t')
        col_keys = ['chrom','pos','mut','ref']
        chrom,pos,mut,ref = [ splt[ indexOf[k] ] for k in col_keys ]

        dinfo = broad.makeInfoDict( splt[ indexOf["info"] ] )
        try :
            strand = dinfo["refseq.transcriptStrand"]
        except KeyError :
            try :
                strand = dinfo["refseq.transcriptStrand_1"]
            except KeyError : #dont have it guess '+'
                strand = '+' 
                ##raise Exception("what the fuck: %s" % splt[ indexOf["info"] ] )

        if strand == "+" : strand = 1
        elif strand == "-" : strand = -1
        else : raise Exception("Strand is not + or - ??")

        isInsertion = len(ref) == 1 and len(mut) > 1
        isDeletion = len(ref) > 1 and len(mut) == 1
        if isInsertion :
            start = int(pos)
            end = start
            allele = mut
        elif isDeletion :
            start = int(pos)
            end = start + (len(ref)-len(mut))
            allele = '-/'
        else : assert False

        fout.write( "%s,%d,%d,%d,%s\n" % (chrom,start,end,strand,allele) )

    fout.close()
    fin.close()

def run( switch, input_filename, run_locally = True ) :
    if switch.lower() == "snp" :
        if run_locally :
            sift = "%s/SIFT_exome_nssnvs.pl" % (globes.SIFT_BIN)

            #-i : input
            #-o : output
            #-d : variation db directory
            #-c : coding info directory
            pop = subprocess.Popen(
                    [sift, \
                     '-i', \
                     '%s/%s.csv' % (globes.SIFT_INPUT, input_filename), \
                     '-o', \
                     '%s/sift_snp_output.csv'% (globes.SIFT_OUTPUT), \
                     '-d', \
                     '%s/db/Human_db_37' % (globes.SIFT_HOME), \
                     '-A', '1', \
                     '-B', '1', \
                     '-C', '1', \
                     '-J', '1', \
                     '-K', '1', \
                     '-L', '1'], \
                     stdout = open( "%s/stdout" % globes.SIFT_OUTPUT, 'wb' ), \
                     stderr = open( "%s/stderr"% globes.SIFT_OUTPUT, 'wb' ) )

            pop.wait()
        else :
            print "remember the filename doesn't matter, you need to be targeting the correct input file @ /home/aheiberg/andrew/batch/sift_snps.csh"
            pop = subprocess.Popen( \
                     ["qsub", "/home/aheiberg/andrew/batch/sift_snps.csh"]
                  )

    elif switch.lower() == "indel" :
        sift_exec = "%s/SIFT_exome_indels.pl" % (globes.SIFT_BIN)
        print sift_exec
        pop = subprocess.Popen(
               [sift_exec, \
                '-i', \
                '%s' % (input_filename), \
                '-o', \
                '%s' % (globes.SIFT_OUTPUT), \
                '-c', \
                    '%s/coding_info/Homo_sapien_37' % (globes.SIFT_HOME), \
                '-d', \
                    '%s/db/Human_db_37' % (globes.SIFT_HOME), \
                ], \
               stdout = open( "%s/stdout" % globes.SIFT_OUTPUT, 'wb' ), \
               stderr = open( "%s/stderr" % globes.SIFT_OUTPUT, 'wb' ))

        pop.wait()
    else :
        assert False

def annotate( SNPList ) :
    print "Inputting SIFT: %s" % sift_file
    name = globes.ROOT_DIR + globes.PLATE_NAME + \
           "intermediate_files/sift/output/%d_predictions.txt.sorted"

    #if cant_find :
        #take the steps

    fin = open( name, 'rb' )

    expected_columns = 24
    header = fin.readline()
    assert( len( header.split("\t") ) == expected_columns )

    indexOf = {"coords" :   0,
               "dbSNP" :    5,
               "SIFT" :     8,
               "SIFTProb" : 9,
               "GeneName" : 13,
               "OMIM" :     21}

    printColumnWarning( sift_file, indexOf )

    def lineAndListMatch( j ) :
        return chrom == SNPList[j].chrom and \
               int(loc) == SNPList[j].genomicLoc

    j = 0
    for line in fin.readlines() :
        splt = line.strip().split('\t')
        [chrom,line,loc,aa] = splt[ indexOf["coords"] ].split(',')

        #move down the SNPList until it matches the current line
        while lineAndListMatch( j ) : j += 1

        if j < len(SNPList) : #and lineAndListMatch(j) :
            dbSNP = splt[ indexOf["dbSNP"] ]

            #what about novel and ENSSNPxxxxx?
            new_dbSNP_info = SNPList[j].dbSNP == "." and \
                             dbsnp.startswith("rs")
            if new_dbSNP_info :
                SNPList[j].dbSNP = justDigitsFromDBSNP( dbSNP )

            SNPList[j].SIFT = splt[ indexOf["SIFT"] ]
            SNPList[j].SIFTProb = float( splt[ indexOf["SIFTProb"] ] )

            geneName = splt[ indexOf["GeneName"] ]
            OMIM = splt[ indexOf["OMIM"] ]
            if geneName :
                if geneName == SNPList[j].geneName :
                    SNPList[j].fullName = geneName
                    SNPList[j].OMIM = OMIM
                else :
                    SNPList[j].OMIM = "%s (%s?)" % (OMIM,geneName)

            j += 1

            while entrySameAsLast( SNPList, j ) :
                SNPList[j].dbSNP = SNPList[j-1].dbSNP;
                SNPList[j].SIFT = SNPList[j-1].SIFT;
                #not SIFTProb as well??
                SNPList[j].fullName = SNPList[j-1].fullName;
                SNPList[j].OMIM = OMIM
                j += 1
        else :
            break

if __name__ == '__main__' :
    #SNPList is confusingly named...
    var_list = snp.SNPList( globes.INDEL_FILE )
    anno = SIFTAnnotator('indel')
    #anno.parse(var_list)
    anno.run()
    #anno.merge( 21757 )
    #print len(anno.var_list.variants)
    #(columns,types) = anno.var_list.getColumnsAndTypes()
    #for variant in anno.var_list.variants :
        #print variant.getFields(columns)
        #print '\n'
    pass
