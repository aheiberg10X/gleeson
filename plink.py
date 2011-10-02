import globes
import os
import subprocess
import memo
import re

#level 1 indentation
def debug( s ) : globes.debug( s, 1 )

#######################################################################
############   Plink      #############################################
#######################################################################
def parse( SNPList ) :
   pass 

def run( input_file = "%s/%s/intermediate_data/plink/input/input" \
                        % (globes.DATA_DIR, globes.PLATE_NAME),
              switch='hom', args=[(1000,1)] ) :

    print "Plink run input: %s" % input_file
    if switch == 'hom' :
        old_cwd = os.getcwd()

        #why plink not being recognized as on PATH idk
        plink_exe = '/projects/gleeson-lab/plink-1.07-x86_64/plink'
        for (runkb,hets) in args:
            print "Working on %d kb run (grouped)" % runkb
            #out_dir = "%s/%s/intermediate_data/plink/output/%d_runkb_%d_hets" % \
                      #(globes.DATA_DIR, globes.PLATE_NAME, runkb, hets)
            out_dir = "%s/%s/intermediate_data/plink/output/custom" % \
                      (globes.DATA_DIR, globes.PLATE_NAME)
            if not os.path.isdir( out_dir ) : os.mkdir( out_dir )
            os.chdir( out_dir )
            pop = subprocess.Popen( [plink_exe, \
                                     '--file', input_file, \
                                     '--homozyg'], \
                                     #'--homozyg-window-kb', "%d"%runkb, \
                                     #'--homozyg-window-het', '%d'%hets, \
                                     #'--homozyg-window-missing', 3000000, \
                                     #'--homozyg-group'],\
                                   stdout = open( "%s/stdout"%out_dir, 'wb' ), \
                                   stderr = open( "%s/stderr"%out_dir, 'wb' ) )

            #pop = subprocess.Popen( [plink_exe, \
                                     #'--file', input_file, \
                                     #'--maf', '0.1', \
                                     #'--geno', '0.01', \
                                     #'--homozyg-window-snp', '5', \
                                     #'--homozyg-snp', '10', \
                                     #'--homozyg-window-kb', '500', \
                                     #'--homozyg-kb', '1000', \
                                     #'--homozyg-gap', '500', \
                                     #'--homozyg-window-missing', '1', \
                                     #'--homozyg-window-threshold', '.05', \
                                     #'--homozyg-window-het', '1', \
                                     #'--homozyg-density', '500'],
                                   #stdout = open( "%s/stdout"%out_dir, 'wb' ), \
                                   #stderr = open( "%s/stderr"%out_dir, 'wb' ) )

            pop.wait()

        print "test exited"
        os.chdir( old_cwd )

#see if a location is in a given plink file
def find( chrom, loc, plink_output ) :
    fin = open( plink_output )
    header = fin.readline()
    ws = re.compile(r'\s+')
    for line in fin.readlines() :
        splt = ws.split( line.strip() )
        c,l,r = int(splt[3]), int(splt[6]), int(splt[7])
        if c == chrom and l <= loc <= r :
            print c
            print "loc: %d is: \n\t%s" % (loc,line)
            

    print "end of file reached"

#####################################################################
########## filter and helpers  ######################################
#####################################################################

def getPatientID( line ) :
    return int( line.split('\t',1)[0][1:] )

def getInterval( line ) :
    splt = line.split('\t')
    chrom = splt[3]
    left = int(splt[6])
    right = int(splt[7])
    return (chrom,left,right)

def sortHelper( interval ) :
    return ( globes.sortableChrom( interval[0] ), interval[1], interval[2] )

#determine where the snp falls relative to the given interval
# -1: left, <
#  0: inside, ==
# +1: right, >
# Orderd by chrom then position
def whichSide( snp, interval ) :
    ichrom,ileft,iright = interval[0],interval[1],interval[2]
    diff = globes.sortableChrom( snp.chrom ) - \
           globes.sortableChrom( ichrom )

    if diff < 0 : return -1
    elif diff == 0 :
        if ileft <= snp.genomicLoc <= iright :  return 0
        elif snp.genomicLoc < ileft : return -1
        else : return 1
    else : return 1

datalines = []
def filter( patSNPs, patientID, runkb ) :
    debug("Filtering w/ run KB = %d" % runkb )
    base_dir = globes.DATA_DIR + globes.PLATE_NAME + \
               "intermediate_data/plink/output"
    fplink = open( "%s/%d_runkb/plink.hom" % (base_dir,runkb) )

    #dont' read the plink files over and over for each patient, memoize
    global datalines
    if not datalines :
        header = fplink.readline()
        datalines = fplink.readlines()

    intervals = [getInterval(l) for l in datalines \
                 if getPatientID(line) == patientID]

    intervals = sorted( intervals, key=sortHelper )

    newSNPs = []
    for snp in patSNPs :
        curint = interval.next()
        side = whichSide( snp, curint )
        if   side < 0 :  continue                  #look at next snp
        elif side == 0 : newSNPs.append( snp )     #inside, look at next snp
        else :           curint = interval.next()  #next interval

    fplink.close()


if __name__ == '__main__' :
    #run()
    #
    find( 11, 67209291, "/projects/gleeson-lab/BroadData/PlateII/intermediate_data/plink/output/custom/plink.hom")
