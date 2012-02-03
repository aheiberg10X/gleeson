from collimator import Collimator, Source
from variant import Variant
import globes
from seattle import SeattleSource
from jointsnvmix import JointSNVMixSource
import db
from math import log
import broad
from plates import Plate_JSM_HCD_1577_2_1

plate = Plate_JSM_HCD_1577_2_1()
switch = 'snp'

jointSource = JointSNVMixSource( plate.varFile(switch), plate.pat_name )
seattleSource = SeattleSource( plate.seattleFile(switch), switch )
sources = [jointSource, seattleSource]

#method needed to construct a Collimator
def comparator(a,b) :
    return globes.compareVariants( a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3] )

#method needed to construct a Collimator
def targetCreator() :
    return Variant([],[])

if __name__ == '__main__' :
    coll = Collimator( sources, comparator, targetCreator )
    print "TODO: effing seattleseq isn't sorter by chrom correctly either now?"
    #for i,var in enumerate(coll) :
        #print var
