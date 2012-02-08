from collimator import Collimator, Source
from variant import Variant
import globes
from seattle import SeattleSource
from jointsnvmix import JointSNVMixSource
import db
from math import log
import broad
from plates import Plate_JSM_HCD_1577_2_1, Plate_JSM_HME_1563_2_1, Plate_JSM_HME_1565_2_4, Plate_JSM_HME_1573_2_1, Plate_JSM_HME_1574_2_2, Plate_JSM_HME_1620_2_2
import importer

dry_run = False 

plate = Plate_JSM_HME_1620_2_2()
plate_id = plate.getPlateID()
switch = 'snp'


#method needed to construct a Collimator
def comparator(a,b) :
    return globes.compareVariants( a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3] )

#method needed to construct a Collimator
def targetCreator() :
    return Variant([],[])

if __name__ == '__main__' :

    conn = db.Conn("localhost", dry_run=dry_run )

    jointSource = JointSNVMixSource( plate.varFile(switch), plate.pat_name )
    seattleSource = SeattleSource( plate.seattleFile(switch), switch )
    sources = [jointSource, seattleSource]
    coll = Collimator( sources, comparator, targetCreator )

    variant_cols = conn.getColumns("Variants")
    
    importer.setColumns("JointCalls")
    #importer.calls_table = "JointCalls"
    #importer.calls_cols = conn.getColumns(importer.calls_table)
    #importer.calls_cols_tograb = importer.calls_cols[3:]

    inserted_count = 0
    updated_count = 0

    #print "TODO: effing seattleseq isn't sorter by chrom correctly either now?"
    for i,var in enumerate(coll) :
        if i % 5000 == 0 :
            print "%d variants processed" % i

        values = var.getPosition()
        vid = importer.getVariantID( values )

        if not vid :
            inserted_count += 1
            importer.insertVariant( conn, var, plate_id, switch )
        else :
            updated_count += 1
            importer.addToCalls( conn, vid, var, plate_id )
