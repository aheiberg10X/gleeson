import re
import subprocess
import os
import globes
import broad

##################################################################
#######         Globals   ########################################
##################################################################

#string identifiers of all patients
patients = []

#######################################################################
##############    dbSNP     ###########################################
#######################################################################

def dbSNPFilter( SNPList ) : pass



if __name__ == "__main__" :
    #inputData( globes.BROAD_FILE, globes.FAMILIES )
    pickOutFamilies( globes.BROAD_FILE, "spoan-dbSNP.txt", colsToUse=[2], families=["SPOAN-1513"] )
    #sortSIFT( "/projects/gleeson-lab/sift4.0.3/tmp/16382/16382_predictions.tsv" )

    #runPlink( args=[10000] )
    #plinkFind( 80874923, "/projects/gleeson-lab/BroadData/PlateI/intermediate_data/plink/output/10000_runkb/plink.hom" )

