import globes

ids =    { "Pilot" :          0, \
           "PlateI" :         1, \
           "PlateII" :        2, \
           "PlateIII" :       3, \
           "CIDR" :           4, \
           "Frazer_ali2" :    5, \
           "Frazer_aligned" : 6, \
           "FrazerII" :       7, \
           "PlateIV" :       8}

# A class to abstract the notion of a plate
# Each plate must specify 5 things:
#   A name (the same as the folder you put all the data in.  
#           This folder should be directly under globes.ROOT_DIR)
#   SNP and INDEL vcf file
#   Corresponding SeattleSeq SNP and INDEL file
class Plate :
    def __init__(self) :
        self.data_dir = "%s/%s" % (globes.ROOT_DIR, self.name)

    def varFile(self,switch) :
        if switch == "snp" :
            return "%s/%s" % (self.data_dir, self.snpfile)
        else :
            return "%s/%s" % (self.data_dir, self.indelfile)

    def seattleFile(self,switch) :
        if switch == 'snp' :
            return "%s/%s" % (self.data_dir, self.seattle_snp_file)
        else :
            return "%s/%s" % (self.data_dir, self.seattle_indel_file)

    def getPlateID( self ) :
        return ids[self.name]

class Pilot(Plate) :
    def __init__(self) :
        self.name = 'Pilot'
        Plate.__init__(self)
        self.snpfile = "pilot_inhouse_hg19_snps.vcf"
        self.indelfile = "pilot_inhouse_hg19_indels.vcf"
        self.seattle_snp_file = "Seattle_pilotsnps.txt"
        self.seattle_indel_file = "Seattle_pilotindels.txt"

class PlateII(Plate) :
    def __init__(self) :
        self.name = "PlateII"
        Plate.__init__(self)
        self.snpfile = "plateII_snps.vcf"
        self.indelfile = "plateII_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf.gz.218613238491.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_INDELS.vcf.216474182240.tsv"
        self.families = []

class PlateI(Plate) :
    def __init__(self) :
        self.name = "PlateI"
        Plate.__init__(self)
        self.snpfile = "plateI_inhouse_hg19_snps.vcf"
        self.indelfile = "plateI_inhouse_hg19_indels.vcf"
        #no cDNA column!
        self.seattle_snp_file = "SeattleSeqAnnotation131.Broad_WES_Data_GVS.218386164983_SNPS.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.plateI_inhouse_hg19_indels.vcf.219046801512.txt"
        self.broadfile = "doesn't exist"
        self.families= []

class PlateIII(Plate) :
    def __init__(self) :
        self.name = "PlateIII"
        Plate.__init__(self)
        self.snpfile = "plateIII_snps.vcf"
        self.indelfile = "plateIII_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.plate3_broad_hg19_snps.tar.gz.219684944227.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.plate3_broad_hg19_indels.vcf.219648601649.txt"

class CIDR(Plate) :
    def __init__(self) :
        self.name = "CIDR"
        Plate.__init__(self)
        self.snpfile = "allcidrsnps.vcf"
        self.indelfile = "allcidrindels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.allcidrsnps.vcf.218621577750.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.allcidrindels.vcf.218621621211.tsv"

class Frazer_ali2(Plate) :
    def __init__(self) :
        self.name = "Frazer_ali2"
        Plate.__init__(self)
        self.snpfile = "fromfrazer_threads2_snps.vcf"
        self.indelfile = "fromfrazer_threads2_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.fromfrazer_threads2_snps.vcf.221292923075.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.fromfrazer_threads2_indels.vcf.221295179021.txt"

class Frazer_aligned(Plate) :
    def __init__(self) :
        self.name = "Frazer_aligned"
        Plate.__init__(self)
        self.snpfile = "fromfrazer_threads_snps.vcf"
        self.indelfile = "fromfrazer_threads_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.fromfrazer_threads_snps.vcf.221296444806.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.fromfrazer_threads_indels.vcf.221298757736.txt"

class FrazerII(Plate) :
    def __init__(self) :
        self.name = "FrazerII"
        Plate.__init__(self)
        self.snpfile = "CBH-348.snps.vcf"
        self.indelfile = "CBH-348.indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.CBH-348.snps.vcf.222970028015.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.CBH-348.indels.vcf.222971913506.txt"

class PlateIV(Plate) :
    def __init__(self) :
        self.name = "PlateIV"
        Plate.__init__(self)
        self.snpfile = "plateIV.snps.vcf"
        self.indelfile = "plateIV.indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.plateIV.snps.vcf.223309049013.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.plateIV.indels.vcf.223306147216.txt"

