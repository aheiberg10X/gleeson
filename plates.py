import globes

# id should match class name
fin = open( "plate_ids.txt" )
ids = {}
for line in fin.readlines() :
    (name,eyed) = line.split(":")
    ids[name] = int(eyed)

#ids =    { "Pilot" :          0, \
           #"PlateI" :         1, \
           #"PlateII" :        2, \
           #"PlateIII" :       3, \
           #"CIDR" :           4, \
           #"Frazer_ali2" :    5, \
           #"Frazer_aligned" : 6, \
           #"FrazerII" :       7, \
           #"PlateIV" :        8, \
           #"PlateIV_1" :      9, \
           #"PlateIV_2" :     10, \
           #"PlateIV_3" :     11, \
           #"PlateV_1" :      12, \
           #"PlateV_2" :      13, \
           #"PlateV_3" :      14, \
           #"PlateV_4" :      15, \
           #"Plate_JSM_HCD_1577_2_1" : 16, \
           #"Plate_JSM_HME_1563_2_1" : 17, \
           #"Plate_JSM_HME_1565_2_4" : 18, \
           #"Plate_JSM_HME_1573_2_1" : 19, \
           #"Plate_JSM_HME_1574_2_2" : 20, \
           #"Plate_JSM_HME_1620_2_2" : 21, \
           #"Plate_frazer2" : 22, \
           #"Plate_nadia" :   23, \
           #"Plate_frazer3" : 24}

# A class to abstract the notion of a plate
# Each plate must specify 5 things:
#   A folder_name (the same as the folder you put all the data in.  
#           This folder should be directly under globes.ROOT_DIR)
#   SNP and INDEL vcf file
#   Corresponding SeattleSeq SNP and INDEL file
class Plate :
    def __init__(self) :
        self.data_dir = "%s/%s" % (globes.ROOT_DIR, self.folder_name)

    def broadFile(self) :
        return "%s/%s" % (self.data_dir, self.broadfile)

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
        h = self.__class__.__name__
        return ids[h]

class Pilot(Plate) :
    def __init__(self) :
        self.folder_name = 'Pilot'
        Plate.__init__(self)
        self.snpfile = "pilot_inhouse_hg19_snps.vcf"
        self.indelfile = "pilot_inhouse_hg19_indels.vcf"
        self.seattle_snp_file = "Seattle_pilotsnps.txt"
        self.seattle_indel_file = "Seattle_pilotindels.txt"

class PlateII(Plate) :
    def __init__(self) :
        self.folder_name = "PlateII"
        Plate.__init__(self)
        self.snpfile = "plateII_snps.vcf"
        self.indelfile = "plateII_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf.gz.218613238491.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_INDELS.vcf.216474182240.tsv"
        self.families = []

class PlateI(Plate) :
    def __init__(self) :
        self.folder_name = "PlateI"
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
        self.folder_name = "PlateIII"
        Plate.__init__(self)
        self.broadfile = "Gleeson_B3_110908_105_Samples.vcf"
        self.snpfile = "plateIII_snps.vcf"
        self.indelfile = "plateIII_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.plate3_broad_hg19_snps.tar.gz.219684944227.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.plate3_broad_hg19_indels.vcf.219648601649.txt"

class CIDR(Plate) :
    def __init__(self) :
        self.folder_name = "CIDR"
        Plate.__init__(self)
        self.snpfile = "allcidrsnps.vcf"
        self.indelfile = "allcidrindels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.allcidrsnps.vcf.218621577750.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.allcidrindels.vcf.218621621211.tsv"

class Frazer_ali2(Plate) :
    def __init__(self) :
        self.folder_name = "Frazer_ali2"
        Plate.__init__(self)
        self.snpfile = "fromfrazer_threads2_snps.vcf"
        self.indelfile = "fromfrazer_threads2_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.fromfrazer_threads2_snps.vcf.221292923075.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.fromfrazer_threads2_indels.vcf.221295179021.txt"

class Frazer_aligned(Plate) :
    def __init__(self) :
        self.folder_name = "Frazer_aligned"
        Plate.__init__(self)
        self.snpfile = "fromfrazer_threads_snps.vcf"
        self.indelfile = "fromfrazer_threads_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.fromfrazer_threads_snps.vcf.221296444806.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.fromfrazer_threads_indels.vcf.221298757736.txt"

class FrazerII(Plate) :
    def __init__(self) :
        self.folder_name = "FrazerII"
        Plate.__init__(self)
        self.snpfile = "CBH-348.snps.vcf"
        self.indelfile = "CBH-348.indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.CBH-348.snps.vcf.222970028015.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.CBH-348.indels.vcf.222971913506.txt"

class PlateIV(Plate) :
    def __init__(self) :
        self.folder_name = "PlateIV"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_120211.vcf"
        self.snpfile = "Complete_Gleeson_120211.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_120211.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.plateIV.snps.vcf.223309049013.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.plateIV.indels.vcf.223306147216.txt"

class PlateIV_1(Plate) :
    def __init__(self) :
        self.folder_name = "PlateIV_1"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_295_120711_1.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_295_120711_1.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_295_120711_1.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.snps.tar.gz.223721330517.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.indels.tar.gz.223897410296.txt"

class PlateIV_2(Plate) :
    def __init__(self) :
        self.folder_name = "PlateIV_2"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_295_120711_2.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_295_120711_2.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_295_120711_2.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.snps.tar.gz.223729941739.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.indels.tar.gz.223726827483.txt"

class PlateIV_3(Plate) :
    def __init__(self) :
        self.folder_name = "PlateIV_3"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_295_120711_3.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_295_120711_3.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_295_120711_3.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.snps.tar.gz.223801286912.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation131.indels.tar.gz.223735315294.txt"

class PlateV_1(Plate) :
    def __init__(self) :
        self.folder_name = "PlateV"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_368_122211_batch_1.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_368_122211_batch_1.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_368_122211_batch_1.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_1.unannotated.vcf_snps.vcf.gz.225899199676.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_1.unannotated.vcf_indels.vcf.gz.225912947699.txt"

class PlateV_2(Plate) :
    def __init__(self) :
        self.folder_name = "PlateV"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_368_122211_batch_2.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_368_122211_batch_2.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_368_122211_batch_2.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_2.unannotated.vcf_snps.vcf.225871687747.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_2.unannotated.vcf_indels.vcf.225888136772.txt"

class PlateV_3(Plate) :
    def __init__(self) :
        self.folder_name = "PlateV"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_368_122211_batch_3.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_368_122211_batch_3.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_368_122211_batch_3.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_3.unannotated.vcf_snps.vcf.225875884690.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_3.unannotated.vcf_indels.vcf.225888255282.txt"

class PlateV_4(Plate) :
    def __init__(self) : 
        self.folder_name = "PlateV"
        Plate.__init__(self)
        self.broadfile = "Complete_Gleeson_368_122211_batch_4.unannotated.vcf"
        self.snpfile = "Complete_Gleeson_368_122211_batch_4.unannotated.vcf_snps.vcf"
        self.indelfile = "Complete_Gleeson_368_122211_batch_4.unannotated.vcf_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_4.unannotated.vcf_snps.vcf.225879594853.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation134.Complete_Gleeson_368_122211_batch_4.unannotated.vcf_indels.vcf.225891223327.txt"

class Plate_JSM_HCD_1577_2_1(Plate) :
    def __init__(self) :
        self.pat_name = "HCD-1577-2-1 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HCD-1577-2-1_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HCD-1577-2-1_SNV_result.txt.filtered.228224114054.txt.sorted"

class Plate_JSM_HME_1563_2_1(Plate) :
    def __init__(self) :
        self.pat_name = "HME-1563-2-1 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HME-1563-2-1_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HME-1563-2-1_SNV_result.txt.filtered.sorted.228553699300.txt"

class Plate_JSM_HME_1565_2_4(Plate) :
    def __init__(self) :
        self.pat_name = "HME-1565-2-4 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HME-1565-2-4_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HME-1565-2-4_SNV_result.txt.filtered.sorted.228554359976.txt"

class Plate_JSM_HME_1573_2_1(Plate) :
    def __init__(self) : 
        self.pat_name = "HME-1573-2-1 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HME-1573-2-1_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HME-1573-2-1_SNV_result.txt.filtered.sorted.228554466106.txt"

class Plate_JSM_HME_1574_2_2(Plate) :
    def __init__(self) :
        self.pat_name = "HME-1574-2-2 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HME-1574-2-2_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HME-1574-2-2_SNV_result.txt.filtered.sorted.228554659298.txt"

class Plate_JSM_HME_1620_2_2(Plate) :
    def __init__(self) :
        self.pat_name = "HME-1620-2-2 brain"
        self.folder_name = "preliminary"
        Plate.__init__(self)
        self.snpfile = "HME-1620-2-2_SNV_result.txt.filtered.sorted"
        self.seattle_snp_file = "SeattleSeqAnnotation134.HME-1620-2-2_SNV_result.txt.filtered.sorted.228648614210.txt"

class Plate_frazer2(Plate) :
    def __init__(self) :
        self.folder_name = "frazer2"
        Plate.__init__(self)
        self.broadfile =""
        self.snpfile = "from_frazer_snps.vcf"
        self.indelfile = "from_frazer_indels.vcf"
        self.seattle_snp_file = "frazer2_snps_seattle.txt"
        self.seattle_indel_file = "frazer2_indels_seattle.txt"

class Plate_nadia(Plate) :
    def __init__(self) :
        self.folder_name = "nadia"
        Plate.__init__(self)
        self.broadfile =""
        self.snpfile = "nadia_all_snps.vcf"
        self.indelfile = "nadia_all_indels.vcf"
        self.seattle_snp_file = "nadia_snps_seattle.txt"
        self.seattle_indel_file = "nadia_indels_seattle.txt"

class frazer3(Plate) :
    def __init__(self) :
        self.folder_name = "frazer3"
        Plate.__init__(self)
        self.snpfile = "from_frazer_snps.vcf"
        self.indelfile = "from_frazer_indels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation134.frazer3_snps.vcf.tar.gz.230368586074.txt"
        self.seattle_indel_file = "SeattleSeqAnnotation134.from_frazer_indels.vcf.230369948975.txt"
