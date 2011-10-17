import globes

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
        return 2^globes.plates[self.name]

class Pilot(Plate) :
    def __init__(self) :
        self.name = 'Pilot'
        Plate.__init__(self)
        self.snpfile = "pilot_full_snps.vcf"
        self.indelfile = "pilot_full_indels.vcf"
        self.seattle_snp_file = ""
        self.seattle_indel_file = ""
        self.broadfile = "doesnt.exist"
        self.families = []

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
        self.seattle_snp_file = ""
        self.seattle_indel_file = ""
        self.broadfile = "doesn't exist"
        self.families= []

class PlateIII(Plate) :
    def __init__(self) :
        self.name = "PlateIII"
        Plate.__init__(self)
        self.snpfile = "plateIII_broad_hg19_snps.vcf"
        self.indelfile = "plateIII_broad_hg19_indels.vcf"
        self.seattle_snp_file = ""
        self.seattle_indel_file = ""

class CIDR(Plate) :
    def __init__(self) :
        self.name = "CIDR"
        Plate.__init__(self)
        self.snpfile = "allcidrsnps.vcf"
        self.indelfile = "allcidrindels.vcf"
        self.seattle_snp_file = "SeattleSeqAnnotation131.allcidrsnps.vcf.218621577750.tsv"
        self.seattle_indel_file = "SeattleSeqAnnotation131.allcidrindels.vcf.218621621211.tsv"

class Frazer(Plate) :
    def __init__(self) :
        self.name = "Frazer"
        Plate.__init__(self)
        self.snpfile = "cbhag_allsnps.vcf"
        self.indelfile = "42"
        self.seattle_snp_file = ""
        self.seattle_indel_file = ""
