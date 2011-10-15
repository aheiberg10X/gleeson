class Pilot :
    def __init__(self) :
        self.folder = "Pilot"
        self.snpfile = "pilot_full_snps.vcf"
        self.indelfile = "pilot_full_indels.vcf"
        self.broadfile = "doesnt.exist"
        self.families = []

class PlateII :
    def __init__(self) :
        self.folder = "PlateII"
        self.snpfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf"
        self.indelfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_INDELS.vcf"
        self.broadfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526.vcf"
        self.families = []

class PlateI :
    def __init__(self) :
        self.folder = "PlateI"
        self.snpfile = "plateI_inhouse_hg19_snps.vcf"
        self.indelfile = "plateI_inhouse_hg19_indels.vcf"
        self.broadfile = "doesn't exist"
        self.families= []

class PlateIII :
    def __init__(self) :
        self.folder = "PlateIII"
        self.snpfile = "plateIII_broad_hg19_snps.vcf"
        self.indelfile = "plateIII_broad_hg19_indels.vcf"
        self.broadfile = "doesn't exist"
        self.families= []

class CIDR :
    def __init__(self) :
        self.folder = "CIDR"
        self.snpfile = "allcidrsnps.vcf"
        self.indelfile = "allcidrindels.vcf"
        self.broadfile = "doesn't exist"
        self.families= []

class Frazer :
    def __init__(self) :
        self.snpfile = "cbhag_allsnps.vcf"
        self.folder = "Frazer"
        self.indelfile = "42"
        self.broadfile = "42"
        self.families = []
