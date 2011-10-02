#############################################################################
###################   Plate II    ###########################################
#############################################################################
class PlateII :
    def __init__(self) :
        self.folder = "PlateII"
        self.snpfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf"
        self.indelfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_INDELS.vcf"
        self.broadfile = "Ciliopathies_Whole_Exome_Gleeson_Indels_20110526.vcf"
        self.families = all_wes_families

###########################
dom = ["AAS-173-III-19","AAS-173-IV-28","FPKD-1482-2-1","FPKD-1482-2-4","MIC-1496-III-4-1"]

rec = ["MTI-231-3-1","CVH-468-I-3-4","CCH-523-4-5","MTI-578-3-4/MTI-1015 ","MTI-649-3-4","CHIME-681-3-4","ACC-590-4-2","LIS-PMG-711-II-4-4","PCH-805-4-2","PCH-819-4-1","HSP-860-4-1","HSP-889-4-2","MSGP-917-4-2","LIS-920-5-3","MRE-938-4-7","EPI-968-I-4-4","MIC-1079-3-2","CCH-1229-4-5","CCH-1234-3-4","PCH-1236-3-2","Ataxia-1239-3-2","Neurop-1243-4-5","CBA-1244-3-2","Hypertonia-1247-3-2","MIC-1250-3-1","MIC-1252-3-2","CCH-1253-3-2","MIC-CCH-1257-3-4","CCH-MIC-1258-5-1","PMG-1260-3-1","MIC-1262-3-3","MIC-1265-6-1","MRE-1266-4-1","MIC-1269-4-1","MRE-1277-5-1","MIC-1279-4-1","MSGP-1281-4-1","MTI-1283-5-5","MIC/PCH-1287-5-1","BFPP-1289-I-4-2","HSP-1290-5-2","PCH-1292-2-2","PCH-1298-4-5","MIC-1300-3-10","EPI-1317-3-1","HSP-1334-4-6","HSP-1349-3-1","MTI-1352-3-3","Dysmorph-1360-4-2","MTI-1367-II-5-2","HSP-1370-III-4-4","CBH/CBA-1373-3-1","AEM-1374-4-4","CCH-1376-4-6","CBA-1377-3-2","MSGP-1383-2-3","MIC-1385-3-1","PCH-1389-5-3","PCH-1391-5-4","MIC-CCH-1392-4-1","HSP-1393-II-4-1","MRPN-1397-4-2","CCH-1405-3-4","CBA-CCH-1410-3-1","AGS-1411-4-2","AGS-1412-2-1","MSGP-1413-2-6","MSGP-1413-2-7","CCH-1414-3-12","MIC-1415-4-4","CCH-MTI-1420-4-2","WMD-1422-4-2","MIC-1426-3-1","MIC-1430-I-4-4","MIC-CVH-1431-I-5-2","MIC-1435-3-8","DWM-1454-3-8","MIC-1456-II-3-4","CBA-CBH-1486-3-1","MTI-1492-3-1","MIC-1496-III-4-1","MIC-1499-II-5-2","MRE-1509-3-3","MR-1510-5-2","MIC-1511-3-1","DRV-1512-1-1","DRV-1512-1-2","DRV-1512-2-1","SPOAN-1513","ALGS-1527-3-1","EPI-1317-3-3","0131-MTI","0242-MTI","0261-MTI","0265-MTI","0267-MTI","0391-MTI","0490-MTI","0492-MTI","0574-MTI","0672-MTI","0865-MTI","1039-MTI","1151-MTI","009-MTI"]

assert len(dom) + len(rec) - len( [f for f in dom if f in rec] )== 109

consang = ["CVH-468-I-3-4","CCH-523-4-5","MTI-578-3-4/MTI-1015","MTI-649-3-4","CHIME-681-3-4","ACC-590-4-2","LIS-PMG-711-II-4-4","PCH-805-4-2","PCH-819-4-1","HSP-860-4-1","HSP-889-4-2","MSGP-917-4-2","LIS-920-5-3","MRE-938-4-7","EPI-968-I-4-4","MIC-1079-3-2","CCH-1229-4-5","CCH-1234-3-4","PCH-1236-3-2","Ataxia-1239-3-2","Neurop-1243-4-5","CBA-1244-3-2","Hypertonia-1247-3-2","MIC-1250-3-1","MIC-1252-3-2","CCH-1253-3-2","MIC-CCH-1257-3-4","CCH-MIC-1258-5-1","PMG-1260-3-1","MIC-1262-3-3","MIC-1265-6-1","MRE-1266-4-1","MIC-1269-4-1","MRE-1277-5-1","MIC-1279-4-1","MSGP-1281-4-1","MTI-1283-5-5","MIC/PCH-1287-5-1","BFPP-1289-I-4-2","HSP-1290-5-2","PCH-1292-2-2","PCH-1298-4-5","MIC-1300-3-10","EPI-1317-3-1","HSP-1334-4-6","HSP-1349-3-1","MTI-1352-3-3","Dysmorph-1360-4-2","MTI-1367-II-5-2","HSP-1370-III-4-4","CBH/CBA-1373-3-1","AEM-1374-4-4","CCH-1376-4-6","CBA-1377-3-2","MSGP-1383-2-3","MIC-1385-3-1","PCH-1389-5-3","PCH-1391-5-4","MIC-CCH-1392-4-1","HSP-1393-II-4-1","MRPN-1397-4-2","CCH-1405-3-4","CBA-CCH-1410-3-1","AGS-1411-4-2","AGS-1412-2-1","CCH-1414-3-12","MIC-1415-4-4","CCH-MTI-1420-4-2","WMD-1422-4-2","MIC-1426-3-1","MIC-1430-I-4-4","MIC-CVH-1431-I-5-2","MIC-1435-3-8","DWM-1454-3-8","MIC-1456-II-3-4","MIC-1496-III-4-1","MIC-1499-II-5-2","MRE-1509-3-3","MR-1510-5-2","SPOAN-1513","ALGS-1527-3-1","EPI-1317-3-3","0131-MTI","0242-MTI","0261-MTI","0265-MTI","0267-MTI","0391-MTI","0490-MTI","0492-MTI","0672-MTI","0865-MTI","1039-MTI","1151-MTI","009-MTI"]

non_consang = ["AAS-173-III-19","AAS-173-IV-28","MSGP-1413-2-6","MSGP-1413-2-7","FPKD-1482-2-1","FPKD-1482-2-4","DRV-1512-1-1","DRV-1512-1-2","DRV-1512-2-1"]

remote_consang = ["0574-MTI","CBA-CBH-1486-3-1","MTI-1492-3-1","MTI-231-3-1","MIC-1511-3-1"]

assert len(consang) + len(non_consang) + len(remote_consang) == 109

no_wes = ["MTI-231-3-1","MIC-1300-3-10","EPI-1317-3-1","HSP-1334-4-6","CBA-1377-3-2","MIC-1385-3-1","CCH-1405-3-4","CBA-CCH-1410-3-1","DRV-1512-1-1","EPI-1317-3-3","0492-MTI"]

all_families = consang + non_consang + remote_consang
all_wes_families = [f for f in all_families if f not in no_wes]
#print len(all_wes_families)
#assert len(all_wes_families) == 97

def whosMissing() :
    plateI = "/projects/gleeson-lab/BroadData/PlateI/raw_data/Ciliopathies_Whole_Exome_Gleeson_WExBatch2.cleaned.annotated.handfiltered.txt"
    plateII = "/projects/gleeson-lab/BroadData/PlateII/raw_data/Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf"
    fin = open( plateII, "rb" )
    nline = 0

    #skip all the header files
    while True :
        nline += 1
        line = fin.readline()
        lineIsHeader = line.find( "CHROM" ) == -1
        if not lineIsHeader : break

    #process the header line w/ header colums
    headers = line.strip().split("\t")
    patients = headers[9:]
    print "%d patients" % len(patients)

    missing = [fam for fam in all_wes_families if fam not in patients]
    print "In wes_families, not in broad patients:", len(missing), missing

    missing = [fam for fam in patients if fam not in all_wes_families]
    print "In broad patients, not in wes_families:", len(missing), missing

    fin.close()

if __name__ == '__main__' :
    whosMissing()
