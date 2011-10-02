import globes
import broad
import pickle

class Gene :
    def __init__(self) :
        self.exonList = []

class Exon :
    def __init__( self, start, end ) :
        self.start = start
        self.end = end

class Indel :
    def __init__(self) :
        pass   # see inputINDELS

    def isCovered( self, thresh = 8 ) :
        return self.dp >= thresh

#right now return True iff only one patient has this indel.
#NOTE: if it is a hanging indel (no patient has this variation in the original
#      vcf, we still return false. However, these are filtered out in 
#      separateSNPSandINDELS()
#TODO: what about if two patients with same disease have the same indel
#      don't want to filter then, right? Extremely rare...
def indelUniqueToDisease( patients, calls, dgroups={} ) :
    #NEED SOME DATA STRUCTURE with patients grouped by disease
    #{patient1 : group1, patient2: group2, patient3 : group1,. ..}
    #but for right now assume every patient is independent of every other
    if not dgroups :
        dgroups = {}
        for patient in patients :
            dgroups[patient] = patient

    groups_with = []

    for patient,call in zip(patients,calls) :
        call_splt = broad.splitCall( call )
        if broad.isMutated( broad.convertGT(call_splt) ) :
            patients_group = dgroups[patient]
            if patients_group not in groups_with :
                groups_with.append( patients_group )

    count = len(groups_with)
    if count == 0 :
        return (False,"No patient has this variant")
    elif count == 1 :
        return (True,"boom")
    else :
        return (False,
                "Indel shows up across unrelated patients: %s" 
                 % (groups_with) )

class IndelList :
    def __init__(self,patients) :
        #list of Indels()
        self.indels = []
        self.nindels = 0
        #maintain the indexes of the Indels that belong to each patient
        self.patient_index = {}
        for pat in patients :
            self.patient_index[pat] = []

    def add( self, indel, patient ) :
        self.indels.append( indel )
        self.patient_index[patient].append( self.nindels )
        self.nindels += 1

    def patientIndels( self, patient ) :
        return [ self.indels[i] for i in self.patient_index[patient] ]


#remember this also filters for unique. 
def inputINDELS( filename ) :
    #k: patient, v: [indels]
    fout = open( "%s/raw_data/input_indel_err.txt" % (globes.DATA_DIR), 'wb' )
    fin = open( filename, 'rb' )

    patients = broad.getPatients( fin )
    il = IndelList(patients)

    indexOf = broad.COLUMN_MAP
    not_unique_count = 0

    for indel_line in fin.readlines() :
        splt = indel_line.strip().split('\t')
        calls = splt[ indexOf["calls"]: ]

        info = splt[ indexOf["info"] ]
        dinfo = broad.makeInfoDict( info )

        ref = splt[ indexOf["ref"] ]
        mut = splt[ indexOf["mut"] ]

        (unique, message) = indelUniqueToDisease( patients, calls )
        if unique :
            for c,call in enumerate(calls) :
                call_splt = broad.splitCall(call)
                GT = broad.convertGT( call_splt )
                if broad.isMutated( GT ) :
                    (gt,AD,DP,GQ,PL) = call_splt
                    indl = Indel()
                    indl.chrom = splt[ indexOf["chrom"] ]
                    indl.start = int( splt[ indexOf["pos"] ] )
                    indl.ref = ref
                    indl.mut = mut
                    indl.ac = int( dinfo['AC'] )
                    indl.dp = DP
                    indl.is_splice, indl.is_utr = False, False
                    indl.ProtStart, indl.ProtEnd = 0,0

                    il.add( indl, patients[c] )
        else :
            not_unique_count += 1
            fout.write( "\n\nloc: %s, %s" % (splt[indexOf['pos']],message) )

    print "not unqiue: %d" % not_unique_count
    fin.close()
    fout.close()

    return il

#redundant, not that valuable anyway as 
def filterGenic( indelList, geneList ) :
    newIndelList = []
    num_genes = len(geneList)

    for patient in indelList :
        k = 0
        for indel in patient :
            start = indel.start+1

            isDeletion = len(indel.ref) > len(indel.alt)
            if isInsertion : end = indel.start + len(indel.ref)
            else :           end = indel.start+1

            if k < num_genes :
                reached_chr = (geneList[k].chrom == indel.chrom)
            else :
                reached_chr = True

            is_inserted = False

            #get to the right chrom
            while k < num_genes and \
                  (not reached_chr or \
                   (geneList[k].chrom == indel.chrom and \
                    geneList[k].end < start) \
                  ) :
                pass

            kprime = k
            cds_start, cds_end = 0,0
            splice_buffer = globes.SPLICE_BUFFER

            for gene in geneList :
                if gene.strand :
                    for exon in gene.exonList : #break when is_inserted
                        if exon.start - splice_buffer <= end and \
                           exon.end + splice_buffer >= start and \
                           exon.end >= start : #why this last one?
                            if indel_in_splice_buffer : pass
                            elif in_utr : pass
                            else : pass
                        else : pass
                else :
                    for exon in reversed( gene.exonList ) :
                        if exon.start - splice_buffer <= end and \
                           exon.end + splice_buffer >= start and \
                           exon.start <= start :
                            if True: pass
                            elif True: pass
                            else : pass
                        else : pass
                kprime += 1

def inputGeneList( gene_list_filename ) :
    geneList = []

    fin = open( gene_list_filename )
    indexOf = {"#bin" :       0,\
               "name" :       1,\
               "chrom" :      2,\
               "strand" :     3,\
               "txStart" :    4,\
               "txEnd" :      5,\
               "cdsStart" :   6,\
               "cdsEnd" :     7,\
               "exonCount" :  8,\
               "exonStarts" : 9,\
               "exonEnds" :   10,\
               "name2" :      12}

    header = fin.readline()
    for line in fin :
        g = Gene()
        splt = line.strip().split('\t')

        #name attribute is fancy, handle it explicitly
        g.name = "%s(%s)" % (splt[indexOf["name2"]], splt[indexOf["name"]])

        #otherwise set attributes of g en masse
        attributes = ["chrom","strand","txStart","txEnd",\
                      "cdsStart","cdsEnd"]
        for attr in attributes :
            setattr( g, attr, splt[ indexOf[attr] ] )

        #turn strand attr into a boolean
        g.strand = (g.strand == '+')

        #for some reason there is a trailing ',' on starts and ends
        exonCount = int(splt[ indexOf["exonCount"] ])
        exonStarts = splt[ indexOf["exonStarts"] ].strip(',').split(',')
        exonEnds = splt[ indexOf["exonEnds"] ].strip(',').split(',')
        assert len(exonStarts) == exonCount == len(exonEnds)

        exonList = [ Exon( int(se[0]),int(se[1]) ) \
                     for se in zip(exonStarts,exonEnds)]
        g.exonList.append( exonList )
        geneList.append( g )

    fin.close()
    return geneList


if __name__ == "__main__" :
    indelList = inputINDELS( globes.INDEL_FILE )
    print len( indelList.indels )
    print len( indelList.patientIndels( "SPOAN-1513" ) )


    #indelPipeline :
        #inputIndels : There are LOTS of indel variations without any representation in the patients.  Broad has no explanation at the moment.  For the moment they will be ignored.
        #filterIndels : apparently determining homozygosity from frequency is not reliable, so deferring this filter till later
        #filterUnique : this is being done in input, helped by new format
        #inputGeneList : done, not tested
        #filterGenic : not started
        #inputDomains : not started
        #outputIndels :

        #also, at some point this needs to go through SIFT, where and when?



#THIS IS OLD, MOVED TO BROAD inputIndel()
#old because the new separateSNPsAndINDELS dedupes by varaint already
#def inputIndels( filename ) :
    ##k: patient, v: [indels]
    #indelList = {}
#
    #fin = open( filename, 'rb' )
    #patients = broad.getPatients( fin )
    #for pat in patients : indelList[pat] = []
    #indexOf = broad.BROAD_COLUMN_MAP
#
    #hanging_indels = []
    #not_unique = []
    #tot_indels = 0
#
    #for ln, indel_line in enumerate(fin.readlines()) :
        #splt = indel_line.strip().split('\t')
        #calls = splt[ indexOf["calls"]: ]
#
        #ref = splt[ indexOf["ref"] ]
        #muts = splt[ indexOf["mut"] ].split(',')
#
        #info = splt[ indexOf["info"] ].split(";")
        #dinfo = broad.makeInfoDict( splt[indexOf["info"]] )
        #allele_counts = dinfo["AC"].split(',')
#
        #assert len(allele_counts) == len(muts)
#
        #pos = int(splt[indexOf["pos"]])
        #for variation in range(len(muts)) :
            #indel_has_unique_patient = False
            #unique_patient = ""
            #indel = Indel()
            #for c,call in enumerate(calls) :
                #call_splt = broad.splitCall(call,variation+1)
                ##if pos == 907630 : print call_splt
                #(GT,AD,DP,GQ,PL) = call_splt
                #if broad.isMutated( GT ) :
                    ##if pos == 907630 : print "is mutated!"
                    #if indel_has_unique_patient :
                        #indel_has_unique_patient = False
                        #not_unique.append( pos )
                        #break
                    #else :
                        ##if pos == 907630 : print "setting unique_patient=42"
                        #indel_has_unique_patient = True
                        #unique_patient = patients[c]
                        #indel = Indel()
                        #indel.chrom = splt[ indexOf["chrom"] ]
                        #indel.start = int( splt[ indexOf["pos"] ] )
                        #indel.ref = ref
                        #indel.mut = muts[variation]
                        #indel.ac = allele_counts[variation]
                        #indel.dp = DP
                        #indel.is_splice, indel.is_utr = False, False
                        #indel.ProtStart, indel.ProtEnd = 0,0
#
            #if unique_patient == "" :
                #hanging_indels.append( (pos,variation) ) 
#
            #if indel_has_unique_patient :
                #indelList[unique_patient].append( indel )
#
            #tot_indels += 1
#
    ##for pat in indelList :
        ##print "-----"
        ##print pat
        ##if len(indelList[pat]) : print indelList[pat][0].start
#
    #fin.close()
#
    ##fpick = open("pickled.pck",'wb')
    ##fpick.write( pickle.dumps( indelList ) )
    ##fpick.close()
#
