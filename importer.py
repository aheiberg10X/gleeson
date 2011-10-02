import db
import snp
from annotator import walker, Walker2
import snpeff
import globes
import sift
import seattle
from math import log
import time

@db.catch
def makeEmptyGene(conn,col,value) :
    nid = conn.getNextID("Genes")
    conn.insert( "Genes", [nid,value], ["id",col] )
    return nid

def reset() :
    conn = db.Conn("gleeson-closet")
    conn.wipe("Variants")
    conn.wipe("Isoforms")
    conn.wipe("Calls")

@db.catch
def populatePatients( patient_names ) :
    conn = db.Conn("gleeson-closet")
    patient_dbix = conn.getNextID("Patients");
    for patient in patient_names :
        query = "select id from Patients where name='%s'" % patient
        pid = conn.queryScalar(query,int)
        if not pid :
            conn.insert( "Patients", [patient_dbix,patient], ["id","name"] )
        else :
            pass #update pid
        patient_dbix += 1

def Walker2Test() :
    t1 = time.time()
    dbconn = db.Conn()
    varList = snp.SNPList( globes.SNP_FILE )

    effAnn = snpeff.SNPEffAnnotator()
    dargs = {'file':"%s/snpeff/output/Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.csv" % globes.INT_DIR, \
             'dbconn' : dbconn}
    effAnn.register(dargs)

    siftAnn = sift.SIFTAnnotator('snp')
    dargs = {"prediction" : 9534}
    siftAnn.register( dargs )

    ssAnn = seattle.SeattleAnnotator()
    dargs = {"file" : "/home/andrew/gleeson/PlateII/intermediate_data/seattle/output/SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf.gz.213169538436.txt", \
             'dbconn' : dbconn}
    ssAnn.register( dargs )

    def callback( variant ) :
        #print variant.fields 
        pass

    w = Walker2( [varList,effAnn,siftAnn,ssAnn] ) 
    for t in w.iterate() :
        callback(t)

    t2 = time.time()
    print "done, took %f s" % (t2-t1)

def updateWithAnnotations() :
    conn = db.Conn()
    cols = ["id","chrom","pos","ref","mut","gene_id"]
    iSQL = conn.iterate("Variants",cols=cols)

    conn2 = db.Conn()
    effAnn = snpeff.SNPEffAnnotator()
    dargs = {'file' : "/home/andrew/gleeson/PlateII/intermediate_data/snpeff/output/Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.csv", \
             'dbconn' : conn2}
    effAnn.register(dargs)
    effAnn.annotate( iSQL )

    #siftAnn = sift.SIFTAnnotator('snp')
    #dargs = {"prediction" : 9534, "dbconn" : conn2}
    #siftAnn.register( dargs )
    #siftAnn.annotate( iSQL )

    #ssAnn = seattle.SeattleAnnotator()
    #dargs = {"file" : "/home/andrew/gleeson/PlateII/intermediate_data/seattle/output/SeattleSeqAnnotation131.Ciliopathies_Whole_Exome_Gleeson_Indels_20110526_SNPS.vcf.gz.213169538436.txt",\
             #"dbconn" : conn2}
    #ssAnn.register( dargs )
    #ssAnn.annotate( iSQL )

@db.catch
def insertVariantList( ) :
    dry_run = False 
    dbconn = db.Conn("gleeson-closet")
    varList = snp.SNPList( globes.INDEL_FILE )
    populatePatients( varList.patients )
    switch = 'indel'

    effAnn = snpeff.SNPEffAnnotator()
    dargs = {'file':"%s/snpeff/output/indels.txt" % globes.INT_DIR, \
             'dbconn' : dbconn, \
             'switch' : switch}
    effAnn.register(dargs)

    #siftAnn = sift.SIFTAnnotator( 'snp' )
    #dargs = {"prediction" : 9534}
    #siftAnn.register(dargs)

    ssAnn = seattle.SeattleAnnotator()
    dargs = {'file' : "%s/seattle/output/SeattleSeqAnnotation131.pilot_full_indels.vcf.216816932058.txt" % globes.INT_DIR, \
             'dbconn' : dbconn, \
             'switch' : switch}
    ssAnn.register(dargs)

    w = Walker2( [varList,effAnn,ssAnn] )
    iWalk = w.iterate()

    fout = open("insert_varlist_output.txt",'w')
    global update_count, insert_count, place
    insert_count = 0
    update_count = 0
    place = int(log(globes.SOURCE,2))

    conn1 = db.Conn("gleeson-closet")
    conn2 = db.Conn("gleeson-closet")
    cols = ["id","chrom","pos","ref","mut","source"]
    iSQL = conn1.iterate("Variants", cols=cols, order_by=cols[1:]) #make sure it is sorted
    #iVar = varList.iterate()

    #make patient name->id lookup
    results = conn2.query( '''select id,name from Patients''' )
    patients_dbix = {}
    for (ID,name) in results :
        patients_dbix[name.lower()] = int(ID)

    #columns of the tables as well as the ones to get out of variant.fields
    variant_cols = conn2.getColumns("Variants")
    #print variant_cols
    #variant_cols_tograb = variant_cols[1:-1]
    call_cols = conn2.getColumns("Calls")
    call_cols_tograb = call_cols[2:]
    iso_cols = conn2.getColumns("Isoforms")
    iso_cols_tograb = iso_cols[2:]

    def comp(sqlrow, variant) :
        #if at the end of iSQL, return 1 so gtfunc gets called on new 'variant'
        if sqlrow == db.endOfIteration : return 1
        (vid,chrom1,pos1,ref1,mut1,source1) = sqlrow
        chrom2,pos2,ref2,mut2 = variant.getPosition()
        r =  globes.compareVariants( chrom1,pos1,ref1,mut1, \
                                       chrom2,pos2,ref2,mut2 )
        #print chrom1,pos1,ref1,mut1,source1,chrom2,pos2,ref2,mut2,globes.SOURCE
        #print r
        #print "\n"
        return r

    def insertVariant(sqlrow, variant) :
        global insert_count
        insert_count += 1
        variant_dbix = conn2.getNextID("Variants")
        values = variant.getFields( variant_cols )
        pos = variant.getPosition()[1]
        values[0] = variant_dbix
        values[11] = globes.SOURCE
        if switch == 'indel': tipe = 2
        else : tipe = 1
        values[12] = tipe
        #print variant_cols
        #print values
        if not dry_run :
            conn2.insert( 'Variants', values, variant_cols )

        #print len(variant.isoforms), len(variant.base_calls)

        #print variant.isoforms
        for iso in variant.isoforms :
            values = [variant_dbix] + iso.getFields( iso_cols_tograb )
            if not dry_run :
                conn2.insert( "Isoforms", values, iso_cols[1:] )

        #now do the calls
        for call in variant.base_calls :
            patient_name = varList.patients[call.pat_ix]
            pat_dbix = patients_dbix[patient_name]
            values = [variant_dbix, pat_dbix] + \
                     call.getFields( call_cols_tograb )
            if not dry_run :
                conn2.insert( 'Calls', values, call_cols) #, skip_dupes=True )


    def addToCalls(sqlrow, variant) :
        global update_count, place
        update_count += 1

        #merge qual,filter,AF (or do AF separately across global pop?)
        #nah just leave and use the first one, can worry about merging later
        #if these columns turn out to be useful
        (vid,chrom,pos,ref,mut,source) = sqlrow
        src = conn2.queryScalar("select source from Variants where id=%d" % vid, int)
        if not src >> place & 0x00000001 :
            #update the source column if this variant is from a novel source
            update = "update Variants set source = source+%d where id=%d" % (globes.SOURCE, vid)
            if not dry_run :
                conn2.cur.execute( update )

        for call in variant.base_calls :
            patient_name = varList.patients[call.pat_ix]
            pat_dbix = patients_dbix[patient_name]
            values = [vid, pat_dbix] + \
                     call.getFields( call_cols_tograb )
            if not dry_run :
                conn2.insert( 'Calls', values, call_cols, skip_dupes=True )


    walker( iSQL, iWalk, comp=comp, eqfunc=addToCalls, gtfunc=insertVariant )
    print "inserted: %d, updated: %d" % (insert_count,update_count)
    fout.close()

def insertVariantList_old( varList ) :
    calls_inserted = 0
    calls_avoided = 0

    conn = db.Conn()
    #create a patient name/ db ix lookup dict
    results = conn.query( '''select id,name from Patients''' )
    patients_dbix = {}
    for (ID,name) in results :
        patients_dbix[name] = int(ID)

    variant_cols = conn.getColumns("Variants")
    variant_cols_tograb = variant_cols[1:-1]

    call_cols = conn.getColumns("Calls")
    call_cols_tograb = call_cols[2:]

    iso_cols = conn.getColumns("Isoforms")
    iso_cols_tograb = iso_cols[2:]

    ##step 1
    #build up patient database

    vars_iter = varList.iterate()
    #get max variant_id, +1, this will be our starting ID
    variant_dbix = conn.getNextID("Variants")
    for snp in vars_iter :
        values = [variant_dbix] + \
                 snp.getFields( variant_cols_tograb ) + \
                 ['plateII']
        try :

            conn.insert( 'Variants', values, variant_cols )
            for iso in snp.isoforms :
                #pk is autoinc
                values = [variant_dbix] + iso.getFields( iso_cols_tograb )
                conn.insert( "Isoforms", values, iso_cols[1:] )
        except db.AlreadyExists :
            assert False#update instead of create?
            continue

        #now do the calls
        for call in snp.base_calls :
            #call_is_hom_ref = call.fields["GT"] == 1
            #if not call_is_hom_ref :
            calls_inserted += 1
            patient_name = varList.patients[call.pat_ix]
            pat_dbix = patients_dbix[patient_name]
            values = [variant_dbix, pat_dbix] + \
                     call.getFields( call_cols_tograb )
            conn.insert( 'Calls', values, call_cols )
            #else : 
                #calls_avoided += 1

        variant_dbix += 1

    #wtf figure out why this isn't valid
    #conn.commit()
    #conn.close()

    print "calls_inserted: ", calls_inserted
    print "calls_avoided: ", calls_avoided

if __name__ == "__main__" :
    #varList = snp.SNPList( frazer1 )
    #effAnn = snpeff.SNPEffAnnotator( frazer1 )
    #effAnn.run()    
    #updateWithAnnotations()
    ####reset()
    insertVariantList( )
    #Walker2Test()
    #makeEmptyGene( db.Conn() )

#ALTER TABLE  `Variants` ADD  `ss_functionGVS` VARCHAR( 30 ) NULL ,
#ADD  `ss_polyPhen` VARCHAR( 30 ) NULL ,
#ADD  `ss_granthamScore` INT UNSIGNED NULL ,
#ADD  `ss_scorePhastCons` FLOAT NULL ,
#ADD  `ss_consScoreGERP` FLOAT NULL ,
#ADD  `ss_distanceToSplice` INT UNSIGNED NULL ,
#ADD  `ss_AfricanHapMapFreq` FLOAT NULL ,
#ADD  `ss_EuropeanHapMapFreq` FLOAT NULL ,
#ADD  `ss_AsianHapMapFreq` FLOAT NULL

#CREATE TABLE  `gleeson`.`calls` (
#`id` INT UNSIGNED NOT NULL ,
#`var_id` INT UNSIGNED NOT NULL ,
#`pat_id` INT UNSIGNED NOT NULL ,
#`genotype` VARCHAR( 3 ) NOT NULL ,
#`read_depth` INT NULL ,
#`allelic_depth_ref` INT NULL ,
#`allelic_depth_mut` INT NULL ,
#`phred_wt` DOUBLE NULL ,
#`phred_het` DOUBLE NULL ,
#`phred_hom` DOUBLE NULL ,
#PRIMARY KEY (  `id` )
#) ENGINE = MYISAM ;


##############################################################################
#CREATE TABLE  `gleeson`.`Patients` (
#`id` INT NOT NULL ,
#`name` VARCHAR( 20 ) NOT NULL ,
#`family` INT NOT NULL ,
#`disease` VARCHAR( 100 ) NOT NULL ,
#PRIMARY KEY (  `id` )
#) ENGINE = MYISAM ;


#####################################################################

#CREATE TABLE  `gleeson`.`Isoforms` (
#`id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY ,
#`effect` VARCHAR( 30 ) NOT NULL ,
#`transcript_id` VARCHAR( 20 ) NOT NULL ,
#`exon_id` VARCHAR( 20 ) NOT NULL ,
#`codon_num(cds)` INT UNSIGNED NOT NULL ,
#`cds_size` INT UNSIGNED NOT NULL
#) ENGINE = MYISAM ;

#ALTER TABLE  `Isoforms` ADD  `exon_rank` INT UNSIGNED NOT NULL AFTER  `exon_id`

#ALTER TABLE  `Isoforms` ADD  `variant_id` INT UNSIGNED NOT NULL AFTER  `id`


#CREATE TABLE  `gleeson`.`Sources` (
    #`id` INT NOT NULL ,
    #`name` VARCHAR( 50 ) NOT NULL ,
    #`size` INT NOT NULL ,
    #`reference` VARCHAR( 20 ) NULL ,
    #`aligner` VARCHAR( 20 ) NULL ,
    #`caller` VARCHAR( 20 ) NULL ,
    #`notes` TEXT NULL ,
    #PRIMARY KEY (  `id` )
#) ENGINE = MYISAM ;


#CREATE TABLE  `gleeson`.`Genes2` (
#`id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY ,
#`ucsc_id` VARCHAR( 20 ) NULL ,
#`mRNA` VARCHAR( 20 ) NULL ,
#`spID` VARCHAR( 20 ) NULL ,
#`spDisplayID` VARCHAR( 20 ) NULL ,
#`geneSymbol` VARCHAR( 20 ) NULL ,
#`refseq` VARCHAR( 20 ) NULL ,
#`protAcc` VARCHAR( 20 ) NULL ,
#`description` VARCHAR( 100 ) NULL
#) ENGINE = MYISAM ;
