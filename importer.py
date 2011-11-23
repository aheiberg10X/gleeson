from collimator import Collimator, Source
from variant import Variant
import globes
from seattle import SeattleSource
import db
from math import log
import broad
from plates import Pilot, PlateI, PlateII, PlateIII, CIDR, Frazer_ali2, Frazer_aligned

##########  CONFIGURE ########################################
# dry_run = True means everything expect the execution of any database update
# queries get run.  The queries that would have been run are printed instead
# Good for ensuring no errors get thrown half way through inserting and for 
# checking that the insert queries make sense (columns are lining up, etc)
dry_run = True

#Can specify what data is to be inserted.  It is a list of (plate,switch) 
#tuples.  Modify plates.py to add a new plate object.
plates_and_switches = [\
                       #(Pilot(),'snp'),(Pilot(),'indel'), \
                       #(PlateI(),'snp') \
                       (PlateI(),'indel') \
                       #(PlateII(),'snp'),(PlateII(),'indel'), \
                       #(PlateIII(),'snp'),\
                       #(CIDR(),'snp'),(CIDR(),'indel'), \
                       #(Frazer_ali2(),'snp'),(Frazer_ali2(),'indel'), \
                       #(Frazer_aligned(),'snp'),(Frazer_aligned(),'indel') \
                      ]

#Run with python importer.py

############ /CONFIGURE  ####################################3

############# GLOBALS ########################################

conn = db.Conn("localhost",dry_run=dry_run)

# X_cols are the columns of the database table we want to fill
# X_cols_tograb is the data we can fetch directly from the Variant object
# the other columns are things like primary and foreign keys that link 
# the rows together
variant_cols = conn.getColumns("Variants")
call_cols = conn.getColumns("Calls")
call_cols_tograb = call_cols[3:]
iso_cols = conn.getColumns("Isoforms")
iso_cols_tograb = iso_cols[2:]

#make patient name->id lookup
results = conn.query( "select id,name from Patients" )
patients_dbix = {}
for (ID,name) in results :
    patients_dbix[name] = int(ID)

############## /Globals  ###########################################


def insertPlate( conn, plate, switch ) :
    inserted_count = 0
    updated_count = 0

    #We will keep track which plate each call comes from.
    plate_id = plate.getPlateID()

    # Construct the Sources for a Collimator to work on
    vcfSource = broad.VCFSource( plate.varFile(switch) )
    seattleSource = SeattleSource( plate.seattleFile(switch), switch )
    sources = [vcfSource, seattleSource]

    populatePatients( conn, vcfSource.patients )
    c = Collimator( sources, comparator, targetCreator )
    for i,var in enumerate(c) :
        if i % 5000 == 0 : print "%d variants processed" % i

        values = var.getPosition()
        vid = getVariantID( values )
        if not vid :
            inserted_count += 1
            insertVariant( conn, var, plate_id )
        else :
            updated_count += 1
            addToCalls( conn, vid, var, plate_id )

    print "Plate: %d" % plate_id
    print "    inserted", inserted_count
    print "    updated", updated_count

#method needed to construct a Collimator
def comparator(a,b) :
    return globes.compareVariants( a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3] )

#method needed to construct a Collimator
def targetCreator() :
    return variant.Variant([],[])

#Would very much like to speed up this process, dont see how
def getVariantID( var_tuple ) :
    query = '''select id
               from Variants as v
               where v.chrom = %s
                 and v.pos = %s
                 and v.ref = '%s'
                 and v.mut = '%s' ''' \
            % (var_tuple[0],var_tuple[1],var_tuple[2],var_tuple[3])

    return conn.queryScalar( query, int )

def variantToDatabase( conn, variant ) :
    global inserted_count, updated_count
    keys = ["chrom","pos","ref","mut"]
    values = [variant.fields[k] for k in keys]
    vid = getVariantID( values )
    if not vid :
        inserted_count += 1
        insertVariant( conn, variant )
    else :
        updated_count += 1
        addToCalls( conn, vid, variant )

#We have a brand new variant, put it into the db
def insertVariant(conn, variant, plate_id) :
    variant_dbix = conn.getNextID("Variants")
    values = variant.getFields( variant_cols )
    #pos = variant.getPosition()[1]
    values[0] = variant_dbix

    if switch == 'indel': tipe = 2
    else : tipe = 1
    typeix = variant_cols.index("type")
    values[typeix] = tipe

    #insert the variants
    conn.insert( 'Variants', values, variant_cols )

    #insert the isoforms
    for iso in variant.isoforms :
        #iso.fields['gene_id'] = geneIDFromAccession( conn, iso.fields["accession"] )
        values = [variant_dbix] + iso.getFields( iso_cols_tograb )
        conn.insert( "Isoforms", values, iso_cols[1:] )

    #insert the calls
    for call in variant.base_calls :
        pat_dbix = patients_dbix[call.pat_name]
        values = [variant_dbix, pat_dbix, plate_id] + \
                 call.getFields( call_cols_tograb )
        conn.insert( 'Calls', values, call_cols) #, skip_dupes=True )

# This variant (identified by vid) has already been added to the system
# By definition so have the isoforms
# The calls, however, are novel
def addToCalls(conn, vid, variant, plate_id) :
    #Do we want to integrate quality score here?  Doing a weighted average
    #by the number of calls?
    for call in variant.base_calls :
        #pat_dbix = lookupPatientID( call )
        pat_dbix = patients_dbix[call.pat_name]
        values = [vid, pat_dbix, plate_id] + \
                 call.getFields( call_cols_tograb )
        conn.insert( 'Calls', values, call_cols, skip_dupes=True )

#Add the patient names to the Patients table
@db.catch
def populatePatients( conn, patient_names ) :
    patient_dbix = conn.getNextID("Patients");
    for patient in patient_names :
        query = "select id from Patients where name='%s'" % patient
        pid = conn.queryScalar(query,int)
        if not pid :
            patient = broad.sanitizePatientName( patient )
            conn.insert( "Patients", [patient_dbix,patient], ["id","name"] )
            patients_dbix[patient] = patient_dbix

        patient_dbix += 1

#Given a BaseCall object (see variant.py), get the actual patient name
def lookupPatientID( call ) :
    patient_name = vcfSource.patients[call.pat_ix]
    return patients_dbix[patient_name]

if __name__ == '__main__' :
    for plate,switch in plates_and_switches :
        insertPlate( conn, plate, switch )


###############################################################################
###############################################################################
#unfinished.  Idea was to remove Calls that were spurious,
#then removed the Variants, Patients, and Isoforms that got orphaned
def removeCalls() :
    plate_id = globes.plates("PlateIII")

    #delete all Calls that are old
    query = '''
delete c
from Calls as c inner join Variants as v on c.var_id = v.id
where c.plate =  and type = '''

    #delete any stranded variants
    query = '''
    delete v
    from Variants as v left join Calls as c on v.id = c.var_id
    where c.var_id is null
    '''

    #delete any stranded patients
    query = '''delete p 
    FROM Patients AS p
    LEFT JOIN (

        SELECT DISTINCT pat_id
        FROM Calls
    ) AS t ON t.pat_id = p.id
    WHERE pat_id IS NULL
    '''


#############################################################################
            #OLd stuff
#############################################################################

#Hold over from when we were maintaining a Genes table.  Not used currently
@db.catch
def makeEmptyGene(conn,col,value) :
    nid = conn.getNextID("Genes")
    conn.insert( "Genes", [nid,value], ["id",col] )
    return nid

#Hold over from when we were maintaining a Genes table.  Not used currently
def geneIDFromAccession( conn, accession ) :
    query = "select id from Genes where refseq = '%s'" % accession
    gene_ids = conn.query( query )
    if not gene_ids :
        gid = makeEmptyGene( conn, 'refseq', accession )
    elif len(gene_ids) == 1 :
        gid = int(gene_ids[0][0])
    else :
        #'what to do when an accession: %s matches multiple Genes %s?' % (accession, str(gene_ids))
        #take the first one
        gid = int(gene_ids[0][0])

    return gid

#quick hack to add in some missing variants
def insertMissingVariants(conn) :
    query = '''
   SELECT id,chrom,pos,ref,mut
   FROM (

       SELECT var_id
       FROM Calls AS c
       INNER JOIN Variants AS v ON c.var_id = v.id
       WHERE c.plate = 3
       AND TYPE =1
       GROUP BY var_id
   ) AS t inner join Variants as v on t.var_id = v.id
   ORDER BY chrom,pos,ref,mut
   '''

    rows = conn.query( query )
    print len(rows)
    it = iter(rows)
    def eqkey( it ) :
        return [int(it[1]), int(it[2]), it[3], it[4]]

    def integrator( target, it ) :
        return target

    def absentHandler( ae ) :
        print str(ae)
        if ae.ix == 2 :
            variantToDatabase( conn, ae.missed_target )
        else :
            assert False

    dbSource = Source( it, eqkey, integrator, allow_absent = False )

    c = Collimator( [vcfSource,seattleSource,dbSource], comparator, targetCreator, absentHandler = absentHandler )
    for i,v in enumerate(c) :
        if i % 5000 == 0 : print i
        pass #waiting for dbSource to be absent, this will trigger absentHandler

#quick hack to fix fill in the ref_aa/mut_aa columns
def seattleToDatabase( conn, seattleSource ) :
    count = 0
    for out_splt in seattleSource.iterate(fast_forward=0) :
        if count % 1000 == 0 : print "%d updated" % count
        count += 1
        vid = getVariantID( seattleSource.getPosition( out_splt ) )
        if vid :
            aas = out_splt[11].split(',')
            if len(aas) == 2 :
                conn.update('Variants',aas,['ref_aa','mut_aa'], vid)
        else :
            print "Seattle Row: %s has no matching variant in the Variants table" % str(out_splt)
            assert False

