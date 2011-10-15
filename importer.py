from collimator import Collimator
import variant
import globes
from seattle import SeattleAnnotator
import db
from math import log
import broad
from plates import PlateI, PlateII, PlateIII, CIDR
#### HIStory ##########

#PlateII indels

#CIDR Indels
#### /HISTORY #########

##########  CONFIGURE ########################################
dry_run = False
switch = 'snp'
plate = CIDR()

############# /CONFIGURE  ####################################3

conn2 = db.Conn("localhost",dry_run=dry_run)

variant_cols = conn2.getColumns("Variants")
call_cols = conn2.getColumns("Calls")
call_cols_tograb = call_cols[3:]
iso_cols = conn2.getColumns("Isoforms")
iso_cols_tograb = iso_cols[2:]
plate_id = plate.getPlateID()

#make patient name->id lookup
results = conn2.query( '''select id,name from Patients''' )
patients_dbix = {}
for (ID,name) in results :
    patients_dbix[name] = int(ID)

varSource = variant.VariantList( plate.varFile(switch) )
seattleSource = SeattleAnnotator( plate.seattleFile(switch), switch )
sources = [varSource, seattleSource]

inserted_count = 0
updated_count = 0

def comparator(a,b) :
    return globes.compareVariants( a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3] )

def targetCreator() :
    return variant.Variant([],[])

def getVariantID( var_tuple ) :
    query = "select id from Variants as v where v.chrom = %s and v.pos = %s and v.ref = '%s' and v.mut = '%s'" \
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

def insertVariant(conn, variant) :
    variant_dbix = conn.getNextID("Variants")
    values = variant.getFields( variant_cols )
    pos = variant.getPosition()[1]
    values[0] = variant_dbix
    values[11] = plate_id
    if switch == 'indel': tipe = 2
    else : tipe = 1
    values[12] = tipe
    #print variant_cols
    #print values
    conn.insert( 'Variants', values, variant_cols )

    #print len(variant.isoforms), len(variant.base_calls)

    #print variant.isoforms
    for iso in variant.isoforms :
        iso.fields['gene_id'] = geneIDFromAccession( conn, iso.fields["accession"] )
        values = [variant_dbix] + iso.getFields( iso_cols_tograb )
        conn.insert( "Isoforms", values, iso_cols[1:] )

    #now do the calls
    for call in variant.base_calls :
        pat_dbix = lookupPatientID( call )
        values = [variant_dbix, pat_dbix, plate_id] + \
                 call.getFields( call_cols_tograb )
        conn.insert( 'Calls', values, call_cols) #, skip_dupes=True )

def addToCalls(conn, vid, variant) :
    place = int(log(plate_id,2))
    #merge qual,filter,AF (or do AF separately across global pop?)
    #nah just leave and use the first one, can worry about merging later
    #if these columns turn out to be useful
    src = conn.queryScalar("select plate from Variants where id=%d" % vid, int)
    if not src >> place & 0x00000001 :
        #update the plate column if this variant is from a novel plate
        update = "update Variants set plate = plate+%d where id=%d" % (plate_id, vid)
        conn.cur.execute( update )

    for call in variant.base_calls :
        pat_dbix = lookupPatientID( call )
        values = [vid, pat_dbix, plate_id] + \
                 call.getFields( call_cols_tograb )
        conn.insert( 'Calls', values, call_cols, skip_dupes=True )

@db.catch
def populatePatients( conn, patient_names ) :
    patient_dbix = conn.getNextID("Patients");
    for patient in patient_names :
        query = "select id from Patients where name='%s'" % patient
        pid = conn.queryScalar(query,int)
        if not pid :
            patient = broad.sanitizePatientName( patient )
            conn.insert( "Patients", [patient_dbix,patient,plate_id], ["id","name","plate"] )
            patients_dbix[patient] = patient_dbix
        else :
            place = int(log(plate_id,2))
            src = conn.queryScalar("select plate from Patients where id=%d" % pid, int)
            if not src >> place & 0x00000001 :
                #update the plate column if this variant is from a novel plate
                update = "update Patients set plate = plate+%d where id=%d" % (plate_id, pid)
                conn.cur.execute( update )
        
        patient_dbix += 1

def lookupPatientID( call ) :
    patient_name = varSource.patients[call.pat_ix]
    return patients_dbix[patient_name]

@db.catch
def makeEmptyGene(conn,col,value) :
    nid = conn.getNextID("Genes")
    conn.insert( "Genes", [nid,value], ["id",col] )
    return nid

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


if __name__ == '__main__' :
    conn = db.Conn("localhost",dry_run=dry_run)
    #seattleToDatabase( conn, seattleSource )    
    
    populatePatients( conn, varSource.patients )
    c = Collimator( sources, comparator, targetCreator )
    for i,v in enumerate(c) :
        if i % 1000 == 0 : print "%d variants processed" % i
        variantToDatabase( conn, v )

    print "inserted", inserted_count
    print "updated", updated_count
