from collimator import Collimator, Source
import variant
import globes
from seattle import SeattleAnnotator
import db
from math import log
import broad
from plates import Pilot, PlateI, PlateII, PlateIII, CIDR
##########  CONFIGURE ########################################
dry_run = False 
switch = 'indel'
plate = CIDR()
############ /CONFIGURE  ####################################3
#### HIStory ##########

#pilot plate: snps, indels

#PlateI snps, indels

#PlateII indels
#PlateII snps

#CIDR Indels
#CIDR snps

#plateIII indels, waiting on seattle for snps
#### /HISTORY #########

conn2 = db.Conn("gleeson-closet",dry_run=dry_run)

variant_cols = conn2.getColumns("Variants")
call_cols = conn2.getColumns("Calls")
call_cols_tograb = call_cols[3:]
iso_cols = conn2.getColumns("Isoforms")
iso_cols_tograb = iso_cols[2:]
plate_id = plate.getPlateID()
#plate_contrib = 2^plate_id

#make patient name->id lookup
results = conn2.query( '''select id,name from Patients''' )
patients_dbix = {}
for (ID,name) in results :
    patients_dbix[name] = int(ID)

#for plateI snps can fast-forward 380000
ff = 0
varSource = variant.VariantList( plate.varFile(switch), fast_forward = ff )
ff = 0
seattleSource = SeattleAnnotator( plate.seattleFile(switch), switch, fast_forward = ff )
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
    #pos = variant.getPosition()[1]
    values[0] = variant_dbix

    if switch == 'indel': tipe = 2
    else : tipe = 1
    typeix = variant_cols.index("type")
    values[typeix] = tipe

    #print variant_cols
    #print values
    conn.insert( 'Variants', values, variant_cols )

    #print len(variant.isoforms), len(variant.base_calls)

    #print variant.isoforms
    for iso in variant.isoforms :
        #iso.fields['gene_id'] = geneIDFromAccession( conn, iso.fields["accession"] )
        values = [variant_dbix] + iso.getFields( iso_cols_tograb )
        conn.insert( "Isoforms", values, iso_cols[1:] )

    #now do the calls
    for call in variant.base_calls :
        pat_dbix = lookupPatientID( call )
        values = [variant_dbix, pat_dbix, plate_id] + \
                 call.getFields( call_cols_tograb )
        conn.insert( 'Calls', values, call_cols) #, skip_dupes=True )

##RETHINKING THE WHOLE UNIQUE SUM THING
##IF WE WANT ALL THE VARS CORRESPONDING TO Plate1 out,
##it will be slow to test all possible sums that have 2 bit set
##rather just select distinct(plate) from Calls where var_id = X
#multi_plate is from the 'plate' column of either Variants, Patients 
#Let P be the set of plates which have this Var/Pat in their data
#Then multi_plate_sum is SUM{x E P}(2^x)
#Return true iff plate_id contributed to this sum
def presentIn( multi_plate_sum, plateID ) :
    return (multi_plate_sum >> plateID) & 0x1

def addToCalls(conn, vid, variant) :
    #merge qual,filter,AF (or do AF separately across global pop?)
    #nah just leave and use the first one, can worry about merging later
    #if these columns turn out to be useful
    #src = conn.queryScalar("select plate from Variants where id=%d" % vid, int)
    #if not presentIn( src, plate_id ) :
        ##update the plate column if this variant is from a novel plate
        #update = "update Variants set plate = plate+%d where id=%d" \
                  #% (plate_contrib, vid)
        #conn.cur.execute( update )

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
            conn.insert( "Patients", [patient_dbix,patient], ["id","name"] )
            patients_dbix[patient] = patient_dbix
        #else :
            #src = conn.queryScalar("select plate from Patients where id=%d" % pid, int)
            #if not presentIn( src, plate_id ) :
                ##update the plate column if this variant is from a novel plate
                #update = "update Patients set plate = plate+%d where id=%d" \
                          #% ( plate_contrib, pid)
                #conn.cur.execute( update )
        
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

    c = Collimator( [varSource,seattleSource,dbSource], comparator, targetCreator, absentHandler = absentHandler )
    for i,v in enumerate(c) :
        if i % 5000 == 0 : print i
        pass #waiting for dbSource to be absent, this will trigger absentHandler

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


if __name__ == '__main__' :

    conn = db.Conn("gleeson-closet",dry_run=dry_run)
    #insertMissingVariants(conn)

    populatePatients( conn, varSource.patients )
    c = Collimator( sources, comparator, targetCreator )
    for i,v in enumerate(c) :
        if i % 1000 == 0 : print "%d variants processed" % i
        variantToDatabase( conn, v )

    print "inserted", inserted_count
    print "updated", updated_count

    ##Check that all variants showing up in the DB
    ##They are, but some calls are not being added
    ##Consider: two plates share the same patient and get the same variant
    ##attempt insert to calls, but since pk is var_id, pat_id, will be skipped
    ##therefore have extended pk to include plate
    #varSource = variant.VariantList( plate.varFile(switch), fast_forward = ff )
    #count = 0
    #for v in varSource.iterator :
        #if count % 1000 ==0 :print count
        #(chrom,pos,ref,mut) = varSource.eqkey( v )
        #chrom = globes.chromNum(chrom)
        #vid = conn.queryScalar( "select id from Variants where chrom='%s' and pos='%s' and ref='%s' and mut='%s'" % (chrom,pos,ref,mut), int )
        #if not vid :
            #print vid, count, chrom, pos, ref, mut
            #assert False
        #count += 1
    #print "count",count

