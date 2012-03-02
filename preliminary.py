import db
import globes
import csv
import queries

conn = db.Conn("localhost")
conn2 = db.Conn("localhost")

patients = ["HCD-1577-2-1", "HME-1563-2-1", "HME-1565-2-4", "HME-1573-2-1", "HME-1574-2-2", "HME-1620-2-2"]


patient = "HME-1620-2-2"
fname = "%s_SNV_result.txt.filtered" % patient
fin = open("../preliminary/%s" % fname)
fout = open("../preliminary/%s.functional.tsv" % patient,'w')
csvout = csv.writer( fout,\
                     delimiter='\t', \
                     quoting=csv.QUOTE_MINIMAL )
csvout.writerow( queries.column_headers )

hits = 0
already_called = 0
total = 0

for line in fin.readlines() :
    splt = line.split()
    chrom, pos = splt[0], splt[1]
    #print chrom, pos
    chrom = globes.chromNum( chrom )
    if not chrom : continue

    query = '''
select %s,%s,%s
from Variants as v inner join Isoforms as i on i.var_id = v.id
                   inner join Genes as g on g.id = i.gene_id
where chrom = %s and pos = %s and AF < .1 and (%s)''' % \
    (queries.vcols_string, queries.icols_string, queries.gcols_string, \
     chrom, pos, queries.gvs)


    rows = conn.query(query)
    num_rows = len(rows)
    if num_rows > 0 :
        hits += 1

    total += 1

    for row in rows :
        outrow = queries.formatQueryRow(row)
        (noinfs, hets, homs) = queries.getPatients( conn2, row[0] )
        hets_string = ", ".join([t[1] for t in hets])
        homs_string = ", ".join([t[1] for t in homs])
        if patient in hets_string or patient in homs_string :
            already_called += 1
        outrow.extend( ["-",len(homs),homs_string,len(hets),hets_string] )
        csvout.writerow( outrow )
        # hom and het shares

print patient
print "%d/%d hits/total" % (hits,total)
print "%d/%d already_called/hits" % (already_called,hits)
fin.close()
fout.close()

