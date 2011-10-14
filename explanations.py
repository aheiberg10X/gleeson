template = '''Report variants in <interval> with >= <include_atleast> <include_genotype> belonging to someone in <group>, excluding those genes with <exclude_genotype> calls belonging to someone in <exclude_group>.

Calls should have depth >= <call_depth> and qual >= <call_qual>

Only count variants where <variant_restrictions>'''

detail = '''Template:

Report variants in <interval> with >= <include_atleast> <include_genotype> belonging to someone in <group>, excluding those genes with <exclude_genotype> calls belonging to someone in <exclude_group>.

Calls should have depth >= <call_depth> and qual >= <call_qual>

Only count variants where <variant_restrictions>

=========================================================================

Syntax :

    All values must be quoted: "value".  Values cannot be empty ("") unless otherwise noted.

    <interval> : can be one of 'Gene', 'Exon', 'Variant'
    <include_atleast> : must be an integer
    <include_genotype> : must be 'Hom','Het','Both'
    <group> : must be a comma separated list of patient_ids, i.e '200,202,203'. See below on how to lookup these ID's
    <exclude_genotype> : same requirements as <include_genotype>
    <exclude_group> : same requirements as <group>
    <call_depth> : must be an integer.  Can be empty
    <call_qual> : must be an integer. Can be empty
    <variant_restrictions> : edit this at your own risk.  A sensible default is to use: 
    " v.dbSNP = '.' and (i.ss_functionGVS = 'missense' or i.ss_functionGVS = 'nonsense' or i.ss_functionGVS = 'splice-3' or ss_functionGVS = 'splice-5' or (ss_functionGVS is null and effect = 'NON_SYNONYMOUS_CODING'))"

provided in config/example.txt.

Don't forget the containing {} around everything!

=========================================================================

Semantics :
 
    Say we want all genes that have >=3 homozygous variants belonging to any one of patients 200,202,203.  But we want to exclude those genes with homozygous variants belonging to the unaffecteds, and we know 201,204 are definitely unaffected and have no phenotypic overlap with the affecteds.  

So we let "interval" : "Gene"
          "include_atleast" = "3"
          "include_genotype" = "Hom"
          "group" = "200,202,203"
          "exclude_genotype" = "Hom"
          "exclude_group" = "201,204"

Now we probably don't want to consider reads with low read depth, so we can let 
          <call_depth> = 8 
To exclude from consideration all reads with depth < 8.  If we wanted to lower bound call quality we could do somthing like <call_qual> = 90.

I did not get around to paramaterizing <variant_restrictions>, if you want to change it you need to enter valid SQL. Email me @ aheiberg@ucsd.edu if you have a n idea of what you want to do but don't know how.'''
