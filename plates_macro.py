def addPlateMacro( folder_name, \
                   snpfile, \
                   indelfile, \
                   seattle_snp_file, \
                   seattle_indel_file ) :

    #add the new folder_name : id pair
    fin = open( "plate_ids.txt", 'r+a' )
    lines = fin.readlines()
    last_id = int(lines[-1].split(":")[1])
    previous_names = []

    for line in lines :
        previous_names.append( line.split(":")[0] )

    if folder_name in previous_names :
        raise Exception("'%s' is already the name of a plate" % folder_name)

    fin.write( "%s:%d\n" % (folder_name, last_id+1) )
    fin.close()

    #now edit this very same file to define a new class
    fh = open( "plates.py", "a" )
    fh.write( \
'''class %s(Plate) :
    def __init__(self) :
        self.folder_name = "%s"
        Plate.__init__(self)
        self.snpfile = "%s"
        self.indelfile = "%s"
        self.seattle_snp_file = "%s"
        self.seattle_indel_file = "%s"\n\n''' % \
             (folder_name, folder_name, \
              snpfile, indelfile, \
              seattle_snp_file, seattle_indel_file ) )

    fh.close()


if __name__ == '__main__' :
    addPlateMacro( "frazer3", "from_frazer_snps.vcf", "from_frazer_indels.vcf", "SeattleSeqAnnotation134.frazer3_snps.vcf.tar.gz.230368586074.txt", "SeattleSeqAnnotation134.from_frazer_indels.vcf.230369948975.txt" )

