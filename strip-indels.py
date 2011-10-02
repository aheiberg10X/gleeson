# strip-indels
# program to remove indel lines from Broad data
# originated by Jack Cassidy - June 2, 2011

# ask user for file names

fin_name = raw_input('Name of input file? ')
fin = open(fin_name, 'r')

fsnp_name = raw_input('Name of snp output file? ')
fsnp = open(fsnp_name, 'w')

findel_name = raw_input('Name of indel output file? ')
findel = open(findel_name, 'w')

# read each line of input

header_section = True

lines_processed = 0
for line in fin:
    #progress output
    if lines_processed % 10000 == 0  :
        print "Lines Processed: %d\n" % lines_processed

    lines_processed += 1
    field_list = line.split()
    #field_list = str.split(line)

    if header_section:
        fsnp.write(line)
        findel.write(line)
        if len(field_list) > 3 and field_list[3] == 'REF':
            print ('End of header section\n')
            header_section = False
    #end header section

    elif len(field_list) > 3 and len(field_list[3]) == 1 and len(field_list[4]) == 1:

        # good line, send to output file
        fsnp.write(line)

    else:
        findel.write(line)


# tell user we are done

print ("All done")



