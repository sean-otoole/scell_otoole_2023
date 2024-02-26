# fasta construction

import os
from Bio import SeqIO

file_path = 'XXXPathXXXX/addgene-plasmid-101060-sequence-194341.gbk'

from Bio import SeqIO
gb_file = file_path
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
    # now do something with the record
    print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
    print(repr(gb_record.seq))
    
### just get the sequence data

campari2_raw_seq = str(getattr(gb_record,'seq'))
print(len(campari2_raw_seq))

### now we'll read the original fasta file
os.chdir('XXXPathXXXX/original_reference')

with open ("genome.fa", "r") as myfile:
    genome_list=myfile.readlines()

### use list comprehension to find the headers

matched_items = [item for item in genome_list if ">" in item]

#### I'll use the mitochondrial genome as a model
template_header = matched_items[20]

virus_header = template_header.replace('X', 'V')
virus_header = virus_header.replace('171031299', str(len(campari2_raw_seq)))

### next step is the break things up into 61 character increments for each line

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

broken_dna = list(chunkstring(campari2_raw_seq, 61))  ### this is where we've stored our list

## adds new line character to end of each new list item

broken_dna = ['{0}\n'.format(element) for element in broken_dna]

### Now let's place our sequence list after the header

fasta_input_list = []
fasta_input_list.extend(broken_dna) 
fasta_input_list.insert(0, virus_header)

### now want to create a copy of the fasta file so the original is untarnished
### in case mistakes are made

from shutil import copyfile

fileDest = 'XXXPathXXXX/genome.custom.fa'

copyfile('XXXPathXXXX/genome.fa', fileDest)

## change directory for new genome

os.chdir('XXXPathXXXX/custom_reference')

### now we write the list to end of our custom genome file

with open("genome.custom.fa", "a") as myfile:
    for item in fasta_input_list:
        myfile.write(item)

## GTF writer

import os
os.chdir('XXXPathXXXX/custom_reference')

def readLine(lineNumber,file):
    with open(file) as f:
        line = f.readlines()
        output = (line[lineNumber])
        return(output)
        f.close()

### This is the last line of the GTF file

x = readLine(6,'genes.gtf')

### This is a list for CaMPARI2 that will be written to the end of the GTF

### first 8 components are seqname,source,feature,start,end,score,strand, and frame, they need
### to be tab delimited

### all other elements are features, which can vary between genes, they will be semi-colon
### delimited

###I've excluded the calmodulin domain from the reads

campari2 = ['V','Campari2','exon','1','3500','.','+','.',  ###first 8 elements, tab delimited
            'gene_id "Campari2"',   ### remaining elements, semicolon delimited
            ' gene_version "1"',
            ' transcript_id "Campari2"',
            ' transcript_version "1"',
            ' gene_name "Campari2"',
            ' gene_source "addgene"',
            ' gene_biotype "protein_coding"',
            ' transcript_name "Campari2"',
            ' transcript_source "custom"',
            ' transcript_biotype "protien_coding"',
            ' protein_id "Campari2"',
            ' tag "basic"']

### use the above list to construct a new line for the custom gtf

campari2_line_part1 = '\t'.join((campari2[0:8]))
campari2_line_part2 = ';'.join((campari2[8:20]))
campari2_line_final = '\t'.join((campari2_line_part1,campari2_line_part2))  ## there is one tab separation

### now we want to write this line to end of a new custom GTF file
### we will still preserve the old file in case any mistakes are made
### and we need to try again

from shutil import copyfile

fileDest = 'XXXPathXXXX/custom_reference/genes.custom.gtf'

copyfile('XXXPathXXXX/original_reference/genes.gtf', fileDest)

## lets change the desintation away from opt for convenience

os.chdir('XXXPathXXXX/campari2_genome/custom_reference')

## Now we'll write to our copied file

with open("genes.custom.gtf", "a") as myfile:
    myfile.write(campari2_line_final)





