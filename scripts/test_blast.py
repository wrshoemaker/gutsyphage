

import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config


import gzip
from Bio import SeqIO



target_votu = 'vOTU-000062'
votu_dict = data_utils.read_phage_metadata()

genomes = votu_dict[target_votu]['uhgv_genome']

target_genome_all = genomes[:2]


fasta_object = data_utils.classFASTA('%suhgv_mgv.fna.gz' % config.data_directory).readFASTA()


fasta_file_path = '%stest.fna' % (config.data_directory) 
fasta_file_open = open(fasta_file_path, 'w')

with gzip.open('%suhgv_mgv.fna.gz' % config.data_directory, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):

        if record.id in target_genome_all:

            fasta_file_open.write('>%s\n' % record.id)

            clean_sites_seq_split = [record.seq[i:i+data_utils.n_fna_characters] for i in range(0, len(record.seq), data_utils.n_fna_characters)]

            for seq in clean_sites_seq_split:
                fasta_file_open.write('%s\n' % seq)

            fasta_file_open.write('\n')



fasta_file_open.close()



in_file='/Users/williamrshoemaker/GitHub/gutsyphage/data/test.fna'
db_name='/Users/williamrshoemaker/GitHub/gutsyphage/data/uhgv_mgv_blast_db/uhgv_mgv'
out_file='/Users/williamrshoemaker/GitHub/gutsyphage/data/test_blast'

# ./ncbi-blast-2.16.0+/bin/makeblastdb -in ../data/test.fna -out ../data/single_db/votu3_example -dbtype nucl
# ./ncbi-blast-2.16.0+/bin/blastn -query $file -db $db_name -out $out_file -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
#./ncbi-blast-2.16.0+/bin/blastn -query $in_file -db $db_name -out $out_file -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90



