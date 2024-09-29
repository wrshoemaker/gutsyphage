
import pickle
import sys
import config

import numpy

import json

import pypangraph


annotation_dict_path = '%sannotation_dict.pickle' % config.data_directory

annotation_dict = pickle.load(open(annotation_dict_path, "rb"))


base_table = {'A':'T','T':'A','G':'C','C':'G'}


def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])




def load_pangraph_data(pangraph_file_path):

    pangraph_file = open(pangraph_file_path, 'r')
    pangraph_data = json.load(pangraph_file)
    pangraph_file.close()

    return pangraph_data



class classFASTA:

    # class to load FASTA file

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.frn') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list
    




def parse_pangraph(votu):

    sys.stderr.write("Identifying fourfold status of mutations in %s....\n" % votu)

    # get fasta files for the votu
    fasta_all_genomes = classFASTA('%suhgv_mgv_otu_fna/%s.fna' % (config.data_directory, votu)).readFASTA()
    fasta_genome_dict = {x[0]:x[1] for x in fasta_all_genomes}

    # get fasta files for the votu
    pangraph_file_path = '%scomplete_minimap2/%s_complete_polished.json' % (config.data_directory, votu)
    pangraph_votu = load_pangraph_data(pangraph_file_path)
    pan = pypangraph.Pangraph.load_json(pangraph_file_path)
    loc = pypangraph.Locator(pan)


    block_size_equal = []
    block_size_not_equal = []

    syn_sites_dict = {}
    # loop through blocks
    for block_i in pangraph_votu['blocks']:

        block_i_id = block_i['id']

        if block_i_id not in syn_sites_dict:   
            syn_sites_dict[block_i_id] = {}
            syn_sites_dict[block_i_id]['data'] = {}
            syn_sites_dict[block_i_id]['keep_block'] = True

        
        bl = pan.blocks[block_i_id]
        aln = bl.alignment
        sequences, occurrences = aln.generate_sequences()
        occurrences_genome_names = [x[0] for x in occurrences]

        # looping through genomes
        syn_sites_dict[block_i_id]['data']['genomes'] = {}
        for j in range(0, len(block_i['positions'])):
            
            positions_i_j = block_i['positions'][j]
            genome_name_i_j = positions_i_j[0]['name']

            # strand refers to the orientation of the fragment of genome j that has block_i
            strand_i_j = positions_i_j[0]['strand']
            # positions refer to the positions of the block in the FASTA file
            min_position_i_j, max_position_i_j = positions_i_j[1]

            fasta_genome_j = fasta_genome_dict[genome_name_i_j]

            # returns the sequence of the block for the genome with insertions and deletions added
            # all items in sequences are ordered with respect to the block sequence
            # AKA the reverse complement has not been calculated
            block_sequence_j = sequences[occurrences_genome_names.index(genome_name_i_j)]

            # reverse complement
            # this sequence should match its sequence in the fasta file
            if strand_i_j == False:
                block_sequence_j = calculate_reverse_complement_sequence(block_sequence_j)


            if max_position_i_j > min_position_i_j:
                len_position_i_j = max_position_i_j - min_position_i_j + 1
            else:
                len_position_i_j = (len(fasta_genome_j) - min_position_i_j) + max_position_i_j + 1


            if (len_position_i_j !=  len(block_sequence_j)):
                #print('Reconstituted sequence not equal to positions!')
                sys.stderr.write("Length in pangraph positions not equal to length in FASTA!\n")
               
                # ignore block if we cant examine it in a single genome
                syn_sites_dict[block_i_id]['keep_block'] = False
                block_size_not_equal.append(len_position_i_j)

            else:
                #ys.stderr.write("Length in pangraph positions equal to length in FASTA!\n")
                syn_sites_dict[block_i_id]['keep_block'] = True
                block_size_equal.append(len_position_i_j)


    sys.stderr.write("%d equal blocks with mean block size = %0.2f\n" % (len(block_size_equal), numpy.mean(block_size_equal)))
    sys.stderr.write("%d unequal blocks with mean block size = %0.2f\n" % (len(block_size_not_equal), numpy.mean(block_size_not_equal)))



parse_pangraph('vOTU-000002')
