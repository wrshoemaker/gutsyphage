import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle
import gzip

import random
#from collections import Counter
from itertools import combinations

random.seed(123456789)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import pypangraph


#votu = 'vOTU-000020'
#target_genome = 'UHGV-1365760'


annotation_dict_path = '%sannotation_dict.pickle' % config.data_directory
syn_sites_dict_path = config.data_directory + 'syn_sites_dict_all/%s.pickle' 



genomes_to_ignore = ['UHGV-1372371']





def parse_annotation_table():

    #annotation_table_path = '%stop40_all_metadata.tbl' % config.data_directory
    annotation_table_path = '%stop40_all_metadata.tbl.gz' % config.data_directory
    
    #with gzip.open('%stop40_all_metadata.tbl.gz' % config.data_directory, "rt") as handle:

    #annotation_table_file = open(annotation_table_path, 'r')
    annotation_table_file = gzip.open('%stop40_all_metadata.tbl.gz' % config.data_directory, "rt")


    #header = annotation_table_file.readline()
    annotation_dict = {}
    for line in annotation_table_file:
        line = line.strip()

        if '>' in line:
            genome = line.split(' ')[-1]
            #genome = line.split('\t')[-1]
            annotation_dict[genome] = {}
            annotation_dict[genome]['start_all'] = []
            annotation_dict[genome]['stop_all'] = []
            annotation_dict[genome]['type_all'] = []
            
            annotation_dict[genome]['product_all'] = []
            annotation_dict[genome]['locus_tag_all'] = []
            annotation_dict[genome]['transl_table_all'] = []

        else:
            line_split = line.split('\t')

            if 'CDS' in line_split:
                annotation_dict[genome]['start_all'].append(int(line_split[0]))
                annotation_dict[genome]['stop_all'].append(int(line_split[1]))
                annotation_dict[genome]['type_all'].append(line_split[2])

            if 'product' in line_split:
                annotation_dict[genome]['product_all'].append(line_split[1])

            if 'locus_tag' in line_split:
                annotation_dict[genome]['locus_tag_all'].append(line_split[1])

            if 'transl_table' in line_split:
                annotation_dict[genome]['transl_table_all'].append(int(line_split[1]))


    # save dictionary
    sys.stderr.write("Saving dictionary...\n")
    with open(annotation_dict_path, 'wb') as handle:
        pickle.dump(annotation_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")





def make_syn_sites_votu_dict(votu):
    
    # annotation dict is of all genomes
    annotation_dict = pickle.load(open(annotation_dict_path, "rb"))

    # all translation code 11
    #for key in annotation_dict.keys():
    #    if set(annotation_dict[key]['transl_table_all']) != {11}:
    #        print(set(annotation_dict[key]['transl_table_all']))

    sys.stderr.write("Identifying fourfold status of mutations in %s....\n" % votu)

    # get fasta files for the votu
    fasta_all_genomes = data_utils.classFASTA('%suhgv_mgv_otu_fna/%s.fna' % (config.data_directory, votu)).readFASTA()
    fasta_genome_dict = {x[0]:x[1] for x in fasta_all_genomes}
    #print(fasta_genome_dict.keys())

    # get fasta files for the votu
    pangraph_file_path = '%scomplete_minimap2/%s_complete_polished.json' % (config.data_directory, votu)
    pangraph_votu = data_utils.load_pangraph_data(pangraph_file_path)
    pan = pypangraph.Pangraph.load_json(pangraph_file_path)
    loc = pypangraph.Locator(pan)

    start_stop_genome_dict = {}
    # numpy arrays for ease
    for genome_ in fasta_genome_dict.keys():

        #if genome_ == 'UHGV-1363520':
        #    continue
        #print(genome_)
        start_stop_genome_dict[genome_] = {}
        start_all = numpy.asarray(annotation_dict[genome_]['start_all'])
        stop_all = numpy.asarray(annotation_dict[genome_]['stop_all'])
        
        start_stop_matrix = numpy.vstack((start_all, stop_all))
        max_position_over_genes = numpy.max(start_stop_matrix, axis=0)
        min_position_over_genes = numpy.min(start_stop_matrix, axis=0)
        
        start_stop_genome_dict[genome_]['start_all'] = start_all
        start_stop_genome_dict[genome_]['stop_all'] = stop_all
        start_stop_genome_dict[genome_]['max_position_over_genes'] = max_position_over_genes
        start_stop_genome_dict[genome_]['min_position_over_genes'] = min_position_over_genes


    syn_sites_dict = {}
    #syn_sites_dict['keep_block'] 
    syn_sites_dict['data'] = {}
    # loop through blocks
    for block_i in pangraph_votu['blocks']:

        block_i_id = block_i['id']

        if block_i_id not in syn_sites_dict:   
            syn_sites_dict[block_i_id] = {}
            syn_sites_dict[block_i_id]['data'] = {}
            syn_sites_dict[block_i_id]['keep_block'] = True

        
        block_i_mutate = block_i['mutate']

        bl = pan.blocks[block_i_id]
        aln = bl.alignment
        sequences, occurrences = aln.generate_sequences()
        occurrences_genome_names = [x[0] for x in occurrences]

        # looping through genomes
        for j in range(0, len(block_i['positions'])):
            
            positions_i_j = block_i['positions'][j]
            genome_name_i_j = positions_i_j[0]['name']

            if genome_name_i_j in genomes_to_ignore:
                continue

            # strand refers to the orientation of the fragment of genome j that has block_i
            strand_i_j = positions_i_j[0]['strand']
            # positions refer to the positions of the block in the FASTA file
            min_position_i_j, max_position_i_j = positions_i_j[1]

            #if genome_name_i_j != target_genome:
            #    continue

            start_all_j = start_stop_genome_dict[genome_name_i_j]['start_all']
            stop_all_j = start_stop_genome_dict[genome_name_i_j]['stop_all']
            max_position_over_genes_j = start_stop_genome_dict[genome_name_i_j]['max_position_over_genes']
            min_position_over_genes_j = start_stop_genome_dict[genome_name_i_j]['min_position_over_genes']
            fasta_genome_j = fasta_genome_dict[genome_name_i_j]

            # returns the sequence of the block for the genome with insertions and deletions added
            # all items in sequences are ordered with respect to the block sequence
            # AKA the reverse complement has not been calculated
            block_sequence_j = sequences[occurrences_genome_names.index(genome_name_i_j)]

            # reverse complement
            if strand_i_j == False:
                block_sequence_j = data_utils.calculate_reverse_complement_sequence(block_sequence_j)


            if max_position_i_j > min_position_i_j:
                len_position_i_j = max_position_i_j - min_position_i_j + 1
            else:
                len_position_i_j = (len(fasta_genome_j) - min_position_i_j) + max_position_i_j + 1


            if (len_position_i_j !=  len(block_sequence_j)):
                #print('Reconstituted sequence not equal to positions!')
                #sys.stderr.write("Reconstituted sequence not equal to positions!\n")
                # ignore block if we cant examine it in a single genome
                syn_sites_dict[block_i_id]['keep_block'] = False
                continue
            #else:
            #    print('Okay!')


            block_i_mutate_genome_j = [x[1] for x in block_i_mutate if x[0]['name'] == genome_name_i_j][0]

            # no mutations
            if len(block_i_mutate_genome_j) == 0:
                syn_sites_dict[block_i_id]['keep_block'] = False
                continue


            # If the block wraps circular we do two searches
            # The one at the end of molecule, and the second at the beginning of molecule
            # when min_position_i_j > max_position_i_j the block goes from the end of the fasta to the beginning
            # we are only looking at **complete** genes within a fragment
            if (min_position_i_j > max_position_i_j):
                
                #new_max_position = max(stop_all)
                new_max_position = len(fasta_genome_j)
                # counting starts at one for gene positions
                new_min_position = 1

                # get genes from min_position_i_j to end of genomes
                genes_1_idx = (max_position_over_genes_j <= new_max_position) & (min_position_over_genes_j >= min_position_i_j)
                genes_2_idx = (max_position_over_genes_j <= max_position_i_j) & (min_position_over_genes_j >= new_min_position)
                idx_to_keep = (genes_1_idx | genes_2_idx)


            else:
                # Identify genes are contained by the block
                idx_to_keep = (max_position_over_genes_j <= max_position_i_j) & (min_position_over_genes_j >= min_position_i_j)
            

            min_position_to_keep = start_all_j[idx_to_keep]
            max_position_to_keep = stop_all_j[idx_to_keep]

            # skip if there are no genes in this block for this genome
            if (len(min_position_to_keep) == 0):
                continue

            # we have the positions of the genes in this block for this genome.
            # find the positions of the synonymous sites            
            # loop through genes in the chunk
            sites_synon_status_all = []
            sites_position_in_fasta_all = []
            for gene_g_idx in range(len(min_position_to_keep)):

                min_position_g = min_position_to_keep[gene_g_idx]
                max_position_g = max_position_to_keep[gene_g_idx]

                # get sequence from fasta
                # reverse complement
                if min_position_g > max_position_g:
                    # check this with real annotated genomes
                    seq_g = fasta_genome_j[max_position_g-1:min_position_g]
                    # reorient sequence
                    gene_sequence_g = data_utils.calculate_reverse_complement_sequence(seq_g)
                    # swap order
                    min_position_g, max_position_g = max_position_g, min_position_g

                else:
                    # gene positions start counting at one
                    # index at zero for python
                    gene_sequence_g = fasta_genome_j[min_position_g-1:max_position_g]

                # loop among positions within the gene
                for position_in_gene in list( range(len(gene_sequence_g))): #calculate codon start
                    # start position of codon
                    codon_start = int((position_in_gene)/3)*3
                    codon = gene_sequence_g[codon_start:codon_start+3] 
                    position_in_codon = position_in_gene%3
                    
                    #data_utils.codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                    sites_synon_status_all.append(data_utils.codon_synonymous_opportunity_table[codon][position_in_codon])
                    sites_position_in_fasta_all.append(min_position_g + position_in_gene)


            # the positions of each fasta position in the current block
            sites_position_in_block_all = [loc.find_position(strain=genome_name_i_j, pos=h)[1] for h in sites_position_in_fasta_all]
            # the positions and allele of each mutation in the current block that are in CDS regions
            block_i_mutate_genome_j_in_cds = [x for x in block_i_mutate_genome_j if x[0] in sites_position_in_block_all ]

            # no mutations in CDS regions! skip.
            if len(block_i_mutate_genome_j_in_cds) == 0:
                continue

            # splitting positions and allele
            block_i_mutate_genome_j_in_cds_block_positions = [x[0] for x in block_i_mutate_genome_j_in_cds]
            block_i_mutate_genome_j_in_cds_block_alleles = [x[1] for x in block_i_mutate_genome_j_in_cds]
            
            # the position of each mutation in the fasta file
            block_i_mutate_genome_j_in_cds_fasta_positions = [sites_position_in_fasta_all[sites_position_in_block_all.index(s)] for s in block_i_mutate_genome_j_in_cds_block_positions]
            block_i_mutate_genome_j_in_cds_syn_status = [sites_synon_status_all[sites_position_in_block_all.index(s)] for s in block_i_mutate_genome_j_in_cds_block_positions]
            
            # sites_synon_status has the positions of the sites in the gene with respect to the fasta file
            syn_sites_dict[block_i_id]['data']['n_sites'] = {}
            # sites_synon_status_all
            sites_synon_status_all = numpy.asarray(sites_synon_status_all)
            # count number of sites in block for each fourfold status
            for fourfold_status in range(4):
                syn_sites_dict[block_i_id]['data']['n_sites'][fourfold_status] = sum(sites_synon_status_all == fourfold_status)
            
            syn_sites_dict[block_i_id]['data']['genomes'] = {}
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j] = {}
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j]['mutate_block_position'] = block_i_mutate_genome_j_in_cds_block_positions
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j]['mutate_allele'] = block_i_mutate_genome_j_in_cds_block_alleles
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j]['mutate_fasta_position'] = block_i_mutate_genome_j_in_cds_fasta_positions
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j]['mutate_syn_status'] = block_i_mutate_genome_j_in_cds_syn_status

            # strand of the block on the fasta file
            syn_sites_dict[block_i_id]['data']['genomes'][genome_name_i_j]['strand'] = strand_i_j


    # save dictionary
    sys.stderr.write("Saving dictionary...\n")
    with open(syn_sites_dict_path % votu, 'wb') as handle:
        pickle.dump(syn_sites_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")





def run_everything():

    #parse_annotation_table()

    votu_all = [x.split('_')[0] for x in os.listdir('%scomplete_minimap2/' % config.data_directory)]
    votu_all.sort()
    #votu_all = numpy.asarray(votu_all)
    proble_votu_idx = votu_all.index('vOTU-000118')
    #data_utils.build_votu_fasta('vOTU-000050')

    # votus_to_skip = ['vOTU-000118', 'vOTU-000016', vOTU-000035, vOTU-000118]
    
    for votu_idx, votu in enumerate(votu_all[proble_votu_idx:]):
        # make fasta file
        if votu == 'vOTU-000118':
            continue
        #data_utils.build_votu_fasta(votu)
        
        # make syn dict
        make_syn_sites_votu_dict(votu)

    #make_syn_sites_votu_dict('vOTU-000016')




#alignment, occurrences = aln.generate_alignments()

run_everything()







