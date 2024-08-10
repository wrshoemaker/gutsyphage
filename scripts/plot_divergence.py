import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle

import random
#from collections import Counter
from itertools import combinations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



#divergence_dict_path = '%sdivergence_dict.pickle' % config.data_directory








def calculate_divergence_all_votus(n_genomes_subsample=6):

    # make dictionary of divergences 
    file_directory = '%scomplete_minimap2/' % config.data_directory
    dict_directory = '%sdivergence_dict_all/' % config.data_directory

    for filename in os.listdir(file_directory):
        f = os.path.join(file_directory, filename)
        # checking if it is a file
        if os.path.isfile(f):

            pangraph_data = data_utils.load_pangraph_data(f)
            pangraph_genome_names = data_utils.get_pangraph_genome_names(pangraph_data)
            pangraph_genome_names.sort()

            genome_pair_all = list(combinations(pangraph_genome_names,2))

            votu = filename.split('/')[-1].split('.')[0].split('_')[0]


            print(votu, len(pangraph_genome_names))
            div_dict_path = '%s%s.pickle' % (dict_directory, votu)

            div_dict = {}
            for genome_pair in genome_pair_all:

                bins, binned_divergence, total_divergence, len_fraction_shared_blocks, len_fraction_shared_blocks_union = data_utils.calculate_divergence_across_pangraph_blocks(genome_pair[0], genome_pair[1], pangraph_data, calculate_binned_divergence=False)

                if total_divergence == None:
                    continue

                div_dict[genome_pair] = {}
                #div_dict[genome_pair]['bins'] = bins.tolist()
                #div_dict[genome_pair]['binned_divergence'] = binned_divergence.tolist()
                div_dict[genome_pair]['total_divergence'] = total_divergence
                div_dict[genome_pair]['len_fraction_shared_blocks'] = len_fraction_shared_blocks
                div_dict[genome_pair]['len_fraction_shared_blocks_union'] = len_fraction_shared_blocks_union

    
            sys.stderr.write("Saving dictionary...\n")
            with open(div_dict_path, 'wb') as handle:
                pickle.dump(div_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            sys.stderr.write("Done!\n")




def plot_divergence_vs_shared_blocks():

    votu_file = '%stop_20_votus_per_lifestyle.csv' % config.data_directory
    votu_open = open(votu_file, 'r')

    header = votu_open.readline()
    header_split = header.strip().split(',')

    lifestyle_all = []
    votu_all = []
    for line in votu_open:
        line_split = line.strip().split(',')
        lifestyle_all.append(line_split[0])
        votu_all.append(line_split[1])

    votu_open.close()

    votu_all.sort()


    votu_dict, vgenome_dict = data_utils.read_uhgv_metadata()
    uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()
    sample_metagenome_dict = data_utils.read_sample_metadata()

    #votu = 'vOTU-000048'

    for votu in votu_all:

        print(votu)

        div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)

        dict_ = pickle.load(open(div_dict_path, "rb"))
        genome_pairs = list(dict_.keys())
        total_divergence_all = numpy.asarray([dict_[g]['total_divergence'] for g in genome_pairs])
        len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks'] for g in genome_pairs])
        lifestyle_votu = uhgv_votu_metadata_dict[votu]['lifestyle']



        #to_keep_idx = [True if (vgenome_dict[k[0]]['original_id'] in sample_metagenome_dict) and (vgenome_dict[k[1]]['original_id'] in sample_metagenome_dict) else False for k in genome_pairs]
        to_keep_idx = [True if (k[0] in vgenome_dict) and (k[1] in vgenome_dict) else False for k in genome_pairs]
        genome_pairs_clean = [genome_pairs[g] for g in range(len(genome_pairs)) if to_keep_idx[g] == True]

        #genome_pairs = genome_pairs[to_keep_idx]
        total_divergence_all = total_divergence_all[to_keep_idx]
        len_fraction_shared_blocks_all = len_fraction_shared_blocks_all[to_keep_idx]


        same_country_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['country_code'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['country_code'] else False for k in genome_pairs_clean]
        same_continent_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['continent'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['continent'] else False for k in genome_pairs_clean]

        same_country_idx = numpy.asarray(same_country_idx)
        same_continent_idx = numpy.asarray(same_continent_idx)


        # different continent ~same_continent_idx
        # same continent different country (same_continent_idx) & (~same_country_idx)
        # same continent same country (same_continent_idx) & (same_country_idx)
        

        same_continent_diff_country = (same_continent_idx) & (~same_country_idx)
        same_continent_same_country = (same_continent_idx) & (same_country_idx)

        fig = plt.figure(figsize = (4, 4))
        fig.subplots_adjust(bottom= 0.15)

        ax = plt.subplot2grid((1, 1), (0, 0))
        ax.scatter(total_divergence_all[~same_continent_idx], len_fraction_shared_blocks_all[~same_continent_idx], s=8, c='#FF6347', alpha=0.3, label='Diff. continent')
        ax.scatter(total_divergence_all[same_continent_diff_country], len_fraction_shared_blocks_all[same_continent_diff_country], s=8, c='#FFA500', alpha=0.3, label='Same continent, diff. country')
        ax.scatter(total_divergence_all[same_continent_same_country], len_fraction_shared_blocks_all[same_continent_same_country], s=8, c='#87CEEB', alpha=0.3, label='Same country')

        ax.set_xlabel('Genome-wide pairwise divergence on shared blocks', fontsize=10)
        ax.set_ylabel('Length fraction of shared blocks', fontsize=10)
        ax.set_title('%s\nLifestyle = %s' % (votu, lifestyle_votu), fontsize=12)

        ax.legend(loc='lower left', fontsize=8)

        fig.subplots_adjust(hspace=0.45, wspace=0.45)
        fig_name = "%sdivergence/%s.png" % (config.analysis_directory, votu)
        fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
        plt.close()




def plot_shared_blocks_dist(min_len=0.05):

    votu_file = '%stop_20_votus_per_lifestyle.csv' % config.data_directory
    votu_open = open(votu_file, 'r')

    header = votu_open.readline()
    header_split = header.strip().split(',')

    lifestyle_all = []
    votu_all = []
    for line in votu_open:
        line_split = line.strip().split(',')
        lifestyle_all.append(line_split[0])
        votu_all.append(line_split[1])

    votu_open.close()
    
    votu_dict, vgenome_dict = data_utils.read_uhgv_metadata()
    uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()
    sample_metagenome_dict = data_utils.read_sample_metadata()

    #votu_all.sort()
    votu_all = numpy.asarray(votu_all)
    lifestyle_votu_all = numpy.asarray([uhgv_votu_metadata_dict[votu]['lifestyle'] for votu in votu_all])

    #votu_all = votu_all[lifestyle_votu_all.argsort()]

    len_fraction_shared_blocks_all_votu_dict = {}
    mean_len_fraction_shared_blocks_all_votu = []

    for votu in votu_all:

        div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
        dict_ = pickle.load(open(div_dict_path, "rb"))
        genome_pairs = list(dict_.keys())
        total_divergence_all = numpy.asarray([dict_[g]['total_divergence'] for g in genome_pairs])
        len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks'] for g in genome_pairs])
        len_fraction_shared_blocks_all = len_fraction_shared_blocks_all[len_fraction_shared_blocks_all >= min_len]

        #len_fraction_shared_blocks_all_votu.append(len_fraction_shared_blocks_all)
        len_fraction_shared_blocks_all_votu_dict[votu] = len_fraction_shared_blocks_all
        mean_len_fraction_shared_blocks_all_votu.append(numpy.mean(len_fraction_shared_blocks_all))


    mean_len_fraction_shared_blocks_all_votu = numpy.asarray(mean_len_fraction_shared_blocks_all_votu)
    
    fig = plt.figure(figsize = (4, 8))
    fig.subplots_adjust(bottom= 0.15)
    ax = plt.subplot2grid((1, 1), (0, 0))

    #print(mean_len_fraction_shared_blocks_all_votu.argsort())


    n_rows = 0 
    votu_all_ordered = []
    
    for l in data_utils.lifestyle_all:

        l_idx = (lifestyle_votu_all==l)

        mean_len_fraction_shared_blocks_all_votu_l = mean_len_fraction_shared_blocks_all_votu[l_idx]
        votu_all_l = votu_all[l_idx]

        votu_all_l = votu_all_l[mean_len_fraction_shared_blocks_all_votu_l.argsort()]

        
        for votu_idx, votu in enumerate(votu_all_l):

            #div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
            #dict_ = pickle.load(open(div_dict_path, "rb"))
            #genome_pairs = list(dict_.keys())
            #total_divergence_all = numpy.asarray([dict_[g]['total_divergence'] for g in genome_pairs])
            #len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks'] for g in genome_pairs])
            #len_fraction_shared_blocks_all = len_fraction_shared_blocks_all[len_fraction_shared_blocks_all >= min_len]

            len_fraction_shared_blocks_all = len_fraction_shared_blocks_all_votu_dict[votu]

            y_jitter = [n_rows]*len(len_fraction_shared_blocks_all)
            y_jitter = y_jitter + numpy.random.randn(len(y_jitter)) * 0.1
            #x = rand_jitter(numpy.asarray([votu_idx]*))
            ax.scatter(len_fraction_shared_blocks_all, y_jitter, s=2, alpha=0.2, c=data_utils.lifestyle_color_dict[l])

            votu_all_ordered.append(votu)
            n_rows += 1


    ax.set_yticks(list(range(len(votu_all))))
    ax.set_yticklabels(votu_all_ordered, fontsize=6, ha="right")
    ax.set_xlabel('Length fraction of shared blocks', fontsize=12)

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['lytic'], label='Lytic', markersize=15),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['temperate'], label='Temperate', markersize=15)]

    ax.legend(handles=legend_elements, loc='lower right')



    fig.subplots_adjust(hspace=0.45, wspace=0.45)
    fig_name = "%sshared_blocks_dist.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



    

#plot_divergence_vs_shared_blocks()

#plot_shared_blocks_dist()




calculate_divergence_all_votus()
#plot_shared_blocks_dist()
#




#votu='vOTU-000002'

#pangraph_data = load_pangraph_data(votu)

#pangraph_genome_names = get_pangraph_genome_names(pangraph_data)

#genome_1 = pangraph_genome_names[0]
#genome_2 = pangraph_genome_names[1]


#fig = plt.figure(figsize = (8, 4))
#fig.subplots_adjust(bottom= 0.15)

#ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)


#ax.plot(bins[:-1], binned_divergence)


#fig.subplots_adjust(hspace=0.3, wspace=0.25)
#fig_name = "%sdivergence%s_%s.png" % (config.analysis_directory, genome_1, genome_2)
#fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
#plt.close()

#print(len(bins), len(binned_divergence))



