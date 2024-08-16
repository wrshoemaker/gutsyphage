import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle
import time
import operator
#import logger


import random
#from collections import Counter
from itertools import combinations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from numpy.random import normal
import matplotlib as mpl
import pylab

from matplotlib.patches import Patch


import multiprocessing



random.seed(123456789)


#divergence_dict_path = '%sdivergence_dict.pickle' % config.data_directory

votu_dict, vgenome_dict = data_utils.read_uhgv_metadata()
uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()
sample_metagenome_dict = data_utils.read_sample_metadata()
mgv_uhgv_species_dict = data_utils.parse_mgv_uhgv_species()


#file_directory = '%scomplete_minimap2/' % config.data_directory
div_dict_path_template = config.data_directory + 'divergence_dict_all/%s.pickle'
#div_dict_directory = '%sdivergence_dict_all/' % config.data_directory
syn_sites_path_template =  config.data_directory + 'syn_sites_dict_all/%s.pickle'

#syn_sites_path = '%s%s.pickle' % (syn_sites_directory, votu)
minimap_path =  config.data_directory + 'complete_minimap2/%s_complete_polished.json'




def calculate_divergence(votu, n_pairs=None):#, n_genomes_subsample=6):


    # make dictionary of divergences 
    
    #for filename in os.listdir(file_directory):
        
    #f = os.path.join(file_directory, filename)
    #minimap_path_votu = os.path.join(minimap_path, votu)
    minimap_path_votu = minimap_path % votu
    # checking if it is a file
    if os.path.isfile(minimap_path_votu):

        pangraph_data = data_utils.load_pangraph_data(minimap_path_votu)
        pangraph_genome_names = data_utils.get_pangraph_genome_names(pangraph_data)
        pangraph_genome_names.sort()

        genome_pair_all = list(combinations(pangraph_genome_names,2))

        #votu = filename.split('/')[-1].split('.')[0].split('_')[0]

        print(votu)
        
        # check if fourfold status is detemrined
        #syn_sites_path = '%s%s.pickle' % (syn_sites_directory, votu)
        syn_sites_path = syn_sites_path_template % votu
        if os.path.exists(syn_sites_path):
            syn_sites_dict = pickle.load(open(syn_sites_path, "rb"))
        else:
            syn_sites_dict = None


        # make it a numpy array
        if syn_sites_dict != None:
            for block_i_id, block_i_id_dict in syn_sites_dict.items():

                if block_i_id == 'data':
                    continue
                            # only consider blocks 
                if syn_sites_dict[block_i_id]['keep_block'] == True:

                    genome_all = list(block_i_id_dict['data']['genomes'].keys())
                    for g in genome_all:

                        syn_sites_dict[block_i_id]['data']['genomes'][g]['site_block_position'] = numpy.asarray(block_i_id_dict['data']['genomes'][g]['site_block_position'])
                        syn_sites_dict[block_i_id]['data']['genomes'][g]['site_syn_status'] = numpy.asarray(block_i_id_dict['data']['genomes'][g]['site_syn_status'])



        div_dict_path = div_dict_path_template % votu
        div_dict = {}
        #n_genome_pair = len(genome_pair_all)
        #updates = 0

        if n_pairs != None:
            random.shuffle(genome_pair_all)
            genome_pair_all = genome_pair_all[:n_pairs]
            
            #genome_pair_subsample_idx_all = range(0, len(genome_pair_all), n_pairs)
        
        for genome_pair_idx, genome_pair in enumerate(genome_pair_all):

            #finished = 100*(genome_pair_idx/n_genome_pair)

            if genome_pair_idx % 5000 == 0:
                print(genome_pair_idx)
            
            #if divmod(finished, 10) == (updates, 0):
            #    updates += 1
            #    #if int(finished) == 0:
            #    #    continue
            #    sys.stderr.write(str(int(finished)) + "%" + " done...\n")

            bins, binned_divergence, total_divergence, cumulative_n_nonsyn, cumulative_n_syn, cumulative_block_len_nonsyn, cumulative_block_len_syn, len_fraction_shared_blocks, len_fraction_shared_blocks_union = data_utils.calculate_divergence_across_pangraph_blocks(genome_pair[0], genome_pair[1], pangraph_data, syn_sites_dict=syn_sites_dict, calculate_binned_divergence=False)
            # ignore    
            if total_divergence == None:
                continue

            #dn = (cumulative_n_nonsyn + 1)/(cumulative_block_len_nonsyn + 1)
            #ds = (cumulative_n_syn + 1)/(cumulative_block_len_syn + 1)

            #if (cumulative_n_nonsyn > 3) and (cumulative_n_syn > 3):

            #    print(ds, dn/ds)

            div_dict[genome_pair] = {}
            #div_dict[genome_pair]['bins'] = bins.tolist()
            #div_dict[genome_pair]['binned_divergence'] = binned_divergence.tolist()
            div_dict[genome_pair]['total_divergence'] = total_divergence
            div_dict[genome_pair]['len_fraction_shared_blocks'] = len_fraction_shared_blocks
            div_dict[genome_pair]['len_fraction_shared_blocks_union'] = len_fraction_shared_blocks_union

            div_dict[genome_pair]['cumulative_n_nonsyn'] = cumulative_n_nonsyn
            div_dict[genome_pair]['cumulative_n_syn'] = cumulative_n_syn
            div_dict[genome_pair]['cumulative_block_len_nonsyn'] = cumulative_block_len_nonsyn
            div_dict[genome_pair]['cumulative_block_len_syn'] = cumulative_block_len_syn



        sys.stderr.write("Saving dictionary...\n")
        with open(div_dict_path, 'wb') as handle:
            pickle.dump(div_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        sys.stderr.write("Done!\n")






def plot_divergence_vs_shared_blocks(votu, syn_divergence=True, min_n_muts=30, min_n_sites=1e3):

    #votu_file = '%stop_20_votus_per_lifestyle.csv' % config.data_directory
    #votu_open = open(votu_file, 'r')

    #header = votu_open.readline()
    #header_split = header.strip().split(',')

    #lifestyle_all = []
    #votu_all = []
    #for line in votu_open:
    #    line_split = line.strip().split(',')
    #    lifestyle_all.append(line_split[0])
    #    votu_all.append(line_split[1])

    #votu_open.close()

    #votu_all.sort()

    #for votu in votu_all:

    print(votu)

    div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
    div_dict = pickle.load(open(div_dict_path, "rb"))
    genome_pairs = list(div_dict.keys())
    lifestyle_votu = uhgv_votu_metadata_dict[votu]['lifestyle']

    divergence_all = []
    len_fraction_shared_blocks_all = []
    genome_pairs_clean = []
    for genome_pair_idx, genome_pair in enumerate(genome_pairs):

        if syn_divergence == True:

            cumulative_n_syn = div_dict[genome_pair]['cumulative_n_syn']
            cumulative_block_len_syn = div_dict[genome_pair]['cumulative_block_len_syn']

            if (cumulative_n_syn == None) or (cumulative_block_len_syn == None):
                continue
            
            # not enough data to estimate divergence
            if (cumulative_n_syn < min_n_muts) or (cumulative_block_len_syn < min_n_sites):
                continue

            genome_pair_div = cumulative_n_syn/cumulative_block_len_syn

        else:
            genome_pair_div = div_dict[genome_pair]['total_divergence']


        genome_pairs_clean.append(genome_pair)
        divergence_all.append(genome_pair_div)
        len_fraction_shared_blocks_all.append(div_dict[genome_pair]['len_fraction_shared_blocks_union'])


    #total_divergence_all = numpy.asarray([dict_[g]['total_divergence'] for g in genome_pairs])
    #len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks_union'] for g in genome_pairs])

    divergence_all = numpy.asarray(divergence_all)
    len_fraction_shared_blocks_all = numpy.asarray(len_fraction_shared_blocks_all)

    #to_keep_idx = [True if (vgenome_dict[k[0]]['original_id'] in sample_metagenome_dict) and (vgenome_dict[k[1]]['original_id'] in sample_metagenome_dict) else False for k in genome_pairs]
    to_keep_idx = [True if (k[0] in vgenome_dict) and (k[1] in vgenome_dict) else False for k in genome_pairs_clean]
    genome_pairs_clean = [genome_pairs_clean[g] for g in range(len(genome_pairs_clean)) if to_keep_idx[g] == True]

    #genome_pairs = genome_pairs[to_keep_idx]
    total_divergence_all = divergence_all[to_keep_idx]
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
    ax.scatter(divergence_all[~same_continent_idx], len_fraction_shared_blocks_all[~same_continent_idx], s=8, c='#FF6347', alpha=0.3, label='Diff. continent')
    ax.scatter(divergence_all[same_continent_diff_country], len_fraction_shared_blocks_all[same_continent_diff_country], s=8, c='#FFA500', alpha=0.3, label='Same continent, diff. country')
    ax.scatter(divergence_all[same_continent_same_country], len_fraction_shared_blocks_all[same_continent_same_country], s=8, c='#87CEEB', alpha=0.3, label='Same country')

    if syn_divergence == True:
        x_label = 'Synonymous divergence on shared blocks, ' + r'$dS$'
    else:
        x_label = 'Divergence on shared blocks'

    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel('Length fraction of shared blocks', fontsize=10)
    ax.set_title('%s\nLifestyle = %s' % (votu, lifestyle_votu), fontsize=12)

    ax.set_xscale('log', base=10)

    ax.legend(loc='lower left', fontsize=8)

    fig.subplots_adjust(hspace=0.45, wspace=0.45)
    fig_name = "%sdivergence_vs_shared_blocks/%s.png" % (config.analysis_directory, votu)
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
        len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks_union'] for g in genome_pairs])
        len_fraction_shared_blocks_all = len_fraction_shared_blocks_all[len_fraction_shared_blocks_all >= min_len]

        #len_fraction_shared_blocks_all_votu.append(len_fraction_shared_blocks_all)
        len_fraction_shared_blocks_all_votu_dict[votu] = len_fraction_shared_blocks_all
        mean_len_fraction_shared_blocks_all_votu.append(numpy.mean(len_fraction_shared_blocks_all))


    mean_len_fraction_shared_blocks_all_votu = numpy.asarray(mean_len_fraction_shared_blocks_all_votu)
    
    fig = plt.figure(figsize = (4, 8))
    fig.subplots_adjust(bottom= 0.15)
    ax = plt.subplot2grid((1, 1), (0, 0))


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




def plot_ds_vs_dnds(pseudocount=0, min_n_muts=10, min_n_sites=1e3):

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
    
    votu_all = numpy.asarray(votu_all)
    #lifestyle_votu_all = numpy.asarray([uhgv_votu_metadata_dict[votu]['lifestyle'] for votu in votu_all])


    votu_all = ['vOTU-000018']

    #n_pairs = 0

    for votu in votu_all:

        div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
        div_dict = pickle.load(open(div_dict_path, "rb"))


        fig = plt.figure(figsize = (4, 4))
        fig.subplots_adjust(bottom= 0.15)

        ax = plt.subplot2grid((1, 1), (0, 0))

        ds_all = []
        dnds_all = []
        genome_pairs_clean_all = []

        for genome_pair, genome_pair_dict in div_dict.items():

            cumulative_n_syn = genome_pair_dict['cumulative_n_syn']
            cumulative_n_nonsyn = genome_pair_dict['cumulative_n_nonsyn']

            cumulative_block_len_syn = genome_pair_dict['cumulative_block_len_syn']
            cumulative_block_len_nonsyn = genome_pair_dict['cumulative_block_len_nonsyn']

            if (cumulative_n_syn == None) or (cumulative_n_nonsyn == None) or (cumulative_block_len_syn == None) or (cumulative_block_len_nonsyn == None):
                continue

            # at least one mutation in each 
            if (cumulative_n_syn == 0) or (cumulative_n_nonsyn == 0):
                continue
            
            # total of five mutations
            if (cumulative_n_syn + cumulative_n_nonsyn) <= min_n_muts:
                continue

            # at least min_n_sites possible sites
            if (cumulative_block_len_syn < min_n_sites) or (cumulative_block_len_nonsyn < min_n_sites):
                continue

            dn = (cumulative_n_nonsyn+pseudocount)/(cumulative_block_len_nonsyn+pseudocount)
            ds = (cumulative_n_syn+pseudocount)/(cumulative_block_len_syn+pseudocount)

            ds_all.append(ds)
            dnds_all.append(dn/ds)
            genome_pairs_clean_all.append(genome_pair)


        ds_all = numpy.asarray(ds_all)
        dnds_all = numpy.asarray(dnds_all)

        same_country_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['country_code'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['country_code'] else False for k in genome_pairs_clean_all]
        same_continent_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['continent'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['continent'] else False for k in genome_pairs_clean_all]

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
        ax.scatter(ds_all[~same_continent_idx], dnds_all[~same_continent_idx], s=8, c='#FF6347', alpha=0.3, label='Diff. continent')
        ax.scatter(ds_all[same_continent_diff_country], dnds_all[same_continent_diff_country], s=8, c='#FFA500', alpha=0.3, label='Same continent, diff. country')
        ax.scatter(ds_all[same_continent_same_country], dnds_all[same_continent_same_country], s=8, c='#87CEEB', alpha=0.3, label='Same country')

        ax.set_xlabel('Synonymous divergence, ' + r'$d_{S}$', fontsize=10)
        ax.set_ylabel('Nonsynonymous ratio, ' + r'$d_{N}/d_{S}$', fontsize=10)
        ax.set_title('%s\nLifestyle = %s' % (votu, uhgv_votu_metadata_dict[votu]['lifestyle']), fontsize=12)
        
        ax.set_xlim([1e-4, 1])

        ax.axhline(y=1, c='k', ls=':', lw=1, label='Neutral')
        ax.legend(loc='lower left', fontsize=5)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        fig.subplots_adjust(hspace=0.45, wspace=0.45)
        fig_name = "%sds_vs_dnds/%s.png" % (config.analysis_directory, votu)
        fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
        plt.close()


        #total_divergence_all = numpy.asarray([dict_[g]['total_divergence'] for g in genome_pairs])
        #len_fraction_shared_blocks_all = numpy.asarray([dict_[g]['len_fraction_shared_blocks_union'] for g in genome_pairs])




def plot_ds_vs_dnds_dist_axis(votu, pseudocount=0, min_n_muts=50, min_n_sites=1e3, n_bins=30, poisson_thinning=True):


    div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
    div_dict = pickle.load(open(div_dict_path, "rb"))
    
    cumulative_n_syn_all = []
    cumulative_n_nonsyn_all = []
    cumulative_block_len_syn_all = []
    cumulative_block_len_nonsyn_all = []
    genome_pairs_clean_all = []

    for genome_pair, genome_pair_dict in div_dict.items():

        cumulative_n_syn = genome_pair_dict['cumulative_n_syn']
        cumulative_n_nonsyn = genome_pair_dict['cumulative_n_nonsyn']

        cumulative_block_len_syn = genome_pair_dict['cumulative_block_len_syn']
        cumulative_block_len_nonsyn = genome_pair_dict['cumulative_block_len_nonsyn']

        if (cumulative_n_syn == None) or (cumulative_n_nonsyn == None) or (cumulative_block_len_syn == None) or (cumulative_block_len_nonsyn == None):
            continue

        # at least one mutation in each 
        if (cumulative_n_syn == 0) or (cumulative_n_nonsyn == 0):
            continue
        
        # total of five mutations
        if (cumulative_n_syn + cumulative_n_nonsyn) <= min_n_muts:
            continue

        # at least min_n_sites possible sites
        if (cumulative_block_len_syn < min_n_sites) or (cumulative_block_len_nonsyn < min_n_sites):
            continue

        #dn = (cumulative_n_nonsyn+pseudocount)/(cumulative_block_len_nonsyn+pseudocount)
        #ds = (cumulative_n_syn+pseudocount)/(cumulative_block_len_syn+pseudocount)

        cumulative_n_syn_all.append(cumulative_n_syn)
        cumulative_n_nonsyn_all.append(cumulative_n_nonsyn)
        cumulative_block_len_syn_all.append(cumulative_block_len_syn)
        cumulative_block_len_nonsyn_all.append(cumulative_block_len_nonsyn)
        genome_pairs_clean_all.append(genome_pair)


    cumulative_n_syn_all = numpy.asarray(cumulative_n_syn_all)
    cumulative_n_nonsyn_all = numpy.asarray(cumulative_n_nonsyn_all)
    cumulative_block_len_syn_all = numpy.asarray(cumulative_block_len_syn_all)
    cumulative_block_len_nonsyn_all = numpy.asarray(cumulative_block_len_nonsyn_all)
    
    ds_all = cumulative_n_syn_all/cumulative_block_len_syn_all
    dn_all = cumulative_n_nonsyn_all/cumulative_block_len_nonsyn_all
    dnds_all = dn_all/ds_all

    if poisson_thinning == True:

        ds_1, ds_2 = data_utils.computed_poisson_thinning(cumulative_n_syn_all, cumulative_block_len_syn_all)
        
        to_plot_idx = (ds_1>0) & (ds_2>0)


    else:
        ds_1 = numpy.copy(ds_all)
        ds_2 = numpy.copy(ds_all)
        to_plot_idx = numpy.asarray([True]*len(ds_1))
    

    ds_1 = ds_1[to_plot_idx]
    ds_2 = ds_2[to_plot_idx]
    dnds_all_scatter = dn_all[to_plot_idx]/ds_2

    #ds_all_log10 = numpy.log10(ds_all)
    #dnds_all_log10 = numpy.log10(dnds_all)


    pylab.figure(figsize=(5,4))
    fig = pylab.gcf()
    outer_grid  = gridspec.GridSpec(2,2, height_ratios=[2,8], width_ratios=[8,2], hspace=0.1, wspace=0.1)

    ds_hist_axis = plt.Subplot(fig, outer_grid[0,0])
    fig.add_subplot(ds_hist_axis)
    dnds_hist_axis = plt.Subplot(fig, outer_grid[1,1])
    fig.add_subplot(dnds_hist_axis)

    scatter_axis = plt.Subplot(fig, outer_grid[1,0])
    fig.add_subplot(scatter_axis)


    same_country_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['country_code'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['country_code'] else False for k in genome_pairs_clean_all]
    same_continent_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['continent'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['continent'] else False for k in genome_pairs_clean_all]

    same_country_idx = numpy.asarray(same_country_idx)
    same_continent_idx = numpy.asarray(same_continent_idx)

    same_continent_diff_country = (same_continent_idx) & (~same_country_idx)
    same_continent_same_country = (same_continent_idx) & (same_country_idx)


    same_continent_idx_scatter = same_continent_idx[to_plot_idx]
    same_continent_diff_country_scatter = same_continent_diff_country[to_plot_idx]
    same_continent_same_country_scatter = same_continent_same_country[to_plot_idx]




    scatter_axis.scatter(ds_1[~same_continent_idx_scatter], dnds_all_scatter[~same_continent_idx_scatter], s=8, c='#FF6347', alpha=0.3, label='Diff. continent')
    scatter_axis.scatter(ds_1[same_continent_diff_country_scatter], dnds_all_scatter[same_continent_diff_country_scatter], s=8, c='#FFA500', alpha=0.3, label='Same continent, diff. country')
    scatter_axis.scatter(ds_1[same_continent_same_country_scatter], dnds_all_scatter[same_continent_same_country_scatter], s=8, c='#87CEEB', alpha=0.3, label='Same country')
    

    ds_bins_log10 = numpy.logspace(min(numpy.log10(ds_all)),max(numpy.log10(ds_all)), n_bins, base=10)
    dnds_bins_log10 = numpy.logspace(min(numpy.log10(dnds_all)), max(numpy.log10(dnds_all)), n_bins, base=10)

    ds_hist_axis.hist(ds_all[~same_continent_idx], bins=ds_bins_log10 , histtype='step', density=False, weights=numpy.ones(len(ds_all[~same_continent_idx])) / len(ds_all[~same_continent_idx]), color='#FF6347')
    ds_hist_axis.hist(ds_all[same_continent_diff_country], bins=ds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(ds_all[same_continent_diff_country])) / len(ds_all[same_continent_diff_country]), color='#FFA500')
    ds_hist_axis.hist(ds_all[same_continent_same_country], bins=ds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(ds_all[same_continent_same_country])) / len(ds_all[same_continent_same_country]), color='#87CEEB')

    dnds_hist_axis.hist(dnds_all[~same_continent_idx], bins=dnds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(dnds_all[~same_continent_idx])) / len(dnds_all[~same_continent_idx]), orientation='horizontal', color='#FF6347')
    dnds_hist_axis.hist(dnds_all[same_continent_diff_country], bins=dnds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(dnds_all[same_continent_diff_country])) / len(dnds_all[same_continent_diff_country]), orientation='horizontal', color='#FFA500')
    dnds_hist_axis.hist(dnds_all[same_continent_same_country], bins=dnds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(dnds_all[same_continent_same_country])) / len(dnds_all[same_continent_same_country]), orientation='horizontal', color='#87CEEB')

    ds_hist_axis.set_xlim([ds_bins_log10[0], ds_bins_log10[-1]])

    #d s_hist_axis.set_xlim([ds_bins[0],ds_bins[-1]])
    dnds_hist_axis.set_ylim([dnds_bins_log10[0],dnds_bins_log10[-1]])
    
    scatter_axis.set_xlim([ds_bins_log10[0], ds_bins_log10[-1]])
    scatter_axis.set_ylim([dnds_bins_log10[0], dnds_bins_log10[-1]])


    #scatter_axis_xticks = scatter_axis.get_xticks()
    #print(scatter_axis_xticks)
    scatter_axis.axhline(y=1, c='k', ls=':', lw=1, label='Neutral')
    dnds_hist_axis.axhline(y=1, c='k', ls=':', lw=1)


    scatter_axis.set_xscale('log', base=10)
    scatter_axis.set_yscale('log', base=10)

    ds_hist_axis.set_xscale('log', base=10)
    dnds_hist_axis.set_yscale('log', base=10)

    scatter_axis.set_xlabel('Synonymous divergence, ' + r'$d_{S}$', fontsize=10)
    scatter_axis.set_ylabel('Nonsynonymous ratio, ' + r'$d_{N}/d_{S}$', fontsize=10)
    ds_hist_axis.set_title('%s\nLifestyle = %s' % (votu, uhgv_votu_metadata_dict[votu]['lifestyle']), fontsize=12)




    dnds_hist_axis.set_yticklabels([])
    dnds_hist_axis.set_xticks([])
    ds_hist_axis.set_xticklabels([])
    ds_hist_axis.set_yticks([])

    scatter_axis.legend(loc='upper left',frameon=True, fontsize=6)
    fig_name = "%sds_vs_dnds_dist_axis/%s.png" % (config.analysis_directory, votu)
    fig.savefig(fig_name, bbox_inches='tight', format='png', pad_inches = 0.3, dpi = 600)
    plt.close()




def multiprocessing_calculate_divergence(votu, n_pairs):

    tic = time.time()
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=6)
    res = pool.apply_async(calculate_divergence, args=(votu, n_pairs))

    pool.close()
    pool.join()

    results = res.get()
    print(results)
    toc = time.time()
    print(f'Completed in {toc - tic} seconds')

    #processes = []
    #for _ in range(4):

    #    process = multiprocessing.Process(target=calculate_divergence, args=('vOTU-000002', ))

    #    process.start()
    ##    processes.append(process)
    #    #process.join()
    

    #for process in processes:
    #    process.join()

    #pool_size = multiprocessing.cpu_count()



def plot_divergence_with_temperate_score(min_n_muts=30, min_n_sites=1e3, rescaled_log=True):

    rescaled_log_label = {True:'_rescaled_log', False:''}

    votu_all = [filename.split('.')[0] for filename in os.listdir(config.data_directory + 'divergence_dict_all/')]
    
    from matplotlib import cm
    

    votu_to_ignore = ['vOTU-000118', 'vOTU-003790', 'vOTU-000035', 'vOTU-000016']
    votu_all = [v for v in votu_all if v not in votu_to_ignore]
    #scores = [mgv_uhgv_species_dict[v]['avg_temperate_complete'] for v in votu_all]

    if rescaled_log == True:
        ds_range = numpy.linspace(-2, 2, num=100, endpoint=True)
        dnds_range = numpy.linspace(-2, 2, num=100, endpoint=True)


    else:

        ds_range = numpy.logspace(-3, 0, num=100, endpoint=True, base=10.0)
        dnds_range = numpy.logspace(-1, 1, num=100, endpoint=True, base=10.0)


    fig = plt.figure(figsize = (12, 8))
    fig.subplots_adjust(bottom= 0.15)

    ax_ds_lifestyle = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    ax_dnds_lifestyle = plt.subplot2grid((2, 2), (1, 0), colspan=1)

    ax_ds_score = plt.subplot2grid((2, 2), (0, 1), colspan=1)
    ax_dnds_score = plt.subplot2grid((2, 2), (1, 1), colspan=1)

    for votu in votu_all:

        votu_score = float(mgv_uhgv_species_dict[votu]['avg_temperate_complete'])

        lifestyle = uhgv_votu_metadata_dict[votu]['lifestyle']
        lifestyle_color = data_utils.lifestyle_color_dict[lifestyle]
   
        if votu_score < 0:
            continue

        div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
        div_dict = pickle.load(open(div_dict_path, "rb"))
        genome_pairs = list(div_dict.keys())

        ds_all = []
        dnds_all = []
        
        for genome_pair_idx, genome_pair in enumerate(genome_pairs):

            cumulative_n_syn = div_dict[genome_pair]['cumulative_n_syn']
            cumulative_n_nonsyn = div_dict[genome_pair]['cumulative_n_nonsyn']

            cumulative_block_len_syn = div_dict[genome_pair]['cumulative_block_len_syn']
            cumulative_block_len_nonsyn = div_dict[genome_pair]['cumulative_block_len_nonsyn']

            if (cumulative_n_syn == None) or (cumulative_n_nonsyn == None) or (cumulative_block_len_syn == None) or (cumulative_block_len_nonsyn == None):
                continue

            # at least one mutation in each 
            if (cumulative_n_syn == 0) or (cumulative_n_nonsyn == 0):
                continue
            
            # total of five mutations
            if (cumulative_n_syn + cumulative_n_nonsyn) <= min_n_muts:
                continue

            # at least min_n_sites possible sites
            if (cumulative_block_len_syn < min_n_sites) or (cumulative_block_len_nonsyn < min_n_sites):
                continue

            ds = cumulative_n_syn/cumulative_block_len_syn
            dn = cumulative_n_nonsyn/cumulative_block_len_nonsyn

            ds_all.append(ds)
            dnds_all.append(dn/ds)


        ds_all = numpy.asarray(ds_all)
        dnds_all = numpy.asarray(dnds_all)
        
        if rescaled_log == True:

            ds_all = numpy.log10(ds_all)
            dnds_all = numpy.log10(dnds_all)

            ds_all = (ds_all - numpy.mean(ds_all))/numpy.std(ds_all)
            dnds_all = (dnds_all - numpy.mean(dnds_all))/numpy.std(dnds_all)

        
        survival_ds = data_utils.make_survival_dist(ds_all, ds_range)
        survival_dnds = data_utils.make_survival_dist(dnds_all, dnds_range)

        
        rgb_ = cm.Reds(votu_score)
    
        ax_ds_lifestyle.plot(ds_range, survival_ds, lw=1, ls='-', c=lifestyle_color, alpha=0.5)
        ax_dnds_lifestyle.plot(dnds_range, survival_dnds, lw=1, ls='-', c=lifestyle_color, alpha=0.5)

        #ax_ds_lifestyle.plot(ds_range, survival_ds, lw=2, ls='-', c=data_utils.lifestyle_color_dict['temperate'], alpha=0.6)
        #ax_dnds_lifestyle.plot(dnds_range, survival_dnds, lw=2, ls='-', c=data_utils.lifestyle_color_dict['temperate'], alpha=0.6)

        ax_ds_score.plot(ds_range, survival_ds, lw=1, ls='-', c=rgb_, alpha=0.5)
        ax_dnds_score.plot(dnds_range, survival_dnds, lw=1, ls='-', c=rgb_, alpha=0.5)



    if rescaled_log == True:

        syn_div_label = 'Rescaled log of synonymous divergence, ' + r'$d_{S}$'
        nonsyn_div_label = 'Rescaled log of nonsynonymous divergence, ' + r'$d_{N}/d_{S}$'


    else:
        syn_div_label = 'Synonymous divergence, ' + r'$d_{S}$'
        nonsyn_div_label = 'Nonsynonymous divergence, ' + r'$d_{N}/d_{S}$'

        ax_ds_lifestyle.set_xscale('log', base=10)
        ax_dnds_lifestyle.set_xscale('log', base=10)
        ax_ds_score.set_xscale('log', base=10)
        ax_dnds_score.set_xscale('log', base=10)


    syn_div_survival_label = "Fraction of pairs " + r'$\geq d_{S}$'
    nonsyn_div_survival_label = "Fraction of pairs " + r'$\geq d_{N}/d_{S}$'


    ax_ds_lifestyle.set_xlabel(syn_div_label, fontsize=10)
    ax_ds_lifestyle.set_ylabel(syn_div_survival_label, fontsize=10)

    ax_dnds_lifestyle.set_xlabel(nonsyn_div_label, fontsize=10)
    ax_dnds_lifestyle.set_ylabel(nonsyn_div_survival_label, fontsize=10)

    ax_ds_score.set_xlabel(syn_div_label, fontsize=10)
    ax_ds_score.set_ylabel(syn_div_survival_label, fontsize=10)

    ax_dnds_score.set_xlabel(nonsyn_div_label, fontsize=10)
    ax_dnds_score.set_ylabel(nonsyn_div_survival_label, fontsize=10)


    ax_ds_lifestyle.set_title('Lifestyle', fontsize=12)
    ax_ds_score.set_title('avg_temperate_complete score', fontsize=12)

    
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['lytic'], label='Lytic', markersize=15),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['temperate'], label='Temperate', markersize=15)]

    #ax_ds_lifestyle.legend(handles=legend_elements, loc='upper right')
    


    #cb_ax = fig.add_axes([.91,.14,.025,.70])
    #fig.colorbar(ax_ds, orientation='horizontal',cax=ax_ds)
    ##cb_ax.tick_params(labelsize=7)
    #cb_ax.set_ylabel('Mean Jaccard distance', rotation=270, labelpad=12)

    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sdivergence_with_temperate_score%s.png" % (config.analysis_directory, rescaled_log_label[rescaled_log])
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()




def plot_dnds_vs_shared_blocks(votu, min_n_muts=50, min_n_sites=1e3, n_bins=30):

    div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
    div_dict = pickle.load(open(div_dict_path, "rb"))
    
    cumulative_n_syn_all = []
    cumulative_n_nonsyn_all = []
    cumulative_block_len_syn_all = []
    cumulative_block_len_nonsyn_all = []
    genome_pairs_clean_all = []

    len_fraction_shared_blocks_all = []

    for genome_pair, genome_pair_dict in div_dict.items():

        cumulative_n_syn = genome_pair_dict['cumulative_n_syn']
        cumulative_n_nonsyn = genome_pair_dict['cumulative_n_nonsyn']

        cumulative_block_len_syn = genome_pair_dict['cumulative_block_len_syn']
        cumulative_block_len_nonsyn = genome_pair_dict['cumulative_block_len_nonsyn']

        if (cumulative_n_syn == None) or (cumulative_n_nonsyn == None) or (cumulative_block_len_syn == None) or (cumulative_block_len_nonsyn == None):
            continue

        # at least one mutation in each 
        if (cumulative_n_syn == 0) or (cumulative_n_nonsyn == 0):
            continue
        
        # total of five mutations
        if (cumulative_n_syn + cumulative_n_nonsyn) <= min_n_muts:
            continue

        # at least min_n_sites possible sites
        if (cumulative_block_len_syn < min_n_sites) or (cumulative_block_len_nonsyn < min_n_sites):
            continue

        cumulative_n_syn_all.append(cumulative_n_syn)
        cumulative_n_nonsyn_all.append(cumulative_n_nonsyn)
        cumulative_block_len_syn_all.append(cumulative_block_len_syn)
        cumulative_block_len_nonsyn_all.append(cumulative_block_len_nonsyn)
        genome_pairs_clean_all.append(genome_pair)

        len_fraction_shared_blocks_all.append(div_dict[genome_pair]['len_fraction_shared_blocks_union'])



    cumulative_n_syn_all = numpy.asarray(cumulative_n_syn_all)
    cumulative_n_nonsyn_all = numpy.asarray(cumulative_n_nonsyn_all)
    cumulative_block_len_syn_all = numpy.asarray(cumulative_block_len_syn_all)
    cumulative_block_len_nonsyn_all = numpy.asarray(cumulative_block_len_nonsyn_all)
    
    ds_all = cumulative_n_syn_all/cumulative_block_len_syn_all
    dn_all = cumulative_n_nonsyn_all/cumulative_block_len_nonsyn_all
    dnds_all = dn_all/ds_all


    len_fraction_shared_blocks_all = numpy.asarray(len_fraction_shared_blocks_all)


    #ds_all_log10 = numpy.log10(ds_all)
    #dnds_all_log10 = numpy.log10(dnds_all)


    pylab.figure(figsize=(5,4))
    fig = pylab.gcf()
    outer_grid  = gridspec.GridSpec(2,2, height_ratios=[2,8], width_ratios=[8,2], hspace=0.1, wspace=0.1)

    dnds_hist_axis = plt.Subplot(fig, outer_grid[0,0])
    fig.add_subplot(dnds_hist_axis)
    length_hist_axis = plt.Subplot(fig, outer_grid[1,1])
    fig.add_subplot(length_hist_axis)

    scatter_axis = plt.Subplot(fig, outer_grid[1,0])
    fig.add_subplot(scatter_axis)


    same_country_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['country_code'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['country_code'] else False for k in genome_pairs_clean_all]
    same_continent_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['continent'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['continent'] else False for k in genome_pairs_clean_all]

    same_country_idx = numpy.asarray(same_country_idx)
    same_continent_idx = numpy.asarray(same_continent_idx)

    same_continent_diff_country = (same_continent_idx) & (~same_country_idx)
    same_continent_same_country = (same_continent_idx) & (same_country_idx)


    to_plot_idx = numpy.asarray([True]*len(dnds_all))

    same_continent_idx_scatter = same_continent_idx[to_plot_idx]
    same_continent_diff_country_scatter = same_continent_diff_country[to_plot_idx]
    same_continent_same_country_scatter = same_continent_same_country[to_plot_idx]


    scatter_axis.scatter(dnds_all[~same_continent_idx_scatter], len_fraction_shared_blocks_all[~same_continent_idx_scatter], s=8, c='#FF6347', alpha=0.3, label='Diff. continent')
    scatter_axis.scatter(dnds_all[same_continent_diff_country_scatter], len_fraction_shared_blocks_all[same_continent_diff_country_scatter], s=8, c='#FFA500', alpha=0.3, label='Same continent, diff. country')
    scatter_axis.scatter(dnds_all[same_continent_same_country_scatter], len_fraction_shared_blocks_all[same_continent_same_country_scatter], s=8, c='#87CEEB', alpha=0.3, label='Same country')
    

    dnds_bins_log10 = numpy.logspace(min(numpy.log10(dnds_all)), max(numpy.log10(dnds_all)), n_bins, base=10)
    length_bins = numpy.linspace(min(len_fraction_shared_blocks_all), 1, n_bins)

    dnds_hist_axis.hist(dnds_all[~same_continent_idx], bins=dnds_bins_log10 , histtype='step', density=False, weights=numpy.ones(len(ds_all[~same_continent_idx])) / len(ds_all[~same_continent_idx]), color='#FF6347')
    dnds_hist_axis.hist(dnds_all[same_continent_diff_country], bins=dnds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(ds_all[same_continent_diff_country])) / len(ds_all[same_continent_diff_country]), color='#FFA500')
    dnds_hist_axis.hist(dnds_all[same_continent_same_country], bins=dnds_bins_log10, histtype='step', density=False, weights=numpy.ones(len(ds_all[same_continent_same_country])) / len(ds_all[same_continent_same_country]), color='#87CEEB')

    length_hist_axis.hist(len_fraction_shared_blocks_all[~same_continent_idx], bins=length_bins, histtype='step', density=False, weights=numpy.ones(len(dnds_all[~same_continent_idx])) / len(dnds_all[~same_continent_idx]), orientation='horizontal', color='#FF6347')
    length_hist_axis.hist(len_fraction_shared_blocks_all[same_continent_diff_country], bins=length_bins, histtype='step', density=False, weights=numpy.ones(len(dnds_all[same_continent_diff_country])) / len(dnds_all[same_continent_diff_country]), orientation='horizontal', color='#FFA500')
    length_hist_axis.hist(len_fraction_shared_blocks_all[same_continent_same_country], bins=length_bins, histtype='step', density=False, weights=numpy.ones(len(dnds_all[same_continent_same_country])) / len(dnds_all[same_continent_same_country]), orientation='horizontal', color='#87CEEB')

    dnds_hist_axis.set_xlim([dnds_bins_log10[0], dnds_bins_log10[-1]])

    #d s_hist_axis.set_xlim([ds_bins[0],ds_bins[-1]])
    length_hist_axis.set_ylim([length_bins[0],length_bins[-1]])
    
    scatter_axis.set_xlim([dnds_bins_log10[0], dnds_bins_log10[-1]])
    scatter_axis.set_ylim([length_bins[0], length_bins[-1]])


    #scatter_axis_xticks = scatter_axis.get_xticks()
    #print(scatter_axis_xticks)
    scatter_axis.axvline(x=1, c='k', ls=':', lw=1, label='Neutral')
    dnds_hist_axis.axvline(x=1, c='k', ls=':', lw=1)


    scatter_axis.set_xscale('log', base=10)
    #scatter_axis.set_yscale('log', base=10)

    dnds_hist_axis.set_xscale('log', base=10)
    #length_hist_axis.set_yscale('log', base=10)
    
    scatter_axis.set_xlabel('Nonsynonymous ratio, ' + r'$d_{N}/d_{S}$', fontsize=10)
    scatter_axis.set_ylabel('Length fraction of shared blocks', fontsize=10)
    dnds_hist_axis.set_title('%s\nLifestyle = %s' % (votu, uhgv_votu_metadata_dict[votu]['lifestyle']), fontsize=12)


    length_hist_axis.set_yticklabels([])
    length_hist_axis.set_xticks([])
    dnds_hist_axis.set_xticklabels([])
    dnds_hist_axis.set_yticks([])

    scatter_axis.legend(loc='upper left',frameon=True, fontsize=6)
    fig_name = "%sdnds_vs_shared_blocks/%s.png" % (config.analysis_directory, votu)
    fig.savefig(fig_name, bbox_inches='tight', format='png', pad_inches = 0.3, dpi = 600)
    plt.close()



    
def plot_divergence_along_genome(votu, genome_1, genome_2, min_n_muts=300, min_n_sites=2e3):

    def find_possible_pairs():

        div_dict_path = '%sdivergence_dict_all/%s.pickle' % (config.data_directory, votu)
        div_dict = pickle.load(open(div_dict_path, "rb"))

        cumulative_n_syn_all = []
        cumulative_n_nonsyn_all = []
        cumulative_block_len_syn_all = []
        cumulative_block_len_nonsyn_all = []
        genome_pairs_clean_all = []

        for genome_pair, genome_pair_dict in div_dict.items():

            cumulative_n_syn = genome_pair_dict['cumulative_n_syn']
            cumulative_n_nonsyn = genome_pair_dict['cumulative_n_nonsyn']

            cumulative_block_len_syn = genome_pair_dict['cumulative_block_len_syn']
            cumulative_block_len_nonsyn = genome_pair_dict['cumulative_block_len_nonsyn']

            if (cumulative_n_syn == None) or (cumulative_n_nonsyn == None) or (cumulative_block_len_syn == None) or (cumulative_block_len_nonsyn == None):
                    continue

            # at least one mutation in each 
            if (cumulative_n_syn == 0) or (cumulative_n_nonsyn == 0):
                continue
            
            # total of five mutations
            if (cumulative_n_syn + cumulative_n_nonsyn) <= min_n_muts:
                continue

            # at least min_n_sites possible sites
            if (cumulative_block_len_syn < min_n_sites) or (cumulative_block_len_nonsyn < min_n_sites):
                continue

            cumulative_n_syn_all.append(cumulative_n_syn)
            cumulative_n_nonsyn_all.append(cumulative_n_nonsyn)
            cumulative_block_len_syn_all.append(cumulative_block_len_syn)
            cumulative_block_len_nonsyn_all.append(cumulative_block_len_nonsyn)
            genome_pairs_clean_all.append(genome_pair)


        cumulative_n_syn_all = numpy.asarray(cumulative_n_syn_all)
        cumulative_n_nonsyn_all = numpy.asarray(cumulative_n_nonsyn_all)
        cumulative_block_len_syn_all = numpy.asarray(cumulative_block_len_syn_all)
        cumulative_block_len_nonsyn_all = numpy.asarray(cumulative_block_len_nonsyn_all)
        
        ds_all = cumulative_n_syn_all/cumulative_block_len_syn_all
        dn_all = cumulative_n_nonsyn_all/cumulative_block_len_nonsyn_all
        dnds_all = dn_all/ds_all

        same_country_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['country_code'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['country_code'] else False for k in genome_pairs_clean_all]
        same_continent_idx = [True if sample_metagenome_dict[vgenome_dict[k[0]]['original_id']]['continent'] == sample_metagenome_dict[vgenome_dict[k[1]]['original_id']]['continent'] else False for k in genome_pairs_clean_all]

        same_country_idx = numpy.asarray(same_country_idx)
        same_continent_idx = numpy.asarray(same_continent_idx)

        same_continent_diff_country = (same_continent_idx) & (~same_country_idx)
        #same_continent_same_country = (same_continent_idx) & (same_country_idx)
        
        boundary_ds = 0.02
        boundary_dnds = 1
        left_block = (ds_all < boundary_ds) & (dnds_all > boundary_dnds)
        right_block = (ds_all > boundary_ds) & (dnds_all < boundary_dnds)

        left_possible_pairs = numpy.where((same_continent_diff_country & left_block) == True)[0]
        right_possible_pairs = numpy.where((same_continent_diff_country & right_block) == True)[0]


        left_pair = genome_pairs_clean_all[left_possible_pairs[10]]
        right_pair = genome_pairs_clean_all[right_possible_pairs[10]]

        return left_pair, right_pair


    left_pair, right_pair = find_possible_pairs()


    # checking if it is a file
    pangraph_data = data_utils.load_pangraph_data(minimap_path % votu)
    syn_sites_dict = pickle.load(open(syn_sites_path_template % votu, "rb"))

    mut_position_syn_left, mut_position_nonsyn_left, block_position_syn_left, block_position_nonsyn_left, block_inter_id_left, block_position_dict_left = data_utils.calculate_annotated_divergence_across_pangraph_blocks(left_pair[0], left_pair[1], pangraph_data, syn_sites_dict=syn_sites_dict)
    mut_position_syn_right, mut_position_nonsyn_right, block_position_syn_right, block_position_nonsyn_right, block_inter_id_right, block_position_dict_right = data_utils.calculate_annotated_divergence_across_pangraph_blocks(right_pair[0], right_pair[1], pangraph_data, syn_sites_dict=syn_sites_dict)

    # blocks shared by both comparisons
    inter_block = numpy.intersect1d(block_inter_id_left, block_inter_id_right)
    
    # sort blocks by position
    #start_block = [(b, block_inter_id_left[b][0]) for b in inter_block]
    #start_block = list(sorted(start_block, key=operator.itemgetter(1)))

    #n_sites_per_window = 500
    #n_
    # get sites for all 
    #inter_block_all_sites = []
    #for s in inter_block:
    #    inter_block_all_sites.extend(list(range( block_position_dict_left[s][0],  block_position_dict_left[s][1])))
    #inter_block_all_sites.sort()
    #inter_block_all_sites = numpy.asarray(inter_block_all_sites)
    
    n_sites_per_window = 1000
    #n_windows = int(len(inter_block_all_sites)/n_sites_per_window)

    #min_len_left = min([min(block_position_syn_left), min(block_position_nonsyn_left)])
    #max_len_left = max([max(block_position_syn_left), max(block_position_nonsyn_left)])

    all_bins_left = list(range(min([min(block_position_syn_left), min(block_position_nonsyn_left)]), max([max(block_position_syn_left), max(block_position_nonsyn_left)])+1, n_sites_per_window))
    all_bins_right = list(range(min([min(block_position_syn_right), min(block_position_nonsyn_right)]), max([max(block_position_syn_right), max(block_position_nonsyn_right)])+1, n_sites_per_window))

    ds_left = []
    dn_left = []
    for i in range(len(all_bins_left)-1):

        n_muts_syn_i = sum((mut_position_syn_left >= all_bins_left[i]) & (mut_position_syn_left < all_bins_left[i+1]))
        n_sites_syn_i = sum((block_position_syn_left >= all_bins_left[i]) & (block_position_syn_left < all_bins_left[i+1]))

        n_muts_nonsyn_i = sum((mut_position_nonsyn_left >= all_bins_left[i]) & (mut_position_nonsyn_left < all_bins_left[i+1]))
        n_sites_nonsyn_i = sum((block_position_nonsyn_left >= all_bins_left[i]) & (block_position_nonsyn_left < all_bins_left[i+1]))

        ds_left.append(n_muts_syn_i/n_sites_syn_i)
        dn_left.append(n_muts_nonsyn_i/n_sites_nonsyn_i)


    ds_right = []
    dn_right = []
    for i in range(len(all_bins_left)-1):

        n_muts_syn_i = sum((mut_position_syn_right >= all_bins_right[i]) & (mut_position_syn_right < all_bins_right[i+1]))
        n_sites_syn_i = sum((block_position_syn_right >= all_bins_right[i]) & (block_position_syn_right < all_bins_right[i+1]))

        n_muts_nonsyn_i = sum((mut_position_nonsyn_right >= all_bins_right[i]) & (mut_position_nonsyn_right < all_bins_right[i+1]))
        n_sites_nonsyn_i = sum((block_position_nonsyn_right >= all_bins_right[i]) & (block_position_nonsyn_right < all_bins_right[i+1]))

        ds_right.append(n_muts_syn_i/n_sites_syn_i)
        dn_right.append(n_muts_nonsyn_i/n_sites_nonsyn_i)

    #ds_left = [sum((mut_position_syn_left >= all_bins_left[i]) & (mut_position_syn_left < all_bins_left[i+1])) for i in range(len(all_bins_left)-1)]

    all_bins_left = numpy.asarray(all_bins_left)
    all_bins_right = numpy.asarray(all_bins_right)

    all_bins_left_midpoint = all_bins_left*0.5
    all_bins_right_midpoint = all_bins_right*0.5


    ds_left = numpy.asarray(ds_left)
    dn_left = numpy.asarray(dn_left)

    ds_right = numpy.asarray(ds_right)
    dn_right = numpy.asarray(dn_right)

    dnds_left = dn_left/ds_left
    dnds_right = dn_right/ds_right


    fig = plt.figure(figsize = (6, 8))
    fig.subplots_adjust(bottom= 0.15)

    ax_ds = plt.subplot2grid((2, 1), (0, 0), colspan=1)
    ax_dn = plt.subplot2grid((2, 1), (1, 0), colspan=1)

    ds_left_idx = (~numpy.isnan(ds_left)) & (ds_left>0)
    ds_right_idx = (~numpy.isnan(ds_right)) & (ds_right>0)

    dnds_left_idx = (~numpy.isnan(dnds_left)) & (dnds_left>0)
    dnds_right_idx = (~numpy.isnan(dnds_right)) & (dnds_right>0)

    # plot shaded regions for shared blocks

    for i in inter_block:
        ax_ds.axvspan(block_position_dict_left[i][0], block_position_dict_left[i][1], facecolor='mediumpurple', alpha=0.3)
        ax_dn.axvspan(block_position_dict_left[i][0], block_position_dict_left[i][1], facecolor='mediumpurple', alpha=0.3)

    
    ax_ds.scatter(all_bins_left[:-1][ds_left_idx], ds_left[ds_left_idx], s=10, c=data_utils.lifestyle_color_dict['lytic'], zorder=2)
    ax_ds.scatter(all_bins_right[:-1][ds_right_idx], ds_right[ds_right_idx], s=10, c=data_utils.lifestyle_color_dict['temperate'], zorder=2)


    ax_dn.scatter(all_bins_left[:-1][dnds_left_idx], dnds_left[dnds_left_idx], s=14, c=data_utils.lifestyle_color_dict['lytic'], zorder=2)
    ax_dn.scatter(all_bins_right[:-1][dnds_right_idx], dnds_right[dnds_right_idx], s=14, c=data_utils.lifestyle_color_dict['temperate'], zorder=2)

    
    ax_ds.set_yscale('log', base=10)
    ax_dn.set_yscale('log', base=10)

    ax_ds.set_xlabel("Cumulative pangraph block position", fontsize=10)
    ax_dn.set_xlabel("Cumulative pangraph block position", fontsize=10)
    

    ax_ds.set_ylabel('Synonymous divergence, ' + r'$d_{S}$', fontsize=10)
    ax_dn.set_ylabel('Nonsynonymous ratio, ' + r'$d_{N}/d_{S}$', fontsize=10)
 
    

    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['lytic'], label='Upper-left quadrant', markersize=12),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=data_utils.lifestyle_color_dict['temperate'], label='Lower-right quadrant', markersize=12),
                       Patch(facecolor='mediumpurple', edgecolor='mediumpurple', label='Block present in both pairs')]


    ax_ds.legend(handles=legend_elements, loc='upper left', fontsize=8)

    ax_dn.axhline(y=1, c='k', lw=1.5, ls=':')


    #scatter_axis.legend(loc='upper left',frameon=True, fontsize=6)
    fig_name = "%sdivergence_comparison.png" % config.analysis_directory
    fig.savefig(fig_name, bbox_inches='tight', format='png', pad_inches = 0.3, dpi = 600)
    plt.close()








if __name__ == "__main__":

    votu_to_ignore = ['vOTU-000118', 'vOTU-003790', 'vOTU-000035', 'vOTU-000016']
    # 300*299/2
    n_pairs=44850

    #plot_shared_blocks_dist()

    #plot_divergence_with_temperate_score()

    #plot_divergence_with_temperate_score(rescaled_log=True)

    #plot_divergence_vs_shared_blocks('vOTU-000002', min_n_muts=50, min_n_sites=1e3)

    plot_divergence_along_genome('vOTU-000002', '1', '2')



    #for filename in os.listdir(config.data_directory + 'divergence_dict_all/'):
        
    #    votu = filename.split('.')[0]

    #    if votu in votu_to_ignore:
    #        continue

    #    plot_dnds_vs_shared_blocks(votu, min_n_muts=50, min_n_sites=1e3)

        #calculate_divergence(votu)
        #plot_divergence_vs_shared_blocks(votu, min_n_muts=30, min_n_sites=1e3)
        #plot_ds_vs_dnds_dist_axis(votu, poisson_thinning=True, min_n_muts=50, min_n_sites=1e3)


    #multiprocessing_calculate_divergence('vOTU-000001', n_pairs)





    




