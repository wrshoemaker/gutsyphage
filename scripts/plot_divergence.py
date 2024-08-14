import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle
import time
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

import multiprocessing



random.seed(123456789)


#divergence_dict_path = '%sdivergence_dict.pickle' % config.data_directory

votu_dict, vgenome_dict = data_utils.read_uhgv_metadata()
uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()
sample_metagenome_dict = data_utils.read_sample_metadata()


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





if __name__ == "__main__":

    votu_to_ignore = ['vOTU-000118', 'vOTU-003790', 'vOTU-000035', 'vOTU-000016']
    # 300*299/2
    n_pairs=44850

    for filename in os.listdir(config.data_directory + 'divergence_dict_all/'):
        
        votu = filename.split('.')[0]

        if votu in votu_to_ignore:
            continue

        #calculate_divergence(votu)

        plot_divergence_vs_shared_blocks(votu, min_n_muts=30, min_n_sites=1e3)
        plot_ds_vs_dnds_dist_axis(votu, poisson_thinning=True, min_n_muts=50, min_n_sites=1e3)


    #multiprocessing_calculate_divergence('vOTU-000001', n_pairs)

    
    #calculate_divergence('vOTU-000001', n_pairs)
    #calculate_divergence('vOTU-000132')
    #calculate_divergence('vOTU-000005')
    #calculate_divergence('vOTU-000039')
    #calculate_divergence('vOTU-000008')

    #plot_ds_vs_dnds_dist_axis('vOTU-000002', poisson_thinning=True)

    #plot_divergence_vs_shared_blocks('vOTU-000002')




    

    

#plot_ds_vs_dnds_dist_axis('vOTU-000007', poisson_thinning=True)
#plot_ds_vs_dnds_dist_axis('vOTU-000003', poisson_thinning=True)




    

#plot_divergence_vs_shared_blocks()

#plot_shared_blocks_dist()

\



