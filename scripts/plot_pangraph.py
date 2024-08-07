
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config
import json




pangraph_path = '%svOTU-000085-polished-pangraph.json' % config.data_directory
pangraph_file = open(pangraph_path, 'r')


pangraph_data = json.load(pangraph_file)

# blocks
### id = str
### sequence = AACCF
### positions = where on the genome the block is from
### 

pangraph_file.close()


def plot_block_size():

    #block_size_all = [len(i['sequence'])*len(pangraph_data['blocks'][0]['positions']) for i in pangraph_data['blocks']]  

    n_genomes = len(pangraph_data['paths'])
    #blocks_all = [len(i['positions'])/n_genomes for i in pangraph_data['blocks']]  

    blocks_all = []
    block_sizes_all = []
    for i in pangraph_data['blocks']:
        genomes_i = [ j[0]['name'] for j in i['positions']]
        n_genomes_i = len(set(genomes_i))

        blocks_all.append(n_genomes_i)
        block_sizes_all.extend([len(i['sequence'])] * n_genomes_i)


    # length notmalized frequency

    fig = plt.figure(figsize = (8, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax_blocks = plt.subplot2grid((1, 2), (0, 0), colspan=1)
    ax_blocks_size = plt.subplot2grid((1, 2), (0, 1), colspan=1)

    ax_blocks.hist(blocks_all, bins=20)
    ax_blocks_size.hist(block_sizes_all, bins=20)

    ax_blocks.set_xlabel('Fraction of genomes\nwhere block is present', fontsize=12)
    ax_blocks.set_ylabel('# blocks', fontsize=12)

    ax_blocks_size.set_xlabel('# bps', fontsize=12)
    ax_blocks_size.set_ylabel('# blocks', fontsize=12)


    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sblock_size_dist_vOTU-000085.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



def plot_block_diversity():

    #diversity = [ len(i['mutate'][0][1]) for i in pangraph_data['blocks'] ]

    n_genomes = len(pangraph_data['paths'])

    frac_genomes_all = []
    diversity_all = []
    for block_i in pangraph_data['blocks']:

        n_muts_i = 0
        for j in block_i['mutate']:
            n_muts_i += len(j[1])

        if n_muts_i == 0:
            continue

        genomes_i = [ j[0]['name'] for j in block_i['positions']]
        n_genomes_i = len(set(genomes_i))


        diversity_i = n_muts_i/(len(block_i['sequence']) * len(block_i['mutate']))
        #diversity_i = diversity_i/

        frac_genomes_all.append(n_genomes_i/n_genomes)
        diversity_all.append(diversity_i)


    fig = plt.figure(figsize = (4, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)

    ax.scatter(frac_genomes_all, diversity_all, alpha=0.7, s=10)

    ax.set_xlabel('Fraction of genomes where block is present', fontsize=12)
    ax.set_ylabel('Avergage diversity in block', fontsize=12)

    ax.set_yscale('log', basey=10)
    
    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sblock_diversity_vOTU-000085.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



#def plot_pairwise_blocks():


#plot_block_diversity()

#print(pangraph_data['blocks'][0]['mutate'])
