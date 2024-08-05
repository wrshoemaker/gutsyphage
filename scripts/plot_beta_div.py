
import pickle
import sys

import numpy
import scipy.stats as stats
import scipy.spatial as spatial

import collections
import data_utils
import config

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
#from matplotlib import colormaps, cm

from itertools import chain, combinations


min_n_country_samples = 50

dist_dict_path = '%sdist_dict.pickle' % config.data_directory


# load vOTU metadata with taxonomy
uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()

# get phage presence/absence dictionary
phage_metadata_dict = data_utils.read_uhgv_metadata(checkv_quality='Complete', checkv_quality_cumulative=True)
# build the presence/absence matrix
votu_pres_abs_array, votu_names, sample_names = data_utils.votu_dict_to_pres_abs(phage_metadata_dict)
# get country IDs from mgv_sample_info.tsv.gz using sample_accession
sample_metagenome_dict = data_utils.read_sample_metadata()





def plot_n_samples_per_country():

    run_accessions_to_country_code_dict = dict( [ (sample_metagenome_dict[k]['sample_accession'], sample_metagenome_dict[k]['country_code']) for k in sample_metagenome_dict.keys()])

    country_code_all = [run_accessions_to_country_code_dict[s] for s in sample_names]

    country_code_counts_dict = collections.Counter(country_code_all)
    country_code_counts_dict = country_code_counts_dict.most_common()

    x_axis = [x[1] for x in country_code_counts_dict]
    y_axis_labels = [data_utils.country_code_dict[x[0]] for x in country_code_counts_dict]

    fig = plt.figure(figsize = (4, 4))
    fig.subplots_adjust(bottom= 0.15)

    ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)


    y_range = list(range(len(country_code_counts_dict)))
    ax.barh(y_axis_labels, x_axis, color='dodgerblue')
    #ax.set_yticklabels(x, fontsize=20)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_xlabel('# hosts', fontsize=12)

    ax.set_title('vOTU analysis sample size', fontsize=12)


    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%scountry_sample_size.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()




def make_votu_dist_dict(min_n_country_samples=min_n_country_samples):

    run_accessions_to_country_code_dict = dict( [ (sample_metagenome_dict[k]['sample_accession'], sample_metagenome_dict[k]['country_code']) for k in sample_metagenome_dict.keys()])

    country_code_all = [run_accessions_to_country_code_dict[s] for s in sample_names]
    country_code_counts_dict = collections.Counter(country_code_all)
    country_code_all_array = numpy.asarray(country_code_all)

    # only look at countries with a sufficient # hosts
    country_code_to_keep = [k for k,v in country_code_counts_dict.items() if v >= min_n_country_samples]
    
    #country_code_to_keep = numpy.asarray(country_code_to_keep)
    country_code_set_sorted = sorted(list(set(country_code_to_keep)))
    country_code_pairs = list(combinations(country_code_set_sorted, 2))
    
    # loop through
    # lifestyles
    # within/between

    lifestyle_votu_all = numpy.asarray([uhgv_votu_metadata_dict[v]['lifestyle'] for v in votu_names])
    lifestyle_all = data_utils.lifestyle_all
    lifestyle_all.append('both')

    dist_dict = {}

    for lifestyle in lifestyle_all:

        if lifestyle != 'both':

            lifestyle_idx = (lifestyle_votu_all == lifestyle)
            votu_names_lifestyle = votu_names[lifestyle_idx]
            votu_pres_abs_array_lifestyle = votu_pres_abs_array[lifestyle_idx,:]

        else:
            votu_names_lifestyle = numpy.copy(votu_names)
            votu_pres_abs_array_lifestyle = numpy.copy(votu_pres_abs_array)


        dist_dict[lifestyle] = {}
        
        for taxonomy_type in ['ictv_taxonomy', 'uhgv_taxonomy']:

            print(lifestyle, taxonomy_type)

            taxonomy_ranks = data_utils.phage_taxonomy_ranks_dict[taxonomy_type]
            dist_dict[lifestyle][taxonomy_type] = {}

            for t in taxonomy_ranks:
            
                votu_pres_abs_cg_array, cg_taxa_labels, cg_n_collapsed_nodes, cg_taxa_to_keep_idx = data_utils.coarse_grain_abundances_by_taxonomy(votu_pres_abs_array_lifestyle, votu_names_lifestyle, uhgv_votu_metadata_dict, taxonomic_level=t, phage_taxonomic_system=taxonomy_type, pres_abs=True)
                dist_dict[lifestyle][taxonomy_type][t] = {}

                for c_pair in country_code_pairs:

                    country_code_all_array_c_pair_idx = (country_code_all_array == c_pair[0]) | (country_code_all_array == c_pair[1])

                    votu_pres_abs_cg_array_c_pair = votu_pres_abs_cg_array[:,country_code_all_array_c_pair_idx]
                    dist_c = spatial.distance.pdist(votu_pres_abs_cg_array_c_pair.T, metric='jaccard')

                    # identify indices from different samples
                    country_code_all_array_c_pair = country_code_all_array[country_code_all_array_c_pair_idx]
                    country_code_all_array_c_pair[country_code_all_array_c_pair == c_pair[0]] = 0
                    country_code_all_array_c_pair[country_code_all_array_c_pair == c_pair[1]] = 1

                    # get indices where the two samples are different 
                    bw_idx = spatial.distance.pdist(country_code_all_array_c_pair.reshape(country_code_all_array_c_pair.shape[0],-1), metric='jaccard') == 1
                    dist_c_bw = dist_c[bw_idx]
                    
                    # mean distance between.
                    mean_dist_c_bw = numpy.mean(dist_c_bw)

                    lower_quantile = numpy.quantile(dist_c_bw, q=0.25)
                    upper_quantile = numpy.quantile(dist_c_bw, q=0.75)

                    dist_dict[lifestyle][taxonomy_type][t][c_pair] = {}
                    dist_dict[lifestyle][taxonomy_type][t][c_pair]['mean_jaccard'] = mean_dist_c_bw
                    dist_dict[lifestyle][taxonomy_type][t][c_pair]['lower_quantile'] = lower_quantile
                    dist_dict[lifestyle][taxonomy_type][t][c_pair]['upper_quantile'] = upper_quantile

                    

                for c in country_code_set_sorted:

                    votu_pres_abs_cg_array_c = votu_pres_abs_cg_array[:,(country_code_all_array == c)]
                    dist_c = spatial.distance.pdist(votu_pres_abs_cg_array_c.T, metric='jaccard')

                    # mean distance within
                    mean_dist_c_win = numpy.mean(dist_c)

                    lower_quantile = numpy.quantile(dist_c, q=0.25)
                    upper_quantile = numpy.quantile(dist_c, q=0.75)

                    dist_dict[lifestyle][taxonomy_type][t][c] = {}
                    dist_dict[lifestyle][taxonomy_type][t][c]['mean_jaccard'] = mean_dist_c_win
                    dist_dict[lifestyle][taxonomy_type][t][c]['lower_quantile'] = lower_quantile
                    dist_dict[lifestyle][taxonomy_type][t][c]['upper_quantile'] = upper_quantile




    sys.stderr.write("Saving dictionary...\n")
    with open(dist_dict_path, 'wb') as handle:
        pickle.dump(dist_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")





def plot_dist_heatmap():

    run_accessions_to_country_code_dict = dict( [ (sample_metagenome_dict[k]['sample_accession'], sample_metagenome_dict[k]['country_code']) for k in sample_metagenome_dict.keys()])

    country_code_all = [run_accessions_to_country_code_dict[s] for s in sample_names]
    country_code_counts_dict = collections.Counter(country_code_all)
    country_code_all_array = numpy.asarray(country_code_all)

    # only look at countries with a sufficient # hosts
    country_code_to_keep = [k for k,v in country_code_counts_dict.items() if v >= min_n_country_samples]
    
    #country_code_to_keep = numpy.asarray(country_code_to_keep)
    country_code_set_sorted = sorted(list(set(country_code_to_keep)))
    country_code_pairs = list(combinations(country_code_set_sorted, 2))

    country_code_set_sorted_label = [data_utils.country_code_dict[c] for c in country_code_set_sorted]

    # load dictionary
    dist_dict = pickle.load(open(dist_dict_path, "rb"))

    ranks_to_plot = ['otu', 'genus', 'subfamily', 'family'] 

    fig = plt.figure(figsize = (16, 12))
    fig.subplots_adjust(bottom= 0.15)

    #cmap = colormaps.get_cmap('Blues')
    #cmap.set_under('k')

    for lifestyle_idx, lifestyle in enumerate(data_utils.lifestyle_all_and_both[::-1]):

        for rank_idx, rank in enumerate(ranks_to_plot):

            dist = numpy.asarray([dist_dict[lifestyle]['uhgv_taxonomy'][rank][c_pair]['mean_jaccard'] for c_pair in country_code_pairs])
            dist = spatial.distance.squareform(dist)

            dist_diagonal = numpy.asarray([dist_dict[lifestyle]['uhgv_taxonomy'][rank][c]['mean_jaccard'] for c in country_code_set_sorted])

            ax = plt.subplot2grid((len(data_utils.lifestyle_all_and_both), len(ranks_to_plot)), (lifestyle_idx, rank_idx))

            # plot heatmap

            dist_tril = numpy.tril(dist)
            dist_tril[numpy.triu_indices(dist_tril.shape[0], 1)] = numpy.nan

            numpy.fill_diagonal(dist_tril, dist_diagonal)
    
            ax_im = ax.imshow(dist_tril, interpolation="nearest", cmap='Blues',  vmin=0.65, vmax=1)

            dummy_range = range(dist.shape[0])
            #pcm = ax.pcolor(dummy_range, dummy_range, dist, cmap='GnBu', norm=colors.Normalize(vmin=0, vmax=1))
            #pcm.cmap.set_under('black')   


            ax.set_xticks(dummy_range)
            ax.set_xticklabels(country_code_set_sorted_label, fontsize=6, rotation=45, ha="right")

            ax.set_yticks(dummy_range)
            ax.set_yticklabels(country_code_set_sorted_label, fontsize=6, rotation=45)


            if lifestyle_idx == 0:
                ax.set_title('UHGV %s scale' % data_utils.uhgv_taxonomy_name_no_cap_dict[rank], fontsize=14)

            if rank_idx == 0:
                ax.text(-0.25, 0.5, lifestyle.capitalize(), fontsize=14, ha='center', va='center',rotation=90, transform=ax.transAxes)



            if (lifestyle_idx==0) and (rank_idx==2):
                cb_ax = fig.add_axes([.91,.14,.025,.70])
                fig.colorbar(ax_im, orientation='vertical',cax=cb_ax)
                cb_ax.tick_params(labelsize=7)
                cb_ax.set_ylabel('Mean Jaccard distance', rotation=270, labelpad=12)



    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sjaccard_heatmap.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()



plot_dist_heatmap()


#make_votu_dist_dict()



# get original_id for all genomes and flatten into list
#original_id_all = list(set(chain.from_iterable([phage_metadata_dict[v]['original_id'] for v in phage_metadata_dict.keys()])))


#country_code_dict = dict( [ (i, sample_metagenome_dict[i]['country_code']) for i in original_id_all])





#plot_n_samples_per_country()

#make_votu_dist_dict()

#country_label = [sample_metagenome_dict[v]['country_code'] for v in votu_names]


