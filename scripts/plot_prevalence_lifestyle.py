import pickle
import sys

import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.spatial as spatial

import collections
import data_utils
import config

from itertools import chain, combinations



prevalence_dict_path = '%sprevalence_dict.pickle' % config.data_directory



def make_prevalence_dict():

    prevalence_dict = {}

    # load vOTU metadata with taxonomy
    uhgv_votu_metadata_dict = data_utils.read_uhgv_votu_metadata()

    # get phage presence/absence dictionary
    phage_metadata_dict = data_utils.read_uhgv_metadata(checkv_quality='Complete', checkv_quality_cumulative=True)
    # build the presence/absence matrix
    votu_pres_abs_array, votu_names, sample_names = data_utils.votu_dict_to_pres_abs(phage_metadata_dict)
    # get country IDs from mgv_sample_info.tsv.gz using sample_accession
    sample_metagenome_dict = data_utils.read_sample_metadata()

    lifestyle_votu_all = numpy.asarray([uhgv_votu_metadata_dict[v]['lifestyle'] for v in votu_names])

    lifestyle_all = data_utils.lifestyle_all
    lifestyle_all.append('both')


    for lifestyle in data_utils.lifestyle_all:

        if lifestyle != 'both':

            lifestyle_idx = (lifestyle_votu_all == lifestyle)
            votu_names_lifestyle = votu_names[lifestyle_idx]
            votu_pres_abs_array_lifestyle = votu_pres_abs_array[lifestyle_idx,:]

        else:
            votu_names_lifestyle = numpy.copy(votu_names)
            votu_pres_abs_array_lifestyle = numpy.copy(votu_pres_abs_array)


        #lifestyle_idx = (lifestyle_votu_all == lifestyle)
        #votu_names_lifestyle = votu_names[lifestyle_idx]
        #votu_pres_abs_array_lifestyle = votu_pres_abs_array[lifestyle_idx,:]

        # remove columns with all zeros
        #votu_pres_abs_array_lifestyle = votu_pres_abs_array_lifestyle[:,~numpy.all(votu_pres_abs_array_lifestyle==0, axis=0)]
        
        prevalence_lifestyle = numpy.sum(votu_pres_abs_array_lifestyle, axis=1)/votu_pres_abs_array_lifestyle.shape[1]

        prevalence_dict[lifestyle] = {}
        prevalence_dict[lifestyle]['votu'] = prevalence_lifestyle.tolist()
        prevalence_dict[lifestyle]['ictv_taxonomy'] = {}
        prevalence_dict[lifestyle]['uhgv_taxonomy'] = {}

        for ictv_taxonomy_ranks_i in data_utils.ictv_taxonomy_ranks:

            votu_pres_abs_cg_array, cg_taxa_labels, cg_n_collapsed_nodes, cg_taxa_to_keep_idx = data_utils.coarse_grain_abundances_by_taxonomy(votu_pres_abs_array_lifestyle, votu_names_lifestyle, uhgv_votu_metadata_dict, taxonomic_level=ictv_taxonomy_ranks_i, phage_taxonomic_system='ictv_taxonomy', pres_abs=True)
            prevalence_lifestyle_i = numpy.sum(votu_pres_abs_cg_array, axis=1)/votu_pres_abs_cg_array.shape[1]
            prevalence_dict[lifestyle]['ictv_taxonomy'][ictv_taxonomy_ranks_i] = prevalence_lifestyle_i.tolist()
            
            print(lifestyle, 'ICTV', ictv_taxonomy_ranks_i, len(prevalence_lifestyle_i))



        for uhgv_taxonomy_ranks_i in data_utils.uhgv_taxonomy_ranks:

            # coarse-grain
            votu_pres_abs_cg_array, cg_taxa_labels, cg_n_collapsed_nodes, cg_taxa_to_keep_idx = data_utils.coarse_grain_abundances_by_taxonomy(votu_pres_abs_array_lifestyle, votu_names_lifestyle, uhgv_votu_metadata_dict, taxonomic_level=uhgv_taxonomy_ranks_i, phage_taxonomic_system='uhgv_taxonomy', pres_abs=True)
            prevalence_lifestyle_i = numpy.sum(votu_pres_abs_cg_array, axis=1)/votu_pres_abs_cg_array.shape[1]
            prevalence_dict[lifestyle]['uhgv_taxonomy'][uhgv_taxonomy_ranks_i] = prevalence_lifestyle_i.tolist()

            print(lifestyle, 'UHGV', uhgv_taxonomy_ranks_i, len(prevalence_lifestyle_i))

    

    sys.stderr.write("Saving dictionary...\n")
    with open(prevalence_dict_path, 'wb') as handle:
        pickle.dump(prevalence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stderr.write("Done!\n")


def plot_prevalence():


    fig = plt.figure(figsize = (8, 8))
    fig.subplots_adjust(bottom= 0.15)

    # get phage presence/absence dictionary
    phage_metadata_dict = data_utils.read_uhgv_metadata(checkv_quality='Complete', checkv_quality_cumulative=True)
    # build the presence/absence matrix
    votu_pres_abs_array, votu_names, sample_names = data_utils.votu_dict_to_pres_abs(phage_metadata_dict)

    prevalence_range = numpy.logspace(-1*numpy.log10(len(sample_names)), 0, num=100, endpoint=True, base=10.0)

    dict_ = pickle.load(open(prevalence_dict_path, "rb"))

    ranks_to_plot = ['otu', 'genus', 'subfamily', 'family'] 

    n_cols = 2
    ranks_chunk_all = [ranks_to_plot[x:x+n_cols] for x in range(0, len(ranks_to_plot), n_cols)]

    for ranks_chunk_idx, ranks_chunk in enumerate(ranks_chunk_all):

        for ranks_idx, ranks in enumerate(ranks_chunk):

            ax = plt.subplot2grid((len(ranks_chunk_all), n_cols), (ranks_chunk_idx, ranks_idx))

            for lifestyle in data_utils.lifestyle_all_and_both:
                
                prevalence = numpy.asarray(dict_[lifestyle]['uhgv_taxonomy'][ranks])

                survival_array = data_utils.make_survival_dist(prevalence, prevalence_range)

                ax.plot(prevalence_range, survival_array, lw=2, ls='-', c=data_utils.lifestyle_color_dict[lifestyle], label=lifestyle.capitalize())


            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)
            ax.set_title('UHGV %s scale' % data_utils.uhgv_taxonomy_name_no_cap_dict[ranks])
            ax.set_xlabel('Prevalence across hosts', fontsize=12)
            ax.set_ylabel('Fraction of taxa ' + r'$\geq$', fontsize=12)


            ax.set_xlim([min(prevalence_range), max(prevalence_range)])


            if ranks_idx + ranks_chunk_idx == 0:
                ax.legend(loc='upper right', fontsize=8)


            # run ks test
            prevalence_lytic = numpy.asarray(dict_['lytic']['uhgv_taxonomy'][ranks])
            prevalence_temperate = numpy.asarray(dict_['temperate']['uhgv_taxonomy'][ranks])

            statistic, p_value = data_utils.permutation_2_sample_ks_test(prevalence_lytic, prevalence_temperate)

            #print(statistic, p_value )

            pbalue_label = data_utils.get_latex_pvalue(p_value)

            ax.text(0.75, 0.7, r'$\mathrm{KS} = $' + str(round(statistic, 3)), fontsize=10, ha='center', va='center', transform=ax.transAxes)
            ax.text(0.75, 0.6, pbalue_label, fontsize=10, ha='center', va='center', transform=ax.transAxes)





    
    fig.subplots_adjust(hspace=0.3, wspace=0.25)
    fig_name = "%sprevalence.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()





#make_prevalence_dict()

plot_prevalence()



#ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)

