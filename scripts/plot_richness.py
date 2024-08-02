
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats

import data_utils
import config


n_cols = 3
taxon_ranks = data_utils.taxon_ranks

sra_phage_dict = data_utils.read_phage_metadata()
sra_microbe_dict = data_utils.read_microbe_metadata()

sample_intersect = numpy.intersect1d(list(sra_phage_dict.keys()), list(sra_microbe_dict.keys()))
phage_richness = numpy.asarray([len(sra_phage_dict[s]) for s in sample_intersect])

# weird outlier
#print(sample_intersect[numpy.argmax(phage_richness)])

taxon_chunk_all = [taxon_ranks[x:x+n_cols] for x in range(0, len(taxon_ranks), n_cols)]

fig = plt.figure(figsize = (12, 8))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for taxon_chunk_idx, taxon_chunk in enumerate(taxon_chunk_all):

    for taxon_idx, taxon in enumerate(taxon_chunk):

        ax = plt.subplot2grid((len(taxon_chunk_all), n_cols), (taxon_chunk_idx, taxon_idx))

        microbe_richness_taxon = [len(set([sra_microbe_dict[s][g]['taxonomy'][taxon] for g in sra_microbe_dict[s].keys()])) for s in sample_intersect]
        microbe_richness_taxon = numpy.asarray(microbe_richness_taxon)

        to_plot_idx = phage_richness < 200

        microbe_richness_taxon_to_plot = microbe_richness_taxon[to_plot_idx]
        phage_richness_to_plot = phage_richness[to_plot_idx]

        rho = numpy.corrcoef(microbe_richness_taxon_to_plot, phage_richness_to_plot)[0,1]

        merged_array = numpy.concatenate([microbe_richness_taxon_to_plot, phage_richness_to_plot])
        #min_ = min(merged_array)
        max_ = max(merged_array)*1.05
        
        max_x = max(microbe_richness_taxon_to_plot)*1.05
        max_y = max(phage_richness_to_plot)*1.05

        ax.scatter(microbe_richness_taxon_to_plot, phage_richness_to_plot, alpha=0.3, s=12, zorder=2)
        ax.plot([0, max_x], [0, max_y], lw=2, ls=':', c='k', zorder=1, label='1:1')

        ax.set_xlim([0, max_x])
        ax.set_ylim([0, max_y])


        ax.text(0.2, 0.75, r'$\rho = $' + str(round(rho, 3)), fontsize=12, ha='center', va='center', transform=ax.transAxes)


        ax.set_xlabel('# microbial genomes', fontsize=12)
        ax.set_ylabel('# phage genomes', fontsize=12)
        ax.set_title('%s-level' % taxon.capitalize(), fontsize=12)

        if taxon_idx + taxon_chunk_idx == 0:
            ax.legend(loc='upper left')
        
        

        #ax.set_xscale('log', basex=10)
        #ax.set_yscale('log', basey=10)



fig.subplots_adjust(hspace=0.3, wspace=0.25)
fig_name = "%sgenome_richness.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()