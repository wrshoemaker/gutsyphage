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


uhgv_votu_metadata_dict, uhgv_genome_metadata_dict = data_utils.read_uhgv_metadata(checkv_quality='Low-quality', checkv_quality_cumulative=True)
uhgv_votu_extended_metadata_dict = data_utils.read_uhgv_votu_metadata()




genome_all = sorted(list(uhgv_genome_metadata_dict.keys()))

#print(uhgv_votu_extended_metadata_dict)

#for g in genome_all:

#    print(uhgv_votu_extended_metadata_dict[ uhgv_genome_metadata_dict[g]['uhgv_votu'] ]['lifestyle']  )
#    print(uhgv_votu_extended_metadata_dict[uhgv_genome_metadata_dict[g]['uhgv_votu']]['lifetyle'])



lifestyle_all = numpy.asarray([uhgv_votu_extended_metadata_dict[uhgv_genome_metadata_dict[g]['uhgv_votu']]['lifestyle'] for g in genome_all ])
checkv_completeness_all = numpy.asarray([uhgv_genome_metadata_dict[g]['checkv_completeness'] for g in genome_all ])

fig = plt.figure(figsize = (4, 4))
fig.subplots_adjust(bottom= 0.15)

ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)


completeness_range = numpy.logspace(numpy.log10(min(checkv_completeness_all)), 0, num=50, endpoint=True, base=10.0)

for l_idx, l in enumerate(data_utils.lifestyle_all_and_both):

    print(l)

    if l != 'both':
        checkv_completeness_l = checkv_completeness_all[lifestyle_all == l]

    else:
        checkv_completeness_l = numpy.copy(checkv_completeness_all)


    survival_array = data_utils.make_survival_dist(checkv_completeness_l, completeness_range)
    ax.plot(completeness_range, survival_array, lw=2, ls='-', c=data_utils.lifestyle_color_dict[l], label=l.capitalize())


#ax.set_xscale('log', basex=10)
#ax.set_yscale('log', basey=10)

ax.set_xlabel('Genome completeness score', fontsize=12)
ax.set_ylabel('Fraction of genomes ' + r'$\geq$', fontsize=12)

ax.legend(loc='lower left', fontsize=8)


fig.subplots_adjust(hspace=0.3, wspace=0.25)
fig_name = "%scompleteness_lifestyle.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()
