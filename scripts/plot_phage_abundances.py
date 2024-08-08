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




file_path = '%sshkoporov_2019_VLP_stool_vOTU_abundance.csv' % config.data_directory
file_open = open(file_path, 'r')
header = file_open.readline()
header_split = header.strip().split(',')
votu_labels = header_split[6:]

is_vlp = numpy.asarray([True if 'VLP'  in v else False for v in votu_labels])

# VLP = virus-like particle

sad_dict = {}
for line in file_open:

    line_split = line.strip().split(',')
    host = line_split[0]

    sad = numpy.asarray([float(c) for c in line_split[6:]])

    sad_vlp = sad[is_vlp]
    sad_stool = sad[~is_vlp]

    #sad_vlp = sad_vlp/sum(sad_vlp)
    #sad_stool = sad_stool/sum(sad_stool)

    sad_dict[host] = {}
    sad_dict[host]['sad_vlp'] = sad_vlp
    sad_dict[host]['sad_stool'] = sad_stool


file_open.close()


hosts_all = numpy.asarray(list(sad_dict.keys()))
n_cols = 3
hosts_chunk_all = [hosts_all[x:x+n_cols] for x in range(0, len(hosts_all), n_cols)]

fig = plt.figure(figsize = (8, 9))
fig.subplots_adjust(bottom= 0.15)


for hosts_chunk_idx, hosts_chunk in enumerate(hosts_chunk_all):

    for host_idx, host in enumerate(hosts_chunk):

        ax = plt.subplot2grid((len(hosts_chunk_all), n_cols), (hosts_chunk_idx, host_idx))

        sad_stool = sad_dict[host]['sad_stool']
        sad_vlp = sad_dict[host]['sad_vlp']

        ax.scatter(sad_stool, sad_vlp, s=10, alpha=0.8, c='k', zorder=2)

        #for lifestyle in data_utils.lifestyle_all_and_both:
            
        #    prevalence = numpy.asarray(dict_[lifestyle]['uhgv_taxonomy'][ranks])

        #    survival_array = data_utils.make_survival_dist(prevalence, prevalence_range)

        #    ax.plot(prevalence_range, survival_array, lw=2, ls='-', c=data_utils.lifestyle_color_dict[lifestyle], label=lifestyle.capitalize())

        sad_merged = numpy.concatenate((sad_stool, sad_vlp))
        min_, max_ = min(sad_merged), max(sad_merged)
        #ax.set_xscale('log', basex=10)
        #ax.set_yscale('log', basey=10)
        ax.set_title(host, fontsize=9)
        #ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='both', labelsize=3)


        if  max(sad_vlp) == 0:
            max_y = max(sad_stool)
        else:
            max_y = max(sad_vlp)
        
        ax.plot([min_, max_ ], [min_, max_], ls=':', lw=1, c='k', zorder=1, label='1:1')

        ax.set_xlim([min(sad_stool), max(sad_stool)* 1.1 ])
        print(sum(sad_vlp))
        ax.set_ylim([min(sad_vlp), max_y*1.1])

        ax.set_xlabel('Relative abundance, stool', fontsize=8)
        ax.set_ylabel('Relative abundance, VLP', fontsize=8)

        if host_idx + hosts_chunk_idx == 0:
            ax.legend(loc='upper left', fontsize=7)

        #ax.set_ylabel('Fraction of taxa ' + r'$\geq$', fontsize=12)



fig.subplots_adjust(hspace=0.45, wspace=0.45)
fig_name = "%sphage_abundance.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()
# plot