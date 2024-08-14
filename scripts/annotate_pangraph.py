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






#genomes_to_ignore = ['UHGV-1372371']






def run_everything():

    #parse_annotation_table()

    votu_all = [x.split('_')[0] for x in os.listdir('%scomplete_minimap2/' % config.data_directory)]
    votu_all.sort()
    #votu_all = numpy.asarray(votu_all)
    #proble_votu_idx = votu_all.index('vOTU-000118')
    #data_utils.build_votu_fasta('vOTU-000050')

    votus_to_skip = ['vOTU-000118', 'vOTU-000016', 'vOTU-000035', 'vOTU-000118']
    
    for votu_idx, votu in enumerate(votu_all):
        # make fasta file
        #f votu == 'vOTU-000118':
        #    continue
        #data_utils.build_votu_fasta(votu)

        if votu in votus_to_skip:
            continue
        
        # make syn dict
       # data_utils.make_syn_sites_votu_dict(votu)

       






#run_everything()


#data_utils.make_syn_sites_votu_dict(votu)







