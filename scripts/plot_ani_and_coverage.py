

import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config



file_path = '%stest_blast_ani.txt' % config.data_directory


def parse_blastani_output(file_path, min_num_alns=0):

    blastani_dict = {}

    file_open = open(file_path, 'r')

    header = file_open.readline()
    header = header.strip().split('\t')

    for line in file_open:
        line = line.strip().split('\t')

        if line[0] not in blastani_dict:
            blastani_dict[line[0]] = {}

        num_alns = int(line[2])
        pid	= float(line[3])
        qcov = float(line[4])
        tcov = float(line[5])

        blastani_dict[line[0]][line[1]] = {}
        blastani_dict[line[0]][line[1]]['num_alns'] = num_alns
        blastani_dict[line[0]][line[1]]['pid'] = pid
        blastani_dict[line[0]][line[1]]['qcov'] = qcov
        blastani_dict[line[0]][line[1]]['tcov'] = tcov


    file_open.close()


    return blastani_dict



blastani_dict = parse_blastani_output(file_path)
votu_dict = data_utils.read_phage_metadata()
genomes_to_plot = list(blastani_dict.keys())

genome_to_votu_dict = {}

#uhgv_genome_all = votu_dict[votu]['uhgv_genome']



#for g in genomes_to_plot:

#for votu in votu_dict.keys():


#    genomes_to_plot

#    print()
