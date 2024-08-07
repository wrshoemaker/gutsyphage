import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config


import json


fasta_votu = data_utils.classFASTA('%svOTU-000085.fna' % config.data_directory).readFASTA()
fasta_votu_dict = {x[0]:x[1] for x in fasta_votu}

pangraph_path = '%svOTU-000085-polished-pangraph.json' % config.data_directory
pangraph_file = open(pangraph_path, 'r')





pangraph_data = json.load(pangraph_file)

# blocks
### id = str
### sequence = AACCF
### positions = where on the genome the block is from
### 

pangraph_file.close()

genome_1 = 'UHGV-1357691'
#genome_1 = 'UHGV-1350565'
genome_2 = 'UHGV-1353779'
#genome_2 = 'UHGV-1350565'
#genome_2 = 'UHGV-1355338'


#votu_dict, vgenome_dict = data_utils.read_uhgv_metadata()

#for g in fasta_votu_dict.keys():

#    print(g, vgenome_dict[g]['original_id'])





genome_1_size = len(fasta_votu_dict[genome_1])
genome_2_size = len(fasta_votu_dict[genome_2])



fig = plt.figure(figsize = (4, 4))
fig.subplots_adjust(bottom= 0.15)

ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)


for block_i in pangraph_data['blocks']:

    positions = block_i['positions']

    genomes_in_positions = [p[0]['name'] for p in positions]
    # start and stop of each item in positions
    #start_stop_in_positions = [p[1] for p in positions]

    if (genome_1 in genomes_in_positions) and (genome_2 in genomes_in_positions):

        position_1_idx = genomes_in_positions.index(genome_1)
        position_2_idx = genomes_in_positions.index(genome_2)

        start_1, stop_1 = positions[position_1_idx][1]
        start_2, stop_2 = positions[position_2_idx][1]

        strand_1 = positions[position_1_idx][0]['strand']
        strand_2 = positions[position_2_idx][0]['strand']

        
        if (strand_1 == False):
            stop_1, start_1 = start_1, stop_1

        if (strand_2 == False):   
            stop_2, start_2 = start_2, stop_2


        if (start_1 > stop_1) and (start_2 <= stop_2):    
            start_stop_all_1 = [(start_1, genome_1_size), (0, stop_1)]
            fraction_1_genome_start = (genome_1_size - stop_2 + start_2)/start_2
            # 
            #start_stop_all_2 = [(start_2, stop_2*fraction_1_genome_start), (stop_2*fraction_1_genome_start, stop_2)]
            #start_stop_all_2 = start_stop_all_2[::-1]
            #print(start_stop_all_1, start_stop_all_2)

            for i in range(len(start_stop_all_1)):
                ax.plot(start_stop_all_1[i], start_stop_all_2[i], lw=2, ls='-', c='dodgerblue')



        elif (start_1 <= stop_1) and (start_2 > stop_2):
            start_stop_all_1 = [(start_1, stop_1 - stop_2), (stop_1 - stop_2, stop_1)]
            start_stop_all_2 = [(start_2, genome_2_size), (0, stop_2)]

        elif (start_1 > stop_1) and (start_2 > stop_2):
            start_stop_all_1 = [(0, stop_1), (start_1, genome_1_size)]
            start_stop_all_2 = [(0, stop_2), (start_2, genome_2_size)]

        else:
            start_stop_all_1 = [(start_1, stop_1)]
            start_stop_all_2 = [(start_2, stop_2)]


        
        ax.set_xlabel('%s position' % genome_1)
        ax.set_ylabel('%s position' % genome_2)



fig.subplots_adjust(hspace=0.3, wspace=0.25)
fig_name = "%sdotplot_%s_%s.png" % (config.analysis_directory, genome_1, genome_2)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()


