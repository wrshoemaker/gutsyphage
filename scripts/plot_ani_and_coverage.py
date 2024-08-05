

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
        pid	= float(line[3])/100
        qcov = float(line[4])/100
        tcov = float(line[5])/100

        blastani_dict[line[0]][line[1]] = {}
        blastani_dict[line[0]][line[1]]['num_alns'] = num_alns
        blastani_dict[line[0]][line[1]]['pid'] = pid
        blastani_dict[line[0]][line[1]]['qcov'] = qcov
        blastani_dict[line[0]][line[1]]['tcov'] = tcov


    file_open.close()


    return blastani_dict



fig = plt.figure(figsize = (8, 4))
fig.subplots_adjust(bottom= 0.15)

ax_pid = plt.subplot2grid((1, 2), (0, 0), colspan=1)
ax_qcov = plt.subplot2grid((1, 2), (0, 1), colspan=1)


blastani_dict = parse_blastani_output(file_path)
#votu_dict = data_utils.read_phage_metadata()
genomes_to_plot = list(blastani_dict.keys())

frame_of_ref_genome = genomes_to_plot[0]

# sort taxa according to frame_of_ref_genome

tuple_frame_pid = [(c, blastani_dict[frame_of_ref_genome][c]['pid']) for c in blastani_dict[frame_of_ref_genome].keys()]
tuple_frame_pid.sort(key=lambda x: x[1], reverse=True)
sorted_genomes = [c[0] for c in tuple_frame_pid]

x_axis_ticks = list(range(len(sorted_genomes)))

colors = ['dodgerblue', '#FF6347']

for g_idx, g in enumerate(genomes_to_plot):
    
    x_axis_to_plot = []
    pid_to_plot = []
    qcov_to_plot = []

    for s_idx, s in enumerate(sorted_genomes):
        
        if g == s:
            continue
        
        x_axis_to_plot.append(s_idx)
        pid_to_plot.append(blastani_dict[g][s]['pid'])
        qcov_to_plot.append(blastani_dict[g][s]['qcov'])


    ax_pid.plot(x_axis_to_plot, pid_to_plot, c=colors[g_idx], lw=1, alpha=0.9, label=g)
    ax_qcov.plot(x_axis_to_plot, qcov_to_plot, c=colors[g_idx], lw=1, alpha=0.9)


ax_pid.set_xticks(x_axis_ticks)
ax_pid.set_xticklabels(sorted_genomes, fontsize=4, rotation=90, ha="right")
#ax_pid.set_yticklabels(sorted_genomes, fontsize=6)

ax_qcov.set_xticks(x_axis_ticks)
ax_qcov.set_xticklabels(sorted_genomes, fontsize=4, rotation=90, ha="right")


ax_pid.set_ylabel('ANI', fontsize=12)
ax_qcov.set_ylabel('Query coverage', fontsize=12)


ax_pid.legend(loc='upper right', fontsize=8)


#print(blastani_dict[frame_of_ref_genome].keys())

# pick reference


#uhgv_genome_all = votu_dict[votu]['uhgv_genome']




fig.subplots_adjust(hspace=0.3, wspace=0.25)
fig_name = "%stest_blast_ani.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()



#for g in genomes_to_plot:

#for votu in votu_dict.keys():


#    genomes_to_plot

#    print()
