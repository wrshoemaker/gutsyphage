


def seq_for_back(s, fb):
    
    if fb==True:
        s1=s
    else:
        s1=s[::-1]
    
    return s1



def plot_syntheny(strain1, strain2):
    isolate_paths={strain1: {}, strain2: {}}
    for i in range(len(G['paths'])):
        if G['paths'][i]['name'] in [strain1, strain2]:
            for p in G['paths'][i]['blocks']:
                isolate_paths[G['paths'][i]['name']][p['id']]=p['strand']
    common_blocks={}
    for idx, block in enumerate(G['blocks']):
        block_name=block['id']
        seq_len=len(block['sequence'])
        name_list=[b[0]['name'] for b in block['positions']]
        if strain1 in name_list and strain2 in name_list:
            common_blocks[block_name]={}
            common_blocks[block_name]['seq_len']=seq_len
            common_blocks[block_name][strain1]={'positions':[], 'mutations':[], 'strand': isolate_paths[strain1][block_name]}
            common_blocks[block_name][strain2]={'positions':[], 'mutations':[], 'strand': isolate_paths[strain2][block_name]}
            for b in block['positions']:
                name=b[0]['name']
                if name in [strain1, strain2]:
                    common_blocks[block_name][name]['positions']=b[1]
            for b in block['mutate']:
                name=b[0]['name']
                if name in [strain1, strain2]:
                    common_blocks[block_name][name]['mutations']=b[1]
            set1=create_set(common_blocks[block_name][strain1]['mutations'])
            set2=create_set(common_blocks[block_name][strain2]['mutations'])
            common_blocks[block_name]['different_muts']= set1.symmetric_difference(set2)
            #snps=np.zeros((seq_len))
            #for pos in common_blocks[block_name]['different_muts']:
            #    snps[pos[0]]=1
            #common_blocks[block_name]['snps']=snps
            common_blocks[block_name]['pairwise_divergence']=len(common_blocks[block_name]['different_muts'])/seq_len
    max1=genome_lengths[strain1]
    max2=genome_lengths[strain2]
    for b in common_blocks:
            x=common_blocks[b][strain1]['positions']
            y=common_blocks[b][strain2]['positions']
            if x[0]<x[1] and y[0]<y[1]:
                x_data=[x]
                y_data=[y]
            elif x[0]>x[1] and y[0]<y[1]:
                xsub1=[x[0], max1]
                xsub2=[1, x[1]]
                ysub1=[y[0], y[0]+max1-x[0]]
                ysub2=[y[0]+max1-x[0]+1, y[1]]
                x_data=[xsub1, xsub2]
                y_data=[ysub1, ysub2]
            elif x[0]<x[1] and y[0]>y[1]:
                ysub1=[y[0], max2]
                ysub2=[1, y[1]]
                xsub1=[x[0], x[0]+max2-y[0]]
                xsub2=[x[0]+max2-y[0]+1, x[1]]
                x_data=[xsub1, xsub2]
                y_data=[ysub1, ysub2]
            else:
                delta1=max1-x[0]
                delta2=max2-y[0]
                if delta1<=delta2:
                    xsub1=[x[0], max1]
                    xsub2=[1,max2-(y[0]+max1-x[0])]
                    xsub3=[max2-(y[0]+max1-x[0])+1, x[1]]
                    ysub1=[y[0], y[0]+max1-x[0]]
                    ysub2=[y[0]+max1-x[0]+1, max2]
                    ysub3=[1, y[1]]
                else:
                    ysub1=[y[0], max2]
                    ysub2=[1,max1-(x[0]+max2-y[0])]
                    ysub3=[max1-(x[0]+max2-y[0])+1, y[1]]
                    xsub1=[x[0], x[0]+max2-y[0]]
                    xsub2=[x[0]+max2-y[0]+1, max1]
                    xsub3=[1, x[1]]
                x_data=[xsub1, xsub2, xsub3]
                y_data=[ysub1, ysub2, ysub3]
            for i in range(len(x_data)):
                plt.plot(seq_for_back(x_data[i], common_blocks[b][strain1]['strand']), seq_for_back(y_data[i], common_blocks[b][strain2]['strand']))
    plt.xlabel(strain1)
    plt.ylabel(strain2)
    plt.show()
    return isolate_paths, common_blocks
#%%
isolate_paths, common_blocks=plot_syntheny('UHGV-1353779', 'UHGV-1357691')