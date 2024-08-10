

def add_delete_and_gaps(block_object, block_seq, target_genome_name, strand):

    #print(block_object.keys())

    block_mutate_target_genome = [x[1] for x in block_object['mutate'] if x[0]['name'] == target_genome_name][0]
    block_delete_target_genome = [x[1] for x in block_object['delete'] if x[0]['name'] == target_genome_name][0]
    #block_gaps_target_genome = [x[1] for x in block_object['gaps'] if x[0]['name'] == target_genome_name][0]
    block_insert_target_genome = [x[1] for x in block_object['insert'] if x[0]['name'] == target_genome_name][0]

    if len(block_delete_target_genome) == 0:
        return block_seq
    
    block_seq_w_delete = copy.copy(block_seq)

    if strand == False:
        block_seq_w_delete = block_seq_w_delete[::-1]

    list_block_seq_w_delete = list(block_seq_w_delete)

    # apply mutations
    for b in block_mutate_target_genome:
        list_block_seq_w_delete[b[0] - 1] = b[1]

    # apply deletions (insert gaps)
    # for pos, L in dels:
    for b in block_delete_target_genome:
        for l in range(b[1]):
            list_block_seq_w_delete[b[0] - 1 + l] = "-"


    # create dictionary of gaps, to later be inserted
    block_gaps = block_object['gaps']
    gap_dict = {}
    for pos, L in block_gaps.items():
        gap_dict[int(pos)] = ["-"] * L


    # fill these gaps with insertions
    #for ins_descr, nts in ins:
    # tag bases that are inserted
    for b in block_insert_target_genome:
        gap_id, gap_pos = b[0]
        gap = gap_dict[gap_id]  # capture gap
        for i, nt in enumerate(b[1]):
            #gap[gap_pos + i] = nt
            gap[gap_pos + i] = '_'


    # add the gaps to the sequence
    for i, gap in gap_dict.items():
        gap_nt = "".join(gap)
        if i > 0:  # add after position i (julia) or after i-1 (python)
            list_block_seq_w_delete[i - 1] = list_block_seq_w_delete[i - 1] + gap_nt
        else:  # add at the beginning
            list_block_seq_w_delete[0] = gap_nt + list_block_seq_w_delete[0]


    #list_block_seq_w_delete = [x for x in list_block_seq_w_delete if x != "-"]
    block_seq_w_delete = ''.join(list_block_seq_w_delete)
    block_seq_w_delete = block_seq_w_delete.replace('-', '')

    # return in the same direction so you do not forget
    if strand == False:
        block_seq_w_delete = block_seq_w_delete[::-1]


    return block_seq_w_delete





def add_gaps_to_block(block_object):

    block_gaps = block_object['gaps']
    block_sequence = block_object['sequence']

    block_i_sequence_w_gaps = copy.copy(block_sequence)
    cumulative_gaps_added = 0
    for gap_position, gap_size in block_gaps.items():
        
        gap_position_int = int(gap_position) + cumulative_gaps_added
        gap_k = "".join(['-']*gap_size) 

        block_i_sequence_w_gaps = block_i_sequence_w_gaps[:gap_position_int] + gap_k + block_i_sequence_w_gaps[gap_position_int:]
        cumulative_gaps_added += gap_size

    # calculate cumulative number 

    cumulative_n_gaps = numpy.cumsum(numpy.asarray(list(block_i_sequence_w_gaps)) == '-' )

    return block_i_sequence_w_gaps, cumulative_n_gaps



