
import gzip
import config
import os
from math import radians, cos, sin, asin, sqrt

from collections import Counter

import numpy
from scipy import stats
import pandas

from matplotlib import colors
from matplotlib import cm
import matplotlib as mpl


from Bio import SeqIO


n_fna_characters = 80

numpy.random.seed(123456789)


# these repeated vOTUs are not an issue if you only use MGV lines
repeated_votu_representative_all = ['vOTU-010833', 'vOTU-014541', 'vOTU-013486']
# ignore 'Not-determined'
checkv_quality_all = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality']

checkv_quality_score_dict = {'Low-quality':0, 'Medium-quality':1, 'High-quality':2, 'Complete':3}

taxon_abbr_dict = {'d':'domain', 'p':'phylum', 'c':'class', 'o':'order', 'f':'family', 'g':'genus', 's':'species'}
taxon_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
#taxon_ranks = taxon_ranks[::-1]

ictv_taxonomy_ranks = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
uhgv_taxonomy_ranks = ['family', 'subfamily', 'genus', 'otu']
phage_taxonomy_ranks_dict = {'ictv_taxonomy': ictv_taxonomy_ranks, 'uhgv_taxonomy': uhgv_taxonomy_ranks}


lifestyle_all = ['lytic', 'temperate']
lifestyle_all_and_both = ['lytic', 'temperate', 'both']

lifestyle_color_dict = {'both':'k', 'temperate':'dodgerblue', 'lytic':'#FF6347'}

uhgv_taxonomy_name_no_cap_dict = {'otu':'vOTU', 'genus': 'genus', 'subfamily':'subfamily', 'family':'family'}
str_to_boolean_dict = {'Yes': True, 'No': False}
country_code_dict = {'RUS': 'Russia', 'MNG':'Mongolia', 'FRA':'France', 'FIN':'Finland', 'AUT':'Austria', 'DEU': 'Denmark', 'GBR':'United Kingdom', 'EST':'Estonia', 'CAN':'Canada', 'PER':'Peru', 'ITA':'Italy', 'JPN':'Japan', 'FJI':'Fiji', 'BGD':'Bangladesh', 'SWE':'Sweeden', 'CHN':'China', 'KAZ':'Kazakistan', 'NLD':'Netherlands', 'ESP':'Spain', 'DNK':'Denmark', 'SLV':'El Salvador', 'ISR':'Israel', 'TZA':'Tanzania', 'USA': 'USA'}





def get_latex_pvalue(p_value):

    if p_value < 0.05:

        label = r'$P < 0.05$'

    else:
        label = r'$P \nleq 0.05$'


    
    return label


def make_colormap(n_entries):
    # some matplotlib tools
    # taxonomic hierarchy colors
    cmap_offset = int(0.2*16)
    # +cmap_offset
    rgb_red_ = cm.Blues(numpy.linspace(0,1,n_entries+5))
    rgb_red_ = mpl.colors.ListedColormap(rgb_red_[cmap_offset:,:-1])

    return rgb_red_



def read_sample_metadata():

    sample_metagenome_dict = {}

    with gzip.open('%smgv_sample_info.tsv.gz' % config.data_directory, 'rt') as f:
        header = f.readline()
        header = header.strip().split('\t')

        # line[0] = contig_id
        # line[1] = assembly_source. Sequence database
        # line[2] = assembly_name (e.g., SRR3740115). May also be in uhgv_metadata.tsv
        # line[3] = study_accession (e.g., SRP077685)
        # line[4] = sample_accession (e.g., SRS1537996). May also be in uhgv_metadata.tsv
        # line[5] = run_accessions (e.g., SRR3740115). May also be in uhgv_metadata.tsv
        # line[6] = continent
        # line[7] = country_code
        # line[8] = sex
        # line[9] = age
        # line[10] = health
        # line[11] = disease

        for line in f:
            line = line.strip().split('\t')

            contig_id = line[0]
            assembly_source = line[1]
            assembly_name = line[2]
            study_accession = line[3]
            sample_accession = line[4]
            run_accessions = line[5]

            continent = line[6]
            country_code = line[7]
            sex = line[8]
            age = line[9]

            health = line[10]
            disease = line[11]

            sample_metagenome_dict[contig_id] = {}
            sample_metagenome_dict[contig_id]['assembly_source'] = assembly_source
            sample_metagenome_dict[contig_id]['assembly_name'] = assembly_name
            sample_metagenome_dict[contig_id]['study_accession'] = study_accession
            sample_metagenome_dict[contig_id]['sample_accession'] = sample_accession

            # run_accessions will be a LIST
            sample_metagenome_dict[contig_id]['run_accessions'] = run_accessions.split(',')

            sample_metagenome_dict[contig_id]['continent'] = continent

            sample_metagenome_dict[contig_id]['country_code'] = country_code
            sample_metagenome_dict[contig_id]['sex'] = sex
            sample_metagenome_dict[contig_id]['age'] = age

            sample_metagenome_dict[contig_id]['health'] = health
            sample_metagenome_dict[contig_id]['disease'] = disease

            #sample_metagenome_dict[]


        f.close()

    
    return sample_metagenome_dict




def read_microbe_metadata(min_completeness=0.9, max_contamination=0.05):

    #study_accession_all = []
    sra_dict = {}

    with gzip.open('%suhgg_genomes-nr_metadata.tsv.gz' % config.data_directory, 'rt') as f:
        header = f.readline()
        header = header.strip().split('\t')

        # line[0] = Genome ID
        # line[1] = Original name (e.g., 11861_6_55)
        # line[2] = Study set (e.g., HBC, NCBI). 
        # line[3] = Genome_type (e.g., Isolate). 
        # line[4] = Length
        # line[5] = N_contigs
        # line[6] = N50
        # line[7] = GC_content
        # line[8] = Completeness
        # line[9] = Contamination
        # line[10] = CMseq (Not sure)
        # line[11] = rRNA_5S (completeness)
        # line[12] = rRNA_16S (completeness)
        # line[13] = rRNA_23S (completeness)
        # line[14] = tRNAs (#)
        # line[15] = Genome_accession (e.g., GCA_003480885)
        # line[16] = Species_rep (e.g., GUT_GENOME000001)
        # line[17] = MGnify_accession (e.g., MGYG-HGUT-00001)
        # line[18] = Lineage (taxonomy)
        # line[19] = Sample_accession (e.g., ERS370061, SRS1992821)
        # line[20] = Study_accession (e.g., ERP105624)
        # line[21] = Country
        # line[22] = Continent
        # line[23] = FTP_download


        for line in f:
            line = line.strip().split('\t')


            genome_id = line[0]
            completeness = float(line[8])/100
            contamination = float(line[9])/100
            sra_sample = line[19]
            lineage = line[18]

            if completeness < min_completeness:
                continue

            if contamination > max_contamination:
                continue


            if sra_sample not in sra_dict:
                sra_dict[sra_sample] = {}

            # genome_id is unique
            sra_dict[sra_sample][genome_id] = {}
            sra_dict[sra_sample][genome_id]['taxonomy'] = {}

            for t in lineage.split(';'):

                t_abbr, t_name = t.split('__')

                if t_name == '':
                    t_name = 'NA'

                sra_dict[sra_sample][genome_id]['taxonomy'][taxon_abbr_dict[t_abbr]] = t_name
                
            #sra_dict[sra_sample][]


        f.close()

    
    return sra_dict



def read_uhgv_metadata(checkv_quality='High-quality', checkv_quality_cumulative=True, viral_confidence_criterion='Confident', viralverify_prediction_criterion='Virus', checkv_completeness=0.9):

    #  subset_votu_representative=False, 

    #sra_dict = {}
    votu_dict = {}
    vgenome_dict = {}

    with gzip.open('%suhgv_metadata.tsv.gz' % config.data_directory, 'rt') as f:
        header = f.readline()
        header = header.strip().split('\t')

        # line[0] = uhgv_genome (e.g., UHGV-1495677)
        # line[1] = uhgv_votu (e.g., vOTU-000135)
        # line[2] = votu_representative (Yes/No)
        # line[3] = original_study_alias (e.g., MGV)
        # line[4] = original_id (e.g., MGV-GENOME-0354054)
        # line[5] = genome_length
        # line[6] = is_jumbo_phage (Yes/No)
        # line[7] = gc_percent
        # line[8] = checkv_quality (e.g., High-quality)
        # line[9] = checkv_completeness (%)
        # line[10] = checkv_completeness_method (e.g., AAI-based)
        # line[11] = checkv_trimmed (Yes/No)
        # line[12] = viral_confidence (e.g., Confident)
        # line[13] = genomad_virus_score 
        # line[14] = genomad_virus_hallmarks (represents if the detection method for the phages (genomad) found a phage hallmark gene (eg. capsid, terminase, etc.))
        # line[15] = genomad_plasmid_hallmarks (?)
        # line[16] = viralverify_prediction (e.g., Virus)
        # line[17] = viralverify_score (numeric)
        # line[18] = checkv_viral_markers (numeric)
        # line[19] = checkv_host_markers (numeric)
        # line[20] = cds_count (numeric)
        # line[21] = cds_density (numeric) ?
        # line[22] = avg_cds_length (numeric)
        # line[23] = genetic_code (numeric)
        # line[24] = is_recoded (Yes/No)
        # line[25] = trna_count_total
        # line[26] = trna_count_suppressor (A suppressor tRNA is a tRNA with a mutation (usually) in the anticodon that allows it to recognize a stop codon and insert an amino acid in its place) 10.1128/JB.186.20.6714-6720.2004
        # line[27] = sra_run (e.g., ERR414592,ERR414311 or SRR5274010)
        # line[28] = sra_sample (e.g., ERS396364 or SRS1992821)
        # line[29] = biosample (e.g., SAMEA2338694)
        # line[30] = country
        # line[31] = latitude
        # line[32] = longitude

        for line in f:
            line = line.strip().split('\t')

            uhgv_genome = line[0]
            uhgv_votu = line[1]

            original_study_alias = line[3]
            original_id = line[4]
            genome_length = line[5]

            quality = line[8]
            checkv_completeness = line[9]

            viral_confidence = line[12]
            viralverify_prediction = line[16]
            sra_sample = line[28]

            country = line[30]
            latitude = line[31]
            longitude = line[32]
            
            if quality == 'Not-determined':
                continue

            if checkv_completeness == 'NULL':
                continue
            
            checkv_completeness = float(checkv_completeness)/100
            
            checkv_quality_score = checkv_quality_score_dict[quality]

            if checkv_quality_cumulative == True:
                if checkv_quality_score < checkv_quality_score_dict[checkv_quality]:
                    continue
            else:
                if (quality != checkv_quality):
                    continue


            if original_study_alias != 'MGV':
                continue
            
            if (viral_confidence != viral_confidence_criterion):
                continue

            if (viralverify_prediction != viralverify_prediction_criterion):
                continue
            
            #if subset_votu_representative == True:
            # make sure it's representative
            #    if str_to_boolean_dict[line[2]] != True:
            #        continue

            if line[28] == 'Null':
                continue


            #checkv_completeness = float(line[9])/100

            #if sra_sample not in sra_dict:
            #    sra_dict[sra_sample] = {}
            #    sra_dict[sra_sample]['latitude'] = latitude
            #    sra_dict[sra_sample]['longitude'] = longitude
            #    sra_dict[sra_sample]['country'] = country
            #    sra_dict[sra_sample]['uhgv_votu_all'] = []


            #sra_dict[sra_sample]['uhgv_votu_all'].append(uhgv_votu)

            if uhgv_votu not in votu_dict:
                votu_dict[uhgv_votu] = {}
                votu_dict[uhgv_votu]['uhgv_genome'] = []
                votu_dict[uhgv_votu]['original_study_alias'] = []
                votu_dict[uhgv_votu]['original_id'] = []
                votu_dict[uhgv_votu]['genome_length'] = []
                votu_dict[uhgv_votu]['quality'] = []
                votu_dict[uhgv_votu]['checkv_completeness'] = []
                votu_dict[uhgv_votu]['viral_confidence'] = []
                votu_dict[uhgv_votu]['viralverify_prediction'] = []
                votu_dict[uhgv_votu]['sra_sample'] = []
                votu_dict[uhgv_votu]['country'] = []


                
            votu_dict[uhgv_votu]['uhgv_genome'].append(uhgv_genome)
            votu_dict[uhgv_votu]['original_study_alias'].append(original_study_alias)
            votu_dict[uhgv_votu]['original_id'].append(original_id)
            votu_dict[uhgv_votu]['genome_length'].append(genome_length)
            votu_dict[uhgv_votu]['quality'].append(quality)
            votu_dict[uhgv_votu]['checkv_completeness'].append(checkv_completeness)
            votu_dict[uhgv_votu]['viral_confidence'].append(viral_confidence)
            votu_dict[uhgv_votu]['viralverify_prediction'].append(viralverify_prediction)
            votu_dict[uhgv_votu]['sra_sample'].append(sra_sample)
            votu_dict[uhgv_votu]['country'].append(country)

            #if line[1] not in sra_dict[sra_sample]:

            # save everything for vgenome_dict

            vgenome_dict[uhgv_genome] = {}
            vgenome_dict[uhgv_genome]['uhgv_votu'] = uhgv_votu
            vgenome_dict[uhgv_genome]['original_study_alias'] = original_study_alias
            vgenome_dict[uhgv_genome]['original_id'] = original_id
            vgenome_dict[uhgv_genome]['genome_length'] = genome_length
            vgenome_dict[uhgv_genome]['quality'] = quality
            vgenome_dict[uhgv_genome]['checkv_completeness'] = checkv_completeness
            vgenome_dict[uhgv_genome]['viral_confidence'] = viral_confidence
            vgenome_dict[uhgv_genome]['viralverify_prediction'] = viralverify_prediction
            vgenome_dict[uhgv_genome]['sra_sample'] = sra_sample
            vgenome_dict[uhgv_genome]['country'] = country

            #    sra_dict[sra_sample][line[1]] = {}
        f.close()


    return votu_dict, vgenome_dict



def read_uhgv_votu_metadata():

    uhgv_genome_metadata_dict = {}

    # one line per genome AND vOTU
    # representative genomes
    with gzip.open('%suhgv_votus_metadata.tsv.gz' % config.data_directory, 'rt') as f:
        
        header = f.readline()
        header = header.strip().split('\t')

        # line[0] = uhgv_genome (e.g., UHGV-1495677)
        # line[1] = uhgv_votu (e.g., vOTU-000135)
        # line[2] = checkv_completeness (0-100)
        # line[3] = checkv_quality (e.g., Complete)
        # line[4] = viral_confidence (e.g., Confident)
        # line[5] = genome_length
        # line[6] = cds_count
        # line[7] = ictv_taxonomy
        # line[8] = uhgv_taxonomy
        # line[9] = host_lineage
        # line[10] = lifestyle
        # line[11] = http_url


        for line in f:
            line = line.strip().split('\t')

            uhgv_genome = line[0]
            uhgv_votu = line[1]
            checkv_completeness = line[2]
            checkv_quality = line[3]
            viral_confidence = line[4]
            genome_length = line[5]
            cds_count = line[6]
            lifestyle = line[10]

            uhgv_genome_metadata_dict[uhgv_votu] = {}
            #uhgv_genome_metadata_dict[uhgv_genome]['uhgv_votu'] = uhgv_votu
            uhgv_genome_metadata_dict[uhgv_votu]['uhgv_genome'] = uhgv_genome
            uhgv_genome_metadata_dict[uhgv_votu]['checkv_completeness'] = checkv_completeness
            uhgv_genome_metadata_dict[uhgv_votu]['checkv_quality'] = checkv_quality
            uhgv_genome_metadata_dict[uhgv_votu]['viral_confidence'] = viral_confidence
            uhgv_genome_metadata_dict[uhgv_votu]['genome_length'] = genome_length
            uhgv_genome_metadata_dict[uhgv_votu]['cds_count'] = cds_count
            uhgv_genome_metadata_dict[uhgv_votu]['lifestyle'] = lifestyle

            ictv_taxonomy_split = line[7].split(';')
            uhgv_taxonomy_split = line[8].split(';')
            host_lineage = line[9]

            uhgv_genome_metadata_dict[uhgv_votu]['ictv_taxonomy'] = {}
            uhgv_genome_metadata_dict[uhgv_votu]['uhgv_taxonomy'] = {}
            uhgv_genome_metadata_dict[uhgv_votu]['host_lineage'] = {}

            if host_lineage == 'NULL':
                for t in taxon_ranks:
                    uhgv_genome_metadata_dict[uhgv_votu]['host_lineage'][t] = 'NA'

            else:
                host_lineage_split = host_lineage.split(';')
                # skip domain, uninformative
                host_lineage_split = [h.split('__')[1] for h in host_lineage_split][1:]

                for taxon_ranks_i_idx, taxon_ranks_i in enumerate(taxon_ranks):
                    if taxon_ranks_i_idx >= len(host_lineage_split):
                        uhgv_genome_metadata_dict[uhgv_votu]['host_lineage'][taxon_ranks_i] = 'NA'

                    else:
                        uhgv_genome_metadata_dict[uhgv_votu]['host_lineage'][taxon_ranks_i] = host_lineage_split[taxon_ranks_i_idx]

            

            for ictv_taxonomy_ranks_i_idx, ictv_taxonomy_ranks_i in enumerate(ictv_taxonomy_ranks):

                if ictv_taxonomy_ranks_i_idx >= len(ictv_taxonomy_split):
                    uhgv_genome_metadata_dict[uhgv_votu]['ictv_taxonomy'][ictv_taxonomy_ranks_i] = 'NA'

                else:
                    uhgv_genome_metadata_dict[uhgv_votu]['ictv_taxonomy'][ictv_taxonomy_ranks_i] = ictv_taxonomy_split[ictv_taxonomy_ranks_i_idx]


            for uhgv_taxonomy_ranks_i_idx, uhgv_taxonomy_ranks_i in enumerate(uhgv_taxonomy_ranks):

                if uhgv_taxonomy_ranks_i_idx >= len(uhgv_taxonomy_split):
                    uhgv_genome_metadata_dict[uhgv_votu]['uhgv_taxonomy'][uhgv_taxonomy_ranks_i] = 'NA'

                else:
                    uhgv_genome_metadata_dict[uhgv_votu]['uhgv_taxonomy'][uhgv_taxonomy_ranks_i] = uhgv_taxonomy_split[uhgv_taxonomy_ranks_i_idx]

            #print(uhgv_genome_metadata_dict[uhgv_votu]['ictv_taxonomy'])

    return uhgv_genome_metadata_dict
        





def votu_dict_to_pres_abs(phage_metadata_dict):

    #uhgv_votu 

    #votu_pres_abs = pandas.DataFrame({k:Counter(v['sra_sample']) for k, v in phage_metadata_dict.items()}).T.fillna(0).astype(int)
    votu_pres_abs = pandas.DataFrame({k:Counter(v['sra_sample']) for k, v in phage_metadata_dict.items()}).T.fillna(0).astype(int)

    votu_pres_abs_array = votu_pres_abs.values   

    sample_names = votu_pres_abs.columns.values
    votu_names = votu_pres_abs.index.values



    # samples are COLUMNS
    # vOTUs are ROWS

    return votu_pres_abs_array, votu_names, sample_names





def coarse_grain_abundances_by_taxonomy(count_array, votus, uhgv_votu_metadata_dict, taxonomic_level='genus', phage_taxonomic_system='uhgv_taxonomy', pres_abs=True):


    taxon_labels_all = [uhgv_votu_metadata_dict[v][phage_taxonomic_system][taxonomic_level] for v in votus]

    taxa_to_keep_idx = numpy.asarray(taxon_labels_all)
    taxa_to_keep_idx = taxa_to_keep_idx != "NA"

    taxon_labels_all_clean = [t for t in taxon_labels_all if t != "NA"]

    # return this list for null coarse-grained site-by-species matrix
    n_collapsed_nodes = list(dict(Counter(taxon_labels_all_clean)).values())

    # remove duplicates from taxon_labels_all_clean while keeping original order
    def unique(sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]

    taxon_labels_all_clean_set = unique(taxon_labels_all_clean)
    taxon_labels_all_clean_set = numpy.asarray(taxon_labels_all_clean_set)
    taxon_labels_all_clean_copy = taxon_labels_all_clean.copy()

    to_remove_idx = [i for i, x in enumerate(taxon_labels_all) if x == "NA"]
    to_remove_idx = numpy.asarray(to_remove_idx)
    count_array_ = numpy.copy(count_array)

    # delete NA
    if len(to_remove_idx) > 0:
        count_array_ = numpy.delete(count_array_, to_remove_idx, axis=0)

    for taxon_i in taxon_labels_all_clean_set:

        taxon_i_idx = [i for i, x in enumerate(taxon_labels_all_clean_copy) if x == taxon_i]
        taxon_i_idx = numpy.asarray(taxon_i_idx)

        # remove taxon from the list
        taxon_labels_all_clean_copy = [s for s in taxon_labels_all_clean_copy if s != taxon_i]

        count_array_to_merge = count_array_[taxon_i_idx,]
        count_array_ = numpy.delete(count_array_, taxon_i_idx, axis=0)

        count_array_sum = numpy.sum(count_array_to_merge, axis=0)
        count_array_ = numpy.vstack((count_array_, count_array_sum))

    
    if pres_abs == True:
        count_array_[count_array_>0] = 1


    return count_array_, taxon_labels_all_clean_set, n_collapsed_nodes, taxa_to_keep_idx






def make_survival_dist(data, range_, probability=True):


    data = data[numpy.isfinite(data)]
    survival_array = [sum(data>=i) for i in range_]
    #survival_array = [sum(data>=i)/len(data) for i in range_]
    survival_array = numpy.asarray(survival_array)

    if probability == True:
        survival_array = survival_array/len(data)

    return survival_array


# estimated host range in uhgv_votus_metadata.tsv.gz



def permutation_2_sample_ks_test(array_1, array_2, n=1000):

    statistic, pvalue = stats.ks_2samp(array_1, array_2, alternative='two-sided')

    array_merged = numpy.concatenate((array_1, array_2), axis=None)
    
    statistic_null = []
    for n_i in range(n):

        numpy.random.shuffle(array_merged)
        statistic_null.append(stats.ks_2samp(array_merged[:len(array_1)], array_merged[len(array_1):], alternative='two-sided')[0])

    statistic_null = numpy.asarray(statistic_null)

    p_value = sum(statistic_null > statistic)/n

    return statistic, p_value




def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r


def build_votu_fasta(votu='vOTU-000085'):

    uhgv_votu_metadata_dict, uhgv_genome_metadata_dict = read_uhgv_metadata(checkv_quality='Complete', checkv_quality_cumulative=True)
    
    target_genome_all = uhgv_votu_metadata_dict[votu]['uhgv_genome']


    #fasta_all_genomes = classFASTA('%suhgv_mgv.fna.gz' % config.data_directory).readFASTA()

    fasta_file_path = '%s%s.fna' % (config.data_directory, votu) 
    fasta_file_open = open(fasta_file_path, 'w')

    with gzip.open('%suhgv_mgv.fna.gz' % config.data_directory, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):

            if record.id in target_genome_all:

                fasta_file_open.write('>%s\n' % record.id)

                clean_sites_seq_split = [record.seq[i:i+n_fna_characters] for i in range(0, len(record.seq), n_fna_characters)]

                for seq in clean_sites_seq_split:
                    fasta_file_open.write('%s\n' % seq)

                fasta_file_open.write('\n')



    fasta_file_open.close()




#build_votu_fasta()



def blast(input_fasta_filename, allowed_contigs, output_directory, single_file=True, max_n=100):

    # NOTE: NEED TO DEFINE A SET OF CONTIG IDS
    
    if input_fasta_filename.endswith('.gz'):
        input_file = gzip.open(input_fasta_filename,"rt")
    else:
        input_file = open(input_fasta_filename,"r")
    
    if single_file:
        output_file = open(output_filename,"w")
    else:
        os.system('mkdir -p %s' % output_directory)
    
    num_recorded = 0
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id in allowed_contigs:
            #print(record.id)
            if single_file:
                SeqIO.write(record, output_file, "fasta")
            else:
                output_filename = "%s/%s.fasta" % (output_directory,record.id)
            
            output_file = open(output_filename,"w")
            SeqIO.write(record, output_file, "fasta")
            output_file.close()
            num_recorded+=1
            
            if (max_n) > 0 and (num_recorded>=max_n):
                break
    
    output_file.close()
    input_file.close()




class classFASTA:

    # class to load FASTA file

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.frn') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list


