
import gzip
import config


repeated_votu_representative_all = ['vOTU-010833', 'vOTU-014541', 'vOTU-013486']
taxon_abbr_dict = {'d':'domain', 'p':'phylum', 'c':'class', 'o':'order', 'f':'family', 'g':'genus', 's':'species'}
taxon_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
taxon_ranks = taxon_ranks[::-1]

str_to_boolean_dict = {'Yes': True, 'No': False}

def read_sample_metadata():

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

        #for line in f:
        #    print(line)

        f.close()




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



def read_phage_metadata(checkv_quality='High-quality', viral_confidence='Confident', viralverify_prediction='Virus', votu_representative=True):

    sra_dict = {}

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
        # line[14] = genomad_virus_hallmarks (?)
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
        # line[26] = trna_count_suppressor
        # line[27] = sra_run (e.g., ERR414592,ERR414311 or SRR5274010)
        # line[28] = sra_sample (e.g., ERS396364 or SRS1992821)
        # line[29] = biosample (e.g., SAMEA2338694)
        # line[30] = country
        # line[31] = latitude
        # line[32] = longitude

        for line in f:
            line = line.strip().split('\t')

            if (line[8] != checkv_quality):
                continue

            if (line[12] != viral_confidence):
                continue

            if (line[16] != viralverify_prediction):
                continue
            
            # make sure it's representative
            if str_to_boolean_dict[line[2]] != votu_representative:
                continue

            if line[28] == 'Null':
                continue

            sra_sample = line[28]

            if sra_sample not in sra_dict:
                sra_dict[sra_sample] = {}

            #if line[1] == 'vOTU-010833':
            #    print(line)

            if line[1] not in sra_dict[sra_sample]:

                sra_dict[sra_sample][line[1]] = {}
        

        f.close()


    return sra_dict






