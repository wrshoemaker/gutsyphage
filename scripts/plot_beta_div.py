
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config



phage_metadata_dict = data_utils.read_phage_metadata(checkv_quality='High-quality', checkv_quality_cumulative=True)

data_utils.votu_dict_to_s_by_s(phage_metadata_dict)


sample_metagenome_dict = data_utils.sample_metagenome_dict()