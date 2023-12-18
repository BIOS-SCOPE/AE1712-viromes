import gzip
import pandas as pd
import glob
import re
import os


def parse_file(file, sample_type, sample_name):

    data = pd.read_csv(file, sep='\t',
                       header=None,
                       skiprows=1,
                       names=['contig', 'read_count', 'contig_length', 'covered_bases'])
    data['sample'] = sample_name
    data['sample_type'] = sample_type
    data['pct_covered'] = data.covered_bases / data.contig_length *100
    data['reads_per_kb_of_genome'] = data.read_count / data.contig_length *1000
    return data

dataframes = []

for f in glob.glob('data/Biller_COVERM_counts/*.txt.gz'):
    print(f'parsing {f}')
    sample_name = os.path.basename(f).split('_', 1)[0]
    dataframes.append(parse_file(f,'Biller', sample_name))

for f in glob.glob('data/orig_bats_data/cellular/*.txt.gz'):
    print(f'parsing {f}')
    sample_name = os.path.basename(f).split('_C', 1)[0]
    dataframes.append(parse_file(f, 'AE1712', sample_name))
complete_df = pd.concat(dataframes)


complete_df = pd.concat(dataframes)

mapping = pd.read_csv('data/phage_mapping.txt', sep='\t')

biller_metadata = pd.read_csv('data/Biller-metadata.tsv', sep='\t')


final_df = complete_df.merge(mapping).merge(biller_metadata)

final_df.to_csv('data/biller_merged_dataframe.tsv.gz', sep='\t', index=False)
