import numpy as np
import pandas as pd
from parse import *


def summary_statistics(data_structure):
    """
    function accepts validated data structure and runs summary statistics as specified: 
    per chromosome: 
    number of records found 
    number of records found per strand
    min, max, and average lengths of the records

    params: 
    data_structure as pandas dataframe

    outputs: 
    summary statistics

    returns: 
    none
    """

    # number of records found
    records_found = data_structure[
        'chrom'].value_counts().to_frame(name="records")
    records_found.index.name = 'chrom'

    # number of records found per strand
    strand_records = pd.crosstab(
        index=data_structure['chrom'], columns=data_structure['strand'])

    # min, max, and average lengths of the records
    data_struct_lengths = data_structure
    data_struct_lengths['lengths'] = data_structure[
        'end position'] - data_structure['start position']

    lengths_stats = data_struct_lengths[['chrom', 'lengths']].groupby(
        ['chrom']).aggregate([min, max, np.mean])

    summary = pd.concat([records_found, strand_records, lengths_stats], axis=1)

    print(summary)


def main():
    data_25 = read_file("data_samples.25rows.bed.txt")

    print('\n' + "4. Report simple data summary statistics" + '\n')
    print("For each chromosome loaded,")
    print("calculate the number of records found, the number of records found per strand,")
    print("and the min,  max, and average lengths of the records" + '\n')

    summary_statistics(data_25)

# EXAMPLE USAGE
if __name__ == '__main__':
    main()
