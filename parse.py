import numpy as np
import pandas as pd


def read_file(in_file):
    """
    function reads input file and parses contents into a data structure

    params:
    in_file: name of input file if local, or full path if otherwise

    outputs:
    none

    returns:
    pandas dataframe containing annotation information
    """
    # import the data
    data = pd.read_csv(in_file, sep='\t', header=None)
    data.columns = ['chrom', 'start position',
                    'end position', 'feature name', 'strand']

    # create return structure
    data_out = data

    # validate the inputs
    # check chrom
    data_out['chrom'] = iter_chrom(data['chrom'])

    # check start and end position together
    data_out[['start position', 'end position']] = iter_pos(
        data[['start position', 'end position']])

    # check feature name
    data_out['feature name'] = iter_feature(data['feature name'])

    # check strand
    iter_strand(data['strand'])

    # dataframe has passed validation, return dataframe
    return(data_out)


def iter_chrom(chrom_series):
    """
    function takes pandas series of chromosomes, checks against specification
    function is a wrapper for the actual test

    params:
    pandas series containing list of chromosome strings eg. "chr6"

    outputs:
    none

    returns:
    pandas series containing list of chromosome ints eg. 6
    """
    chrom_out = chrom_series.apply(parse_chrom)

    return(chrom_out)


def parse_chrom(chrom_entry):
    """
    function takes single chrom entry, checks against specification:
    first three chars should be 'chr'
    remainder of string should be one or two chars specifying chromosome number (1-22)
    after valdating input string, convert to int and return it

    params:
    string containing chromosome description

    outputs:
    none

    returns:
    int containing chromsome number
    """
    # 'chr'
    if (chrom_entry[0:3] == 'chr'):
        # 'X or 'XX'
        chrom_number = chrom_entry[3:]
        if (1 <= len(chrom_number) <= 2):
            # throws ValueError if input is strange like '.9'
            chrom_int = int(chrom_number)
            # evaluates false in padded zero edge case '06'
            if (str(chrom_int) == chrom_number):
                # 1-22
                if (1 <= chrom_int <= 22):
                    return(chrom_int)
                else:
                    raise ValueError(
                        "Chromosome number should have a value of 1 through 22")
            else:
                raise ValueError("first X in chrXX should not be 0")
        else:
            raise ValueError(
                "chrom variable should be in the format 'chrX' or 'chrXX'")
    else:
        raise ValueError("chrom variable should begin with the prefix 'chr'")


def iter_pos(pos_df):
    """
    function takes pandas dataframe of start and end positions, checks against specification
    function is a wrapper for the actual test

    params:
    pandas dataframe containing start and end positions

    outputs:
    none

    returns:
    pandas dataframe containing start and end positions as unsigned 32-bit ints
    """
    pos_df_out = pos_df.apply(parse_pos, axis=1, result_type='expand')

    return(pos_df_out)


def parse_pos(pos_series):
    """
    function takes series or array containing a start and an end position,
    checks against specification:
    positions are integer values from 1 to 2**32
    end position is greater than start position

    params:
    pandas series containing a start and an end position

    outputs:
    none

    returns:
    array containing start and end positions as 64 bit ints

    note: maximum value in spec is 2**32,
            consider using unsigned 32-bit ints if max value was 2*32-1
            possible optimization exists using uint32 anyway, since min value is 1, not 0
    """
    try:
        valid_start = valid_pos(pos_series[0])
    except ValueError:
        print('Start position invalid')
        raise
    try:
        valid_end = valid_pos(pos_series[1])
    except ValueError:
        print('End position invalid')
        raise

    if (valid_end > valid_start):
        return (np.array([valid_start, valid_end]))
    else:
        raise ValueError('end position must be greater than start position')


def valid_pos(pos_entry):
    """
    function takes singular position value and checks against specification:
    positions are integer values from 1 to 2**32

    params:
    singular position value, likely already int64 if input file has no issues,
    but may be alternate dtype if not

    outputs:
    none

    returns:
    valid position values as 64 bit int
    """
    # convert to string before casting to int, protects against decimals in
    # a float (will throw ValueError)
    pos_out = np.int64(str(pos_entry))

    if (0 < pos_out <= 2**32):
        return (pos_out)
    else:
        raise ValueError('Position should have a value of 1 through 2**32')


def iter_feature(feature_series):
    """
    function takes pandas series of features, checks against specification
    function is a wrapper for the actual test

    params:
    pandas series containing list of features

    outputs:
    none

    returns:
    pandas series containing list of valid feature strings
    """
    feature_out = feature_series.apply(valid_feature)

    return (feature_out)


def valid_feature(feature_entry):
    """
    function takes singular feature value and checks against specification:
    string value from character set of: alphanumeric, underscore, hyphen, parentheses

    params:
    singular feature value, likely already a string
    may be a numeric value in specific cases, so iter_feature and valid_feature are not void

    outputs:
    none

    returns:
    valid feature string

    note: pure python implementation, preferred solution requires regex or string libraries
    """
    feature_string = str(feature_entry)

    allowed = set(
        'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ' + '0123456789' + '_-()')

    if (set(feature_string) <= allowed):
        return(feature_string)
    else:
        raise ValueError(
            'feature name should only contain alphanumeric, underscore, hyphen, parentheses')


def iter_strand(strand_series):
    """
    function takes pandas series of strands, checks against specification
    function is a wrapper for the actual test

    params:
    pandas series containing list of strands

    outputs:
    none

    returns:
    none
    """
    strand_series.apply(valid_strand)


def valid_strand(strand_entry):
    """
    function takes singular strand value and checks against specification:
    string value of '-' or '+'

    params:
    singular strand value

    outputs:
    none

    returns:
    none
    """
    if not (strand_entry == "+" or strand_entry == "-"):
        raise ValueError("strand should be a string value of '-' or '+'")


def main():
    print('\n' + "1. Load/parse the data into a structure" + '\n')
    data_25 = read_file("data_samples.25rows.bed.txt")
    print(data_25)

# EXAMPLE USAGE
if __name__ == '__main__':
    main()
