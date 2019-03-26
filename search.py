import numpy as np
import pandas as pd
from parse import *


def search_by_position(data_structure, chrom_value, *args, search_range=True):
    """
    function takes data structure, chrom value, and optional position(s) inputs
    returns data structure with records matching chrom value and position, if provided

    params:
    data_structure as pandas dataframe
    chrom_value as chromosome string
    args as optional position integers
    search_type with range as default

    outputs:
    none

    returns:
    pandas dataframe containing matching records, or
    None if no matches
    """

    # validate chromosome value
    chrom = parse_chrom(chrom_value)

    # validate positional inputs
    # must do this first, else risk returning (technically correct)
    # None when given invalid positions
    position_inputs = pd.Series(args)
    positions_valid = position_inputs.apply(valid_pos)

    # match records on chromosome
    chrom_match = data_structure.loc[data_structure['chrom'] == chrom]
    if (len(chrom_match) == 0):
        return None

    # match records on positions if given
    if (len(positions_valid) == 0):
        # no positions
        records_out = chrom_match
    elif (len(positions_valid) == 2 and search_range):
        # two positions given, default behavior assumes they define a range
        records_out = pos_range_search(chrom_match, positions_valid)
    else:
        # any number of positions given represents individual locations to be checked
        # this implementation prioritizes efficiency, as the called function is vectorized
        # as opposed to a loop within a loop
        # downside is results may be scrambled relative to original input
        records_list = []
        for pos_indiv in positions_valid:
            records_list.append(indiv_pos_search(chrom_match, pos_indiv))
        # combine together
        records_found = pd.concat(records_list)
        # remove potential duplicates
        records_out = records_found.drop_duplicates()

    if (len(records_out) == 0):
        return None
    else:
        return(records_out)


def pos_range_search(data_structure, position_series):
    """
    function takes data structure, and a series containing two validated positions
    returns data structure with records matching the positional range defined

    assumption: the entire record must exist within the range provided,
    deduced by relative sizes given in example

    note: spec specifies 'start position <= position < end position'
    this is assumed to refer to the singular position search
    range search implementation is 'start position <= position <= end position'
    because otherwise specifying a range exactly equal to an entry in the data
    returns no result, which is illogical

    params:
    data_structure as pandas dataframe
    position_series as series containing two position integers
    Note: does not care about order of positions

    outputs:
    none

    returns:
    pandas dataframe containing matching records, may be empty dataframe
    """
    # check edge case
    if (position_series[0] == position_series[1]):
        # pure implementation, prefer to use warnings.warn()
        print("Warning: Position range given is length 0, defaulting to single position search")
        return (indiv_pos_search(data_structure, position_series[1]))
    else:
        range_start = min(position_series)
        range_end = max(position_series)

        data_struct_check_start = data_structure.loc[
            data_structure['start position'] >= range_start]
        data_structure_out = data_struct_check_start.loc[
            data_struct_check_start['end position'] <= range_end]

        return(data_structure_out)


def indiv_pos_search(data_structure, pos_value):
    """
    function takes data structure and single validated position input
    returns data structure with records matching position spec:
    start position <= position < end position

    params:
    data_structure as pandas dataframe
    pos_value as position integer

    outputs:
    none

    returns:
    pandas dataframe containing matching records, may be empty dataframe
    """
    data_struct_check_start = data_structure.loc[
        data_structure['start position'] <= pos_value]
    data_structure_out = data_struct_check_start.loc[
        data_struct_check_start['end position'] > pos_value]

    return(data_structure_out)


def search_by_feature(data_structure, feature_value):
    """
    function takes data structure and feature string as input
    returns data structure with records that are an exact match of the feature string

    params:
    data_structure as pandas dataframe
    feature_value as string

    outputs:
    none

    returns:
    pandas dataframe containing matching records, or
    None if no matches
    """
    data_structure_out = data_structure.loc[
        data_structure['feature name'] == feature_value]

    if (len(data_structure_out) == 0):
        return None
    else:
        return (data_structure_out)


def main():
    data_25 = read_file("data_samples.25rows.bed.txt")

    print('\n' + "2. Perform search query" + '\n')
    print("select records on chromosome 6 between positions 169000000 and 170000000")
    print(search_by_position(data_25, 'chr6', 169000000, 170000000))

    print('\n' + 'select records with feature name "asdf1234qwer0987"')
    print(search_by_feature(data_25, "asdf1234qwer0987"))

# EXAMPLE USAGE
if __name__ == '__main__':
    main()
