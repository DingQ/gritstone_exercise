"""
coding_quiz.py
parses a data file and performs search
Michael Ding 03/2019
"""
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


def search_by_position_test():
    """
    function runs unit tests for search_by_position() and subroutines

    params:
    none

    outputs:
    test results

    returns:
    boolean reporting whether function passed all tests
    """
    # test dependent subroutines
    tests = []

    tests.append(parse_chrom_test())
    tests.append(valid_pos_test())
    tests.append(pos_range_search_test())
    tests.append(indiv_pos_search_test())

    # test whole function
    test_data = {'chrom': [6, 9, 6],
                 'start position': [169717822, 136953847, 154931457],
                 'end position': [169717906, 136954255, 174931577],
                 'feature name': ["470_368746_55274(PHF10)_4", "464_319955_286256(LCN12)_5", "462_296184_25929(GEMIN5)_5"],
                 'strand': ["-", "+", "-"]}
    test_df = pd.DataFrame(test_data)

    tests.append(valid_test("chrom value only search match",
                            search_by_position, test_df.loc[[0, 2]], test_df, "chr6"))
    tests.append(valid_test("chrom value only search no match",
                            search_by_position, None, test_df, "chr10"))
    tests.append(valid_test("single position match", search_by_position,
                            test_df.loc[[1]], test_df, "chr9", 136953850))
    tests.append(valid_test("single position no match",
                            search_by_position, None, test_df, "chr9", 500))
    tests.append(valid_test("two position match", search_by_position, test_df.loc[
                 [0, 2]], test_df, "chr6", 1549314560, 169717850, search_range=False))
    tests.append(valid_test("two position no match", search_by_position,
                            None, test_df, "chr9", 300, 800, search_range=False))
    tests.append(valid_test("range search no matches",
                            search_by_position, None, test_df, "chr9", 300, 800))
    tests.append(valid_test("range search single match", search_by_position,
                            test_df.loc[[0]], test_df, "chr6", 169717822, 169717906))
    tests.append(valid_test("range search two matches", search_by_position,
                            test_df.loc[[0, 2]], test_df, "chr6", 152931400, 175000000))

    return(all(tests))


def valid_test(test_name, method, answer, *args, **kwargs):
    """
    helper function for running a test that should work

    params:
    test_name as string describing test
    method as function being tested
    answer as correct response
    *args as arguments to be passed to method()

    outputs:
    test results

    returns:
    boolean reporting whether the function passed or not
    """
    print(test_name, end="...")

    try:
        if (isinstance(answer, pd.DataFrame)):
            assert(answer.equals(method(*args, **kwargs)))
        else:
            assert(method(*args, **kwargs) == answer)
        print("PASSED")
        return True
    except AssertionError:
        print("FAILED")
        return False


def invalid_test(test_name, method, *args):
    """
    helper function for running a test that should throw a ValueError

    params:
    test_name as string describing test
    method as function being tested
    *args as arguments for function being tested

    outputs:
    test results

    returns:
    boolean reporting whether the function passed or not
    """
    print(test_name, end="...")
    try:
        method(*args)
        print("FAILED")
        return False
    except ValueError:
        print("PASSED")
        return True


def parse_chrom_test():
    """
    function runs unit tests for parse_chrom()

    params:
    none

    outputs:
    test results

    returns:
    boolean reporting whether function passed all tests
    """
    tests = []

    tests.append(valid_test(
        "min valid chromosome number", parse_chrom, 1, "chr1"))
    tests.append(valid_test("max valid chomosome number",
                            parse_chrom, 22, "chr22"))
    tests.append(invalid_test(
        "noninteger chromosome number", parse_chrom, "chr.3"))
    tests.append(invalid_test(
        "padded zero chromosome number", parse_chrom, "chr06"))
    tests.append(invalid_test(
        "greater than max chromosome number", parse_chrom, "chr25"))
    tests.append(invalid_test(
        "fewer than min chromosome number", parse_chrom, "chr0"))
    tests.append(invalid_test(
        "missing prefix chromosome number", parse_chrom, "5"))

    return(all(tests))


def valid_pos_test():
    """
    function runs unit tests for valid_pos()

    params:
    none

    outputs:
    test results

    returns:
    boolean reporting whether function passed all tests
    """
    tests = []

    tests.append(valid_test("min valid position number", valid_pos, 1, 1))
    tests.append(valid_test("max valid position number",
                            valid_pos, 4294967296, "4294967296"))
    tests.append(invalid_test(
        "fewer than min invalid position number", valid_pos, "0"))
    tests.append(invalid_test(
        "greater than max invalid position number", valid_pos, 4294967297))
    tests.append(invalid_test("not a number position number",
                              valid_pos, "asdf1234qwer0987"))
    tests.append(invalid_test("noninteger position number", valid_pos, 2.53))

    return(all(tests))


def pos_range_search_test():
    """
    function runs unit tests for pos_range_search()

    params:
    none

    outputs:
    test results

    returns:
    boolean reporting whether function passed all tests
    """
    test_data = {'chrom': [4, 4, 4],
                 'start position': [169717822, 136953847, 154931457],
                 'end position': [169717906, 136954255, 154931577],
                 'feature name': ["470_368746_55274(PHF10)_4", "464_319955_286256(LCN12)_5", "462_296184_25929(GEMIN5)_5"],
                 'strand': ["-", "+", "-"]}
    test_df = pd.DataFrame(test_data)

    tests = []

    tests.append(valid_test("single result range search", pos_range_search,
                            test_df.loc[[0]], test_df, pd.Series([160000000, 170000000])))
    tests.append(valid_test("two results range search", pos_range_search, test_df.loc[
                 [0, 2]], test_df, pd.Series([150000000, 170000000])))
    tests.append(valid_test("exact range search", pos_range_search, test_df.loc[
                 [1]], test_df, pd.Series([136953847, 136954255])))
    tests.append(valid_test("overlap range search", pos_range_search, test_df.iloc[
                 0:0], test_df, pd.Series([169717850, 170000000])))
    tests.append(valid_test("nonoverlap range search", pos_range_search,
                            test_df.iloc[0:0], test_df, pd.Series([3, 500])))

    return(all(tests))


def indiv_pos_search_test():
    """
    funciton runs unit tests for indiv_pos_search()

    params:
    none

    outputs:
    test results

    returns:
    boolean reporting whether function passed all tests
    """
    test_data = {'chrom': [4, 4, 4],
                 'start position': [169717822, 136953847, 169717821],
                 'end position': [169717906, 136954255, 169717907],
                 'feature name': ["470_368746_55274(PHF10)_4", "464_319955_286256(LCN12)_5", "462_296184_25929(GEMIN5)_5"],
                 'strand': ["-", "+", "-"]}
    test_df = pd.DataFrame(test_data)

    tests = []

    tests.append(valid_test("single result indiv search",
                            indiv_pos_search, test_df.loc[[1]], test_df, 136953850))
    tests.append(valid_test("two result indiv search",
                            indiv_pos_search, test_df.loc[[0, 2]], test_df, 169717850))
    tests.append(valid_test("no result indiv search",
                            indiv_pos_search, test_df.iloc[0:0], test_df, 300))

    return(all(tests))


# EXAMPLE USAGE
if __name__ == '__main__':
        # 1. Load/parse the data into a structure
    print('\n' + "1. Load/parse the data into a structure")
    data_25 = read_file("data_samples.25rows.bed.txt")
    print(data_25)

    #input('\n' + "Press Enter to continue")

    # 2. Perform search query
    print('\n' + "2. Perform search query" + '\n')
    print("select records on chromosome 6 between positions 169000000 and 170000000")
    print(search_by_position(data_25, 'chr6', 169000000, 170000000))

    print('\n' + 'select records with feature name "asdf1234qwer0987"')
    print(search_by_feature(data_25, "asdf1234qwer0987"))

    #input('\n' + "Press Enter to continue")

    # 3. Test code from part 2
    print('\n' + "3. Test code from part 2" + '\n')
    if(search_by_position_test()):
        print('\n' + "All tests...PASSED")
    else:
        print('\n' + 'FAILED')

    #input('\n' + "Press Enter to continue")

    # 4. Report simple data summary statistics
    print('\n' + "4. Report simple data summary statistics" + '\n')
    print("For each chromosome loaded, calculate the number of records found")
    print(data_25['chrom'].value_counts())

    #input('\n' + "Press Enter to continue")

    print('\n' + "For each chromosome loaded, calculate the number of records found per strand")
    print(pd.crosstab(index=data_25['chrom'], columns=data_25['strand']))

    #input('\n' + "Press Enter to continue")

    print('\n' + "For each chromosome loaded, calculate the min, max, and average lengths of the records")
    data_tmp = data_25
    data_tmp['lengths'] = data_tmp['end position'] - data_tmp['start position']
    print(data_tmp)

    print(data_tmp[['chrom', 'lengths']].groupby(
        ['chrom']).aggregate([np.min, np.max, np.mean]))
