import numpy as np
import pandas as pd
from parse import *
from search import *


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


def main():
    print('\n' + "3. Test code from part 2" + '\n')
    if(search_by_position_test()):
        print('\n' + "All tests...PASSED")
    else:
        print('\n' + 'FAILED')

# EXAMPLE USAGE
if __name__ == '__main__':
    main()
