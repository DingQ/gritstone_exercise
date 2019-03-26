import numpy as np
import pandas as pd
import parse
import search
import test
import summary


def wait_enter():
	"""
	function waits for evaluator to press enter before running the next part of the assignment

	params:
	none

	outputs:
	enter prompt

	returns:
	none
	"""
    input('\n' + "Press Enter to continue")


def main():
    parse.main()
    wait_enter()

    search.main()
    wait_enter()

    test.main()
    wait_enter()

    summary.main()


if __name__ == '__main__':
    main()
