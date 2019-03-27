# Gritstone Coding Exercise

[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)

This code repository contiains a solution to the Gritstone Oncology interview coding exercise. 

## Getting Started 

### Prerequisites 

This project was developed in `Python 3.6.5` with the following external libraries: 
```
numpy 1.16.2
pandas 0.24.2
```
Also tested in a second environment: 
```
python 3.6.6
numpy 1.14.5
pandas 0.24.2
```

### Installation

Clone this repo to your local machine from `https://github.com/DingQ/gritstone_exercise`

### Usage

Open a terminal in the local repo directory and run: 
```
python quiz_run_all.py
```

Individual components (the four files discussed in the next section) may also be run independently.


## Discussion of solution

Code for the four parts of the assignment are contained in four files: 
```
parse.py
search.py
test.py
summary.py
```

### parse.py

Individual functions validate and parse inputs according to the specification. 
Although the process of iterating over the data with these functions is only one line in this implementation, this is abstracted into individual methods to enable easier debugging of improperly formed inputs. 

### search.py 

`search_by_position()` was an interesting function, which required different behaviors based on a varying number of positional inputs. If no positions are provided, the function returns all records that match the chromosome given. 

If position(s) are given, according to the second bullet point under 'search by position' in the spec, the method returns matching records that occur on the relevant chromosome that match the pattern: 
```
start position <= given position(s) < end position
```
This functionality is preserved regardless of the number of position arguments, EXCEPT: If exactly two positions are given, according to the first bullet point under 'perform search query' in the spec, the method returns matching records that occur on the relevant chromosome that match the pattern: 
```
smaller given position <= start position <= end position <= larger given position
```
As this is the desired default behavior, it is executed via a function call that resembles: 
```python
search_by_position(data_structure, chromosome_number, position1, position2)
```
If the user wishes to search each of the two provided positions individually, an additional keyword parameter `search_range` must be set to False (default value is True): 
```python
search_by_position(data_structure, chromosome_number, position1, position2, search_range = False)
```
The two behavior modes are implemented using two additional individual functions. 

In contrast, `search_by_feature()` was relatively straightforward. 

### test.py

test.py started out as a series of assert statements verifying correct answer outputs from functions given valid inputs. However, the spec requires raising specific exceptions when provided incorrect inputs. Testing these without the use of `unittest.assertRaises()` involves using try/except blocks. These were bundled into a helper function `invalid_test()`. For consistency, a similar helper function for correct result tests, `valid_test()`, was also written and used. 

These helper functions were used to implement unit tests for `search_by_position()` and all of the functions it depends on, resulting in tests for: 
 - chromosome number parsing
 - position number valdiation
 - searching position by range 
 - searching position by single location
 - the aggregate search_by_position function

### summary.py

These summary statistics were computed with the help of some `pandas` aggregate functions. 

### Miscellaneous notes 

An object-oriented solution was briefly considered early in the development cycle. However, this project is essentially the transformation of an input to an output - it is not particularly interactive, and the specification did not request or mention a GUI or API. As part of a larger system, creating an object to hold the data and carry associated methods might make sense, but with no indication of this bigger picture, such an optimization was not pursued at this time. 

