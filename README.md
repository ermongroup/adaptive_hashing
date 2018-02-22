# F2

F2 is a fast and flexible probabilistic approximate model counter for
CNF formulas. It implements the algorithms in our SAT-18 paper.

The code is currently considered beta quality but it is already pretty usable.

## Requirements

### Operating System
The code has been developed in Linux. It should work in any
Unix-derived system including MacOS.

It has not been tested in any version of Windows.

### Python
The following python packages must be available to the Python 3 interpreter.
Probably the code will work correctly with earlier versions but it has
not been tested.

* networkx >= 2.1
* sympy >= 1.1.1
* scipy >= 1.0.0

### Cryptominisat
An executable of [Cryptominisat](https://github.com/msoos/cryptominisat)
 (version 5) should be in the same directory
as `f2.py` or otherwise it should be specified in command line via the
'--cmsat-exe' option.

### SharpSAT
An executable of [sharpSAT](https://github.com/marcthurley/sharpSAT)
should be in the same directory
as `f2.py` or otherwise it should be specified in command line via the
'--sharpsat-exe' option.

## Usage

### Independent Support and Sampling Sets
Can be specified in the CNF files using lines starting with `c ind`.

### Modes of Operation
* Lower Bound (lb)
* Upper Bound (ub)
* Probabilistic Approximation (appr)


### Command line
    usage: python3 f2.py [options] formula_file.cnf

The acceptable options are:

    positional arguments:
      input_filename        The formula in extended DIMACS form

    optional arguments: (default values in [square brackets])
      -h, --help            show a help message and exit
      --working-dir WORKING_DIR
                            Where to create auxiliary files [current directory]
      --keep-cnf            Keep generated auxiliary formulas? [False]
      --random-seed RANDOM_SEED
                            Initialize the random generator. Negative
                            for random initialization. [42]
      --var-degree VAR_DEGREE
                            Average variable degree of LDPC XORs [12]
      --cmsat-exe CMSAT_EXE
                            The path to a cryptominisat v5 executable
                            [Same dir as f2.py]
      --sharpsat-exe SHARPSAT_EXE
                            The path to a SharpSAT executable
                            [Same dir as f2.py]
      --mode {lb,ub,appr}   Mode of operation. Allowed values: lb, ub, appr
                            for lower bound, upper bound or probabilistic
                            approximation respectively.
      --lower-bound LOWER_BOUND
                            Binary logarithm of lower bound (modes: ub, appr)
                            If not specified it will be calculated.
      --upper-bound UPPER_BOUND
                            Binary logarithm of upper bound (only appr mode).
                            If not specified it will be calculated.
      --error-prob ERROR_PROB
                            Probability of error [0.05]
      --confidence-radius CONFIDENCE_RADIUS
                            Tolerance of error as a ratio of correct value [0.2]
      --max-time MAX_TIME   Maximum system running time (in seconds) [3600]
      --lb-n-iter LB_N_ITER
                            (Optional) Override num of iterations in algorithm 1
      --ub-n-iter UB_N_ITER
                            (Optional) Override num of iterations in algorithm 2

### Examples
* Calculate a rigorous lower bound for `formula.cnf` with probability of
error equal to 10%

        python3 f2.py --mode lb --error-prob 0.1 formula.cnf

* Calculate a rigorous upper bound for `formula.cnf` with probability of
error equal to 5% (which is the default) and using for lower bound the
number 30 (log base2 of the lower bound actually)

        python3 f2.py --mode ub formula.cnf

TODO (add more examples)

### Verbosity
TODO


## Code Repository
F2 can be downloaded from the Github repository at
https://github.com/ptheod/f2

## License
To be written.

## Contact
* Panos Theodoropoulos (p + <7 first letters of last name>)@di.uoa.gr
* Dimitris Achlioptas <last 5 letters of last name>@soe.ucsc.edu

