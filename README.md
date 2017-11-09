# Shorter Linear SLPs for MDS Matrices

SLPs for our ToSC Volume 2017 Issue 3 paper 'Shorter Linear Straight-Line Programs for MDS Matrices'

- The makefile compiles the slp_heuristic C++ program.
- The makefile compiles the LinOpt C++ program.
- The makefile compiles for each matrix header in the 'rebuttal_paar_header' folder
  two executables: one for Paar1 algorithm and one for Paar2. Each executable prints
  the xor count for this matrix and the corresponding SLP for implementing the matrix
  with the said xor count.
- For each matrix there is a second file in the 'rebuttal_bp_format' folder
  that contains the correct format for the slp_heuristic and LinOpt executables. You can
  call the programs on, e.g. the M_4_4 matrix, as follows:
  $ ./slp_heuristic < rebuttal_bp_format/M_4_4.txt
  $ ./LinOpt rebuttal_bp_format/M_4_4.txt 20000 $RANDOM
  where $RANDOM is a seed for the internal random number generator.
  It will then output the corresponding xor count and the SLP.
- The LinOpt program is nondeterministic and actually prints many attempts, not only the
  best result. One could retrieve the lowest xor count by issuing, e.g.:
  $ grep count <resultfile> | sort -n | head -n1
- All matrices and functions for converting them etc. are additionally in the code.sage
  file. This was also used to construct the above headers and text files. The sage code
  can be loaded into a sage REPL and then used to compute Paar1, slp_heuristic or linopt
  xor counts (for the latter two the programs need to be available as a compiled
  executable).

- The code in the slp_heuristic.cpp is based on the implementation of Boyar and Peralta,
  available at http://www.imada.sdu.dk/~joan/xor/. We adjusted the program such that it
  prints the corresponding SLP program.
- The LinOpt.cpp code is not public and was provided by Boyar and Peralta through personal
  communication. We have permission to include it here, but we ask the reviewers to not
  share it elsewhere. We adjusted it such that its input format matches the one required
  for slp_heuristic.
- Three SLPs where provided by Andrea Visconti through personal communication. The code
  used to generate these SLPs is not yet available for us.

- For every matrix listed in the paper there is a corresponding file in the
  'slp_implementations' directory. These files contain a comment with what program the
  SLP was computed, the XOR Count and the SLP itself.
- All SLPs are also contained in the 'check_implementations.sage' file, which can be
  loaded in the sage REPL and all SLPs can be verified with the 'run_tests()' function.
