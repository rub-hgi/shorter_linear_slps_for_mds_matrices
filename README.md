# Shorter Linear SLPs for MDS Matrices

SLPs for our ToSC Volume 2017 Issue 4 paper 'Shorter Linear Straight-Line Programs for MDS Matrices'

- The makefile compiles the slp_heuristic C++ program.
- The makefile compiles for each matrix header in the 'paar_header' folder
  two executables: one for Paar1 algorithm and one for Paar2. Each executable prints
  the XOR count for this matrix and the corresponding SLP for implementing the matrix
  with the said XOR count.
- For each matrix there is a second file in the 'bp_format' folder
  that contains the correct format for the slp_heuristic executable. You can
  call the programs on, e.g. the M_4_4 matrix, as follows:
  $ ./slp_heuristic < bp_format/M_4_4.txt
  It will then output the corresponding XOR count and the SLP.
- All matrices and functions for converting them etc. are additionally in the code.sage
  file. This was also used to construct the above headers and text files. The SageMath
  code can be loaded into a SageMath REPL and then used to compute Paar1 or slp_heuristic
  XOR counts (for the latter the program need to be available as a compiled executable).

- The code in the slp_heuristic.cpp is based on the implementation of Boyar and Peralta,
  available at http://www.imada.sdu.dk/~joan/xor/. We adjusted the program such that it
  prints the corresponding SLP program.
- Three SLPs where provided by Andrea Visconti through personal communication. The code
  used to generate these SLPs is not available to us.
- Twenty-five SLPs were generated with code kindly provided by René Peralta, based on
  recent work of Joan Boyar, Magnus Gausdal Find, and him. This program is called LinOpt
  in this repository. Its source code is not yet public.

- For every matrix listed in the paper there is a corresponding file in the
  'slp_implementations' directory. These files contain a comment with what program the
  SLP was computed, the XOR count and the SLP itself.
- All SLPs are also contained in the 'check_implementations.sage' file, which can be
  loaded in the SageMath REPL and all SLPs can be verified with the 'run_tests()'
  function.
