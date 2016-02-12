
Source code for pattern matching with mismatches and pattern matching with 
k mismatches algorithms.
- A naive O(nm)
- Abrahamson's O(n\sqrt{m \log m})
- Subset k-mismatch (O(nk))
- Knapsack k-mismatch (O(n\sqrt{k \log k})) 

To compile:
===========
Make sure you have g++, make and the fftw3 library installed. Open a command line and type:

make

To run:
=======

./kmismatch inputFile 

where inputFile contains the text T. Various patterns of different lengths will be chosen 
from T and various algorithms will be tested.


