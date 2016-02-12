Copyright (C) 2012 Marius Nicolae

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


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


