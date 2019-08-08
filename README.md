Code by Derek Lim accompanying paper "The Doubly Stochastic Single Eigenvalue Problem: A Computational Approach" by Derek Lim, Charles R. Johnson.

Contact Derek Lim at dl772@cornell.edu or at cptq on Github with any inquiries.

## Make
`make init` uses pip to install requirements as found in `requirements.txt`

`make test` runs the test suite

## Scripts
The functions to compute DSn or search for exceptions to the Perfect-Mirsky conjecture are included in `permhull.py`.

`computeInequivPairs.g` contains the basic GAP script to compute the inequivalents pairs and output in a `.txt` format. `inequivpairs.py` contains scripts to convert the inequivalent pairs into formats that the python scripts can use to compute with.

## Plots 
Some plots are included in Plots/

`genplots.py` contains code to generate some of the plots from the paper.

To generate plots, import genplots, then call the corresponding plot function.

`data` contains some required data for generating plots.
