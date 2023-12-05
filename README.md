# testnh
Suite of programs for the analysis of non-homogeneous substitution processes

## What does TestNH do?

TestNH is a package dedicated to the testing of non-homogeneous process in sequence evolution.

## What does it provide?

It currently contains the following programs:

- *testnh* implements the Bowker test for sequence non-stationarity, as described in (3)
- *mapnh* performs substitution mapping and cluster branches according to their underlying substitution processes, as described in (1) and (2)
- *partnh* fits a non-homogeneous model of evolution according to branch partitions, as defined from a clustering tree. It can test different sets of partitions and use a model selection criterion to select the appropriate number of clusters as described in (2)
- *randnh* generates random non-homogeneous models using two models of non-homogeneity, corresponding to the clustering algorithms implemented in mapnh.

## How can i get it?

The TestNH programs are command-line driven. The programs depend on
the Bio++ libraries. 

<!-- Pre-compiled executables are statically linked
(see http://bioweb.me/testnh/releases) and therefore already include
all required code from the libraries.


You can get executable files pre-compiled for your system (if there
are any), use pre-compiled packages (if there are any) or compile the
programs yourself (should work on any system with a decent C++
compiler). Pre-compiled packages will ask for all required
dependencies, which can be found in the same download directory. The
needed libraries are [bpp-core](https://github.com/BioPP/bpp-core),
[bpp-seq](https://github.com/BioPP/bpp-seq),
[bpp-phyl](https://github.com/BioPP/bpp-phyl).
-->


For compiling the programs yourself, from the downloaded sources or
from the git repository, please follow the instructions from [bppsuite]
(https://github.com/BioPP/bppsuite#readme) and up to date [mapnh](https://pbil.univ-lyon1.fr/bpp-doc/testnh/testnh.html).



## How do I use it?

Several example data sets are distributed along with the source code
of the package. A reference manual is also available at
http://bioweb.me/testnh/releases.

## References

- (1) Dutheil JY, Galtier N, Romiguier J, Douzery EJ, Ranwez V, Boussau B. Efficient selection of branch-specific models of sequence evolution. Mol Biol Evol. 2012 Jul;29(7):1861-74.
- (2) Romiguier J, Figuet E, Galtier N, Douzery EJ, Boussau B, Dutheil JY, Ranwez V. Fast and robust characterization of time-heterogeneous sequence evolutionary processes using substitution mapping. PLoS One. 2012;7(3):e33852.
- (3) Dutheil J, Boussau B. Non-homogeneous models of sequence evolution in the Bio++ suite of libraries and programs. BMC Evol Biol. 2008 Sep 22;8:255.
- (4) Guéguen L., Duret L., Unbiased estimate of synonymous and non-synonymous substitution rates with non-stationary base composition, Mol Biol Evol. 2017.
