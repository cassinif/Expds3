# Expds3

This is a companion software to the manuscript

[F. Cassini. Efficient third-order tensor-oriented directional splitting
for exponential integrators, Electron. Trans. Numer. Anal., 60:520--540, 2024](https://doi.org/10.1553/etna_vol60s520)

In particular, it contains the scripts needed to 
reproduce the MATLAB numerical examples of the article.
The C++/CUDA code is available upon request to the author.

For the execution, some external functions are required. They can be found,
for instance, in the following repositories:

* [KronPACK](https://github.com/caliarim/KronPACK)
* [Phipm_simul_iom](https://github.com/drreynolds/Phipm_simul_iom)

Notice that the function ```phipm_simul_iom``` has been 
slightly modified in order to work with actions of 
matrices instead of matrices themselves.
