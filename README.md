# fbstab-matlab
A MATLAB implementation of the Proximally Stabilized Fischer-Burmeister (FBstab) quadratic programming solver. A C++ implementation is under development and is available here: https://github.com/dliaomcp/fbstab.

## Instructions
Navigate to the root folder then run setpath.m. The directory is arranged assuming that everything is run from the root.

Unit tests can be found in the tests/ folder, MPC examples can be found under mpc_examples/

Please see the entry point functions for documentation. In MATLAB:
- ```help fbstab_mpc```
- ```help fbstab_dense```
- ```help fbstab_sparse```

## Details
The mathematical details regarding FBstab can be found in the following paper: 
https://arxiv.org/abs/1901.04046

## Sponsors
- National Science Foundation Award # CMMI 1562209
- Toyota Research Institute 
