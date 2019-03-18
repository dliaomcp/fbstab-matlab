# fbstab-matlab
A MATLAB implementation of the Proximally Stabilized Fischer-Burmeister (FBstab) quadratic programming solver. 

## Instructions
Navigate to the root folder then run setpath.m. The directory is arranged assuming that evening is run from the root.

Unit tests can be found in the tests/ folder, MPC examples can be found under mpc_examples/

Please see the entry point functions fbstab_dense and fbstab_mpc for documentation. In MATLAB:
'''help fbstab_mpc'''
'''help fbstab_dense'''

## Details
The mathematical details regarding FBstab can be found in the following paper: 
https://arxiv.org/abs/1901.04046

## Sponsors
- National Science Foundation Award # CMMI 1562209
- Toyota Research Institute 
