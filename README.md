# KDD2018 Discovering Models from Structural and Behavioral Brain Imaging Data
This repository contains matlab source code to solve the Block Modeling from Structural and Behavioral Data (BMSB) for cluster indicator matrix [F], mixing matrix [Ms] for the structural pathways between blocks, and a set of mixing matrices [MfT] for the behavioral interactions over structural pathways in the set of behavioral graphs XfT.
% I.e., Algorithm 1 BMSB Discovery. in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis

The details of the model and the solver are described in our paper [Discovering Models from Structural and Behavioral Brain Imaging Data](https://dl.acm.org/citation.cfm?doid=3219819.3220080).
The details of deductions for our method as well as the techniques from related work that are involved in developing our method can also be found in our paper. See the references of our paper for more information on related work.

File(s) in this repository:

Solver.m: Solves the Block Modeling from Structural and Behavioral Data (BMSB) for cluster indicator matrix [F], mixing matrix [Ms] for the structural pathways between blocks, and a set of mixing matrices [MfT] for the behavioral interactions over structural pathways in the set of behavioral graphs XfT.
% I.e., Algorithm 1 BMSB Discovery in our paper.

updateF.m : Implements the multiplicative update rule for updating cluster indicator matrix [F], i.e., equation (6) in our paper.

updateMs.m : Implements a multiplicative update rule for mixing matrix [Ms] for the structural pathways between blocks, i.e, equation (7) in our paper.

updateMfT.m : Implements the multiplicative update rule for updating a set of mixing matrices [MfT] for the behavioral interactions over structural pathways in the set of behavioral graphs XfT. 
I.e., equation (11) in our paper.

RescalMixingM.m : serves to rescale the mixing matrix without changing the overall reconstruction error performance in the current iteration of Algorithm 1 BMSM Discovery in our paper.

RescalBlocks.m : serves to rescale the block indicator matrix by setting the maximum membership within each block to 1 in the current iteration of Algorithm 1 BMSM Discovery in our paper.

errs.m : serves to compute the relative reconstruction errors of the structural graph and the set of behavioral graphs in the current iteration of Algorithm 1 BMSM Discovery in our paper.
