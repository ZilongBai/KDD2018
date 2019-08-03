function [Ms] = updateMs(Xs,XfT,F,Ms,gamma,epsilon)
%
% Implements a multiplicative update rule for mixing matrix [Ms] for the structural pathways between blocks, i.e, equation (7) in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis
% In this function we achieve asymmetric dependency by ruling out the impact of M^l_b on Ms.
%
%% Input
%       Xs: n x n adjacency matrix for the structural graph.
%       XfT: n x n x L tensor of behavioral multigraph. Each frontal slice XfT(:,:,j) is the affinity/similarity matrix of one behavioral graph. This corresponds to X^l_b in our paper.
%       F: n x k, the cluster indicator matrix.
%       Ms: n x n, the mixing matrix for the structural pathways between blocks. 
%       gamma: nonnegative scalar, denoting the relative weight for the sparsity regularization on the mixing matrix [Ms] for the structural pathways between blocks.
%       epsilon: an ignorable positive value to avoid issues in multiplicative update rules caused by zero elements in intermediate updating steps.
%
%% Output
%       Ms: n x n, the mixing matrix for the structural pathways between blocks. 

[k,k1] = size(Ms);

numerator = F'*Xs*F;

FtF = F'*F;

denominator = FtF*Ms*FtF;

denominator = denominator + (ones(size(denominator)).*gamma);

Ms = Ms.*(numerator./denominator);

Ms
