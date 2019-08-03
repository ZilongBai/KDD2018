function [F] = updateF(Xs,XfT,F,Ms,MfT,alpha,epsilon)
%
% Implements the multiplicative update rule for updating cluster indicator matrix [F], i.e., equation (6) in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis
%
% The cluster indicator matrix F is sovled with multiplicative update rules in updateF.m
% The mixing matrix [Ms] for the structural pathways between blocks is solved with mulitplicative update rules in updateMs.m
% The set of mixing matrices [MfT] for the behavioral interactions over structural pathways is solved with multiplicative update rules in updateMfT_sparseMs.m
%% The detailed information and deductions can be found in our paper  "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018, https://dl.acm.org/citation.cfm?doid=3219819.3220080 Please refer to the related work for techniques involved in developing our method.
%
%% Input
%       Xs: n x n adjacency matrix for the structural graph.
%       XfT: n x n x L tensor of behavioral multigraph. Each frontal slice XfT(:,:,j) is the affinity/similarity matrix of one behavioral graph. This corresponds to X^l_b in our paper.
%       F: n x k, the cluster indicator matrix before this interation.
%       Ms: n x n, the mixing matrix for the structural pathways between blocks before this iteration.
%       MfT: n x n x L, the mixing matrices for the behavioral interations over structural pathways before this iteration. 
%       alpha: nonnegative scalar, denoting the relative weight for structural graph reconstruction error.
%       epsilon: an ignorable positive value to avoid issues in multiplicative update rules caused by zero elements in intermediate updating steps.
%
%% Output
%       F: n x k, the cluster indicator matrix updated in the current iteration.

[n,n1,L] = size(XfT);
[k,k1] = size(Ms);

P = (Xs'*F*Ms).*alpha;

for l = 1:L

	P = P + XfT(:,:,l)'*F*(MfT(:,:,l).*Ms);

end

denominator = F*(F'*P);
denominator(find(denominator<=0)) = epsilon;

SQF = P./denominator;
SQF(find(SQF<=0)) = epsilon;

Fnew = F.*sqrt(SQF);

F = Fnew;
