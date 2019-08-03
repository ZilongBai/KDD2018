function [F, Ms, MfT,errstr,errfun] = Solver(Xs,XfT,alpha,gamma,k,MAXI,epsilon)
%
% Solves the Block Modeling from Structural and Behavioral Data (BMSB) for cluster indicator matrix [F], mixing matrix [Ms] for the structural pathways between blocks, and a set of mixing matrices [MfT] for the behavioral interactions over structural pathways in the set of behavioral graphs XfT.
% I.e., Algorithm 1 BMSB Discovery. in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis
%
% The cluster indicator matrix F is sovled with multiplicative update rules in updateF.m
% The mixing matrix [Ms] for the structural pathways between blocks is solved with mulitplicative update rules in updateMs.m
% The set of mixing matrices [MfT] for the behavioral interactions over structural pathways is solved with multiplicative update rules in updateMfT_sparseMs.m
%% The detailed information and deductions can be found in our paper  "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018, https://dl.acm.org/citation.cfm?doid=3219819.3220080 Please refer to the related work for techniques involved in developing our method.
%
%% Input
%	Xs: n x n adjacency matrix for the structural graph.
%	XfT: n x n x L tensor of behavioral multigraph. Each frontal slice XfT(:,:,j) is the affinity/similarity matrix of one behavioral graph. This corresponds to X^l_b in our paper.
% 	alpha: nonnegative scalar, denoting the relative weight for structural graph reconstruction error.
%	gamma: nonnegative scalar, denoting the relative weight for the sparsity regularization on the mixing matrix [Ms] for the structural pathways between blocks.
%	k: The number of blocks (i.e., clusters of nodes).
%	MAXI: The maximum number of iteractions.
%	epsilon: an ignorable positive value to avoid issues in multiplicative update rules caused by zero elements in intermediate updating steps.
%
%% Output
%	F: n x k, the cluster indicator matrix.
%	Ms: n x n, the mixing matrix for the structural pathways between blocks. 
%	MfT: n x n x L, the mixing matrices for the behavioral interations over structural pathways. This corresponds to M^l_b in our paper.	  
%	errstr: length-MAXI vector indicating the relative reconstruction error for the structural graph after each iteration.
%	errfun: length-MAXI vector indicating the relative reconstruction error for the behavioral graphs after each iteraction.

[n,n1,L] = size(XfT);

%% Initialization

F = rand(n,k);
Ms = rand(k,k); Ms = Ms'*Ms;
MfT = zeros(k,k,L);
for l = 1:L
	Mf = rand(k,k); Mf = Mf'*Mf;
	MfT(:,:,l) = Mf;
end

%% Iterative Alternative Least Squares framework. Terminates after MAXI iterations.

for i = 1:MAXI
		
	F = updateF(Xs,XfT,F,Ms,MfT,alpha,epsilon); % Update F with equation (6)

	Ms = updateMs(Xs,XfT,F,Ms,gamma,epsilon); % Update Ms with equation (7)	

	Msr = RescalMixingM(Ms,F); % Rescale Ms
	Fr = RescalBlocks(F); % Rescale F. The rescaling steps serve to facilitate numerical stability in our solver without changing reconstruction performance in each iteraction.
	FxFr = kron(Fr,Fr);
        FtFr = Fr'*Fr;
        FtFxFtFr = kron(FtFr,FtFr);
	MfT = updateMfT_sparseMs(Xs,XfT,FxFr,FtFxFtFr,Msr,MfT,alpha,gamma,epsilon); % Update MfT with equation (11) then reshape each vectorized m^l_b into matrix form.

	[errstr(i),errfun(i)] = errs(Xs,XfT,F,Ms,MfT); % Record reconstruction errors for each iteraction.

end
