function [MfT] = updateMfT(Xs,XfT,FxF,FtFxFtF,Ms,MfT,alpha,gamma,epsilon)
%
% Implements the multiplicative update rule for updating a set of mixing matrices [MfT] for the behavioral interactions over structural pathways in the set of behavioral graphs XfT.
% I.e., equation (11) and the follow-up matricization step in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis
%
%% Input
%       Xs: n x n adjacency matrix for the structural graph.
%       XfT: n x n x L tensor of behavioral multigraph. Each frontal slice XfT(:,:,j) is the affinity/similarity matrix of one behavioral graph. This corresponds to X^l_b in our paper.
%	FxF, FtFxFtF: precomputed intermediate variables based on F to avoid redundancy in computation.
%       Ms: n x n, the mixing matrix for the structural pathways between blocks. 
%       MfT: n x n x L, the mixing matrices for the behavioral interations over structural pathways. This corresponds to M^l_b in our paper.     
%       alpha: nonnegative scalar, denoting the relative weight for structural graph reconstruction error.
%       gamma: nonnegative scalar, denoting the relative weight for the sparsity regularization on the mixing matrix [Ms] for the structural pathways between blocks.
%       epsilon: an ignorable positive value to avoid issues in multiplicative update rules caused by zero elements in intermediate updating steps.
%
%% Output
%       MfT: n x n x L, the mixing matrices for the behavioral interations over structural pathways. This corresponds to M^l_b in our paper.      

[n,n1,L] = size(XfT);
[k,k1] = size(Ms);

xs = reshape(Xs,n*n,1);
xfm = zeros(n*n,L);

ms = reshape(Ms,k*k,1);
mfm = zeros(k*k,L);

mfmnew = zeros(size(mfm));

for l = 1:L
        xfm(:,l) = reshape(XfT(:,:,l),n*n,1);
        mfm(:,l) = reshape(MfT(:,:,l),k*k,1);
	mfmnew(:,l) = mfm(:,l) .* (( (FxF'*xfm(:,l)).*ms) ./( (FtFxFtF*(ms .* mfm(:,l))).*(ms + epsilon.*ones(size(ms))) ));
%	mfmnew(:,l) = mfm(:,l) .* sqrt( (FxF'*xfm(:,l)) ./ (FxF'*FxF*(ms .* mfm(:,l)))  );
end

for l = 1:L

	MfT(:,:,l) = reshape( mfmnew(:,l) , k , k );

end
