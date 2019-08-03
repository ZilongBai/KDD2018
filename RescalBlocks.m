function [RB] = RescalBlocks(B)
% This function serves to rescale the block indicator matrix by setting the maximum membership within each block to 1 in the current iteration of Algorithm 1 BMSM Discovery. in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis 
% This function should follow the RescalMixingM.m function to guarantee that the overall reconstruction performance of the current iteration is not altered.

L = size(B,2);

RB = zeros(size(B));

for l = 1:L

	RB(:,l) = B(:,l)./max(B(:,l));

end
