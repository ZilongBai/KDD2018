function [RM] = RescalMixingM(M,B)
% This function serves to rescale the mixing matrix without changing the overall reconstruction error performance in the current iteration of Algorithm 1 BMSM Discovery. in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis .

k = size(B,2);

RM = zeros(size(M));

for i = 1:k
	for j = 1:k
		RM(i,j) = M(i,j)*max(B(:,i))*max(B(:,j));
	end
end
