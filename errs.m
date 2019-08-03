function [errstr,errfun] = errs(Xs,XfT,F,Ms,MfT)
% This function serves to compute the relative reconstruction errors of the structural graph and the set of behavioral graphs in the current iteration of Algorithm 1 BMSM Discovery. in "Discovering Models from Structural and Behavioral Brain Imaging Data", KDD 2018. Implemented by Zilong Bai, KDD Lab @ University of California, Davis 

recXs = F*Ms*F';

recXfT = zeros(size(XfT));
L = size(XfT,3);

size(recXfT,3)

for l = 1:L
%	size(F)
%	size(MfT(:,:,l))
	recXfT(:,:,l) = F*(MfT(:,:,l).*Ms)*F';

end

errstr = norm(Xs-recXs,'fro')/norm(Xs,'fro');

diff_func = 0;
ori_func = 0;
for l = 1:L

	diff_func = diff_func + norm(recXfT(:,:,l)-XfT(:,:,l),'fro');
	ori_func = ori_func + norm(XfT(:,:,l),'fro');
	
end

errfun = diff_func/ori_func;

