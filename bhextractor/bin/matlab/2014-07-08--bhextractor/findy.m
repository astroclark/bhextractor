function y=findy(V, b1, d, numPCs)
%   function to reconstruct waveform
%   Arguments
%      V -- PC vectors
%      b1 -- vector of PC coefficients
%      d -- scaling factor for distance, d=1 for 10 kpc
%      numPCs -- number of PCs to use
y=V(:,1:numPCs)*(b1(:,1:numPCs))';
%y=y.*d;

