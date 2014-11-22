function y=findy(V, b1, d, mass, numPCs)
% V_ft, betanew, dnew, mnew, numPCs
%   function to reconstruct waveform
%   Arguments
%      V -- PC vectors
%      b1 -- vector of PC coefficients
%      d -- scaling factor for distance, d=1 for 10 kpc
%      numPCs -- number of PCs to use
y=V(:,1:numPCs)*(b1(:,1:numPCs))';

% ADD MASS OPERATIONS TO MODEL HERE
%y=y.*(mass/250);

scale_ratio = mass/250;
f_vec = (0:1:length(y)-1)';

wave_unscaled=y.*(scale_ratio);

y_mass_scaled = interp1(f_vec,wave_unscaled,f_vec*scale_ratio,'nearest');

% y_unscaled = (abs(y)).^2;
% y_scaled = (abs(y_mass_scaled)).^2;
% figure(3);
% semilogy(f_vec,y_unscaled,f_vec,y_scaled);
% xlim([0,200]);
% ylim([10^(-50),10^(-39)]);
% mleg=legend('M=250',sprintf('M=%i (scaled)',round(mass)));
% set(mleg,'Location','SouthWest')

y=y_mass_scaled;

%y=y.*d;