function logL = like_gauss_tspace(wave_ft, noise, sigma2, deltaF, len, freqs, ...
lowfreq_index, highfreq_index, T_shift, dis, ra, dec, psi, tmp, detector, model,varargin)
% This function calculates the likelihood, assuming a Gaussian
% distribution, given some data (containing a signal and noise), the
% PSD of the noise, and a model function.
%
% The code will return the likelihood value, and the natural logarithm of
% the likelihood value.

% check whether model is a function handle (or string that is a function
% name)
if ischar(model)
    f = str2func(model);
elseif isa(model, 'function_handle')
    f = model;
else
    error('Error... Expecting a model function!');
end
td=timedelay(ra,dec,tmp,detector);
%%tshift=(exp(-2*pi*1i*freqs*td));
tshift = 1;
%res=2*(real(wave_ft).*real(noise) + imag(wave_ft).*imag(noise));
wave_ft=wave_ft.*tshift;
%wave_ft=wave_ft*(1/4096);
%%wave_ft=(wave_ft+noise);

% evaluate the model, including a time shift
%md_ft = abs(feval(f, varargin{:}).*tshift);
Fp = det_response_wrapper(ra, dec,detector,tmp, psi);
%Fp=1;
md_ft = feval(f, varargin{:});
%md_ft=fft(md_ft);
%md_ft=md_ft(2:6143+1,:)*(1/4096);
md_ft=md_ft*(Fp);% +res;
%md_ft=md_ft*4*deltaF;
% check that the length of the data and the model are the same
if length(wave_ft) ~= length(md_ft)
    disp(sprintf('length(wave_ft)=%d',length(wave_ft)))
    disp(sprintf('length(md_ft)=%d',length(md_ft)))
    error('Error... Length of data and model are not the same!');
end

%sigma=sigma2(1:end,:);
sigma=sigma2(lowfreq_index:highfreq_index,:);

%D=wave_ft(1:end,:).^2;
%M=md_ft(1:end,:).^2;
%logterm2=-(D+M)./(2*sigma);

% works out log of the bessel term in likelihood
%x=sqrt((D.*M))./sigma;
%I=besseli(0,x);
%logI=log(I);

% For some waveforms logI can have points equal to infinity. For these
% points we use an approximation of the bessel function and replace
% infinities
%  F=find(logI==max(logI));
%   if logI(F(1))==Inf
%      logterm3=-log(sqrt(2*pi*x));
%   logterm4=x;
%   logplu=logterm3+logterm4;
%   for i=1:length(F);
%      logI(F(i))=logplu(F(i));
%   end
%  end
%save Lterms x logterm1 logterm2 logterm3 logterm4 logI D P sigma2

% get the log likelihood
%logL=-sum(logterm2+logI);

%duration=length(wave_ft(:,1))*deltaT;
%deltaF=1/duration;
% get the log likelihood

%logL = -2*deltaF*sum(((abs(wave_ft(1:end) - md_ft(1:end))).^2)./(sigma));
logL = -2*deltaF*sum(((abs(wave_ft(lowfreq_index:highfreq_index) - md_ft(lowfreq_index:highfreq_index))).^2)./(sigma));

%logL = -sum(abs((wave_ft(1:end) - md_ft(1:end)).^2)./(2*sigma));
%logL = -2*((1/4096)/len)*sum(((abs(wave_ft(1:end) - md_ft(1:end))).^2)./(sigma));
%logL = -sum(((abs(wave_ft - md_ft)).^2)./(abs(noise).^2));
%logL = -sum(((abs(wave_ft - md_ft)).^2)./noise_PSD); ...
    %- 0.5*sum(log(2*deltaT/(pi*len*Pf1)));

% get likelihood
%L = exp(logL);
