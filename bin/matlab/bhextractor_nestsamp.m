function bhextractor_nestsamp(run_name, noiseType, model, catalogue, ...
    wv, lowfreq, highfreq, seed, detno,numPCs,doPlot,...
    typeofscaling, scaling, doDistance, doTimeshift)

maxNumCompThreads(1)

BHEX_PREFIX=getenv('BHEX_PREFIX');
FRGETVECT_PATH=getenv('FRGETVECT_PATH');
addpath(genpath(BHEX_PREFIX))
addpath(genpath(FRGETVECT_PATH))

%nested sampling code with coherent analysis of multiple detectors,
%coloured noise, and the possibility of a time shift between the location
%of signal in the data and the reconstruction of it

% a nested sampling code used to compare two supernova models. It does this
% by calculated a Bayes Factor for a waveform using the equation
% y=V*beta+epsilon, over data just containing noise.
% V is matrix containing the principal components and
% epsilon is a gaussian noise vector. beta is a vector of principal
% component coefficients and is what the nested sample is set up to find.
% This is repeated for 2 sets of principal components and so a Bayes Factor
% comparing the two models can be calculated.

% This is set up to compare waveforms from Dimmelmeier et al. 2008 and
% Abdikamalov et al. 2010 but can be adapted to use whatever catalogue you
% want. To create a set of principal components for a catalogue use
% pca_sne_waveform1.m after first using resample.m to resample the
% waveforms in the catalogue to a rate of 4096 Hz.

% This code uses 6 PC's but can be reduced/increased to increase
% speed/accuracy

% Have included functions to the antenna response factor for
% detectors of user's choosing.

% Arguments:
%           run_name -- name of folder in which the store results
%           noiseType -- 'aligo', 'ligo3'. 
%           model -- model used to reconstruct waveform
%               = 'Q', 'HR', 'RO3'
%           catalogue -- catalogue from which to draw waveform
%               = 'Q', 'HR','RO3'
%           wv -- specific waveform to use
%               = 1...13 for Q
%               = 1...15 for HR
%               = 1...20  for RO3
%           lowfreq -- low frequency cut-off
%           highfreq -- high frequency cut-off
%           seed -- random seed used to control noise generation etc.
%           typeofscaling -- type of scaling to be applied to waveform 
%               = 'SNR'  for SNR scaling
%               = 'd' 'dist' 'dis' 'distance'  for Distance scaling
%               = all else  for no scaling (default option)
%               note: case insensitive
%           scaling -- SNR scaling: waveform is scaled to have this SNR;
%                      Distance scaling: waveform is scaled to be at this 
%                      distance, in kpc; default in each case: no scaling
%           detno -- no. of detectors in network; default: 1
%           numPCs -- no. of principal components to use when
%           reconstructing the waveform; default: 6
%               = 1...10 when model ~= Ott
%               = 1...6 when model == Ott
%           doDistance -- flag to use distance as additional parameter;
%           default: not on
%           doTimeshift -- flag to use bounce location in time-series as
%           additional parameter; default: not on
%           doPlot -- flag to plot convergence of algorithm; default: not
%           on

%Normalise inputs
if (ischar(wv)), wv=str2num(wv); end %#ok<*ST2NM>
if (ischar(seed)), seed=str2num(seed); end
if (~exist('typeofscaling','var')), typeofscaling='none'; end
if (~exist('scaling','var')), scaling='';
elseif (ischar(scaling)), scaling=str2num(scaling); end
if (~exist('detno','var')), detno=1;
elseif (ischar(detno)), detno=floor(abs(str2num(detno))); end
if detno < 1 || mod(detno,1)
    fprintf(1,'Error: detno must be positive integer\n');
    detno = floor(abs(detno)+1);
    fprintf(1,'detno changed to %i\n',detno);
end
if (~exist('numPCs','var')), numPCs=6;
elseif (ischar(numPCs)), numPCs=str2num(numPCs); end
if numPCs < 1 || mod(numPCs,1)
    fprintf(1,'Error: numPCs must be positive integer\n');
    numPCs = floor(abs(numPCs)+1);
    fprintf(1,'numPCs changed to %i\n',numPCs);
end
if (~exist('doTimeshift','var')), doTimeshift=0;
elseif (ischar(doTimeshift)), doTimeshift=str2num(doTimeshift); end
if doTimeshift > 0, doTimeshift = 1; else doTimeshift = 0; end
if (~exist('doDistance','var')), doDistance=0;
elseif (ischar(doDistance)), doDistance=str2num(doDistance); end
if doDistance > 0, doDistance = 1; else doDistance = 0; end
if (~exist('doPlot','var')), doPlot=0;
elseif (ischar(doPlot)), doPlot=str2num(doPlot); end

%noise curve to use
if strcmpi(noiseType, 'ligo3')
   noise_curve = 'LIGO3_PSD.txt';
elseif strcmpi(noiseType, 'aligo')
   noise_curve = 'ZERO_DET_high_P.txt';
end

display(sprintf('SMEE: %s %s %s %i', noise_curve, model, catalogue, wv));

reset(RandStream.getDefaultStream,seed);
warning('off','MATLAB:RandStream:GetDefaultStream')
tic

clearvars -except run_name catalogue wv model seed lowfreq highfreq typeofscaling scaling ...
    detno numPCs doTimeshift doDistance doPlot noise_curve BHEX_PREFIX

evnoise = zeros(detno,1);

maxits = 15000;% max no. of iterations expected from nested sampling loop
Z = zeros(maxits,1);
H = zeros(maxits,1);
logw = zeros(maxits,1);
Lw = zeros(maxits,1);
betas = zeros(maxits,numPCs);
Ts = zeros(maxits,1);
Distance = zeros(maxits,1);

%sampling frequency of waveforms
% fs = 4096;
fs = 16384;

%time step to be used
sample_deltaT = 1/fs;

% load the catalogues you want to compare
load(sprintf('%s/data/signal_data/%s_catalogue_theta-0',BHEX_PREFIX,catalogue))

% load the set of eigenvectors for each catalogue
load(sprintf('%s/data/PCA_data/%s_PCs_theta-0',BHEX_PREFIX,model));


% sets up the priors and initial chain values, will need to adjust these to
% include other catalogues. Can use findbetas.m to find max and mins.
maxbeta(1)=50; % maximum value for first coeffiecient
minbeta(1)=-50; % minimum value for first coeffiecient
betarange(1)=maxbeta(1)-minbeta(1); % range of values
betaprior(1)=-log(betarange(1)); % log of the first flat prior

maxbeta(2)=50;
minbeta(2)=-50;
betarange(2)=maxbeta(2)-minbeta(2);
betaprior(2)=-log(betarange(2));

maxbeta(3)=50;
minbeta(3)=-50;
betarange(3)=maxbeta(3)-minbeta(3);
betaprior(3)=-log(betarange(3));

maxbeta(4)=50;
minbeta(4)=-50;
betarange(4)=maxbeta(4)-minbeta(4);
betaprior(4)=-log(betarange(4));

maxbeta(5)=50;
minbeta(5)=-50;
betarange(5)=maxbeta(5)-minbeta(5);
betaprior(5)=-log(betarange(5));

maxbeta(6)=50;
minbeta(6)=-50;
betarange(6)=maxbeta(6)-minbeta(6);
betaprior(6)=-log(betarange(6));

maxbeta(7)=50;
minbeta(7)=-50;
betarange(7)=maxbeta(7)-minbeta(7);
betaprior(7)=-log(betarange(7));

maxbeta(8)=50;
minbeta(8)=-50;
betarange(8)=maxbeta(8)-minbeta(8);
betaprior(8)=-log(betarange(8));

maxbeta(9)=50;
minbeta(9)=-50;
betarange(9)=maxbeta(9)-minbeta(9);
betaprior(9)=-log(betarange(9));

maxbeta(10)=50;
minbeta(10)=-50;
betarange(10)=maxbeta(10)-minbeta(10);
betaprior(10)=-log(betarange(10));

if numPCs > 10;
    for Num=1:(numPCs-10);
        maxbeta(Num+10)=100;
        minbeta(Num+10)=-110;
        betarange(Num+10)=maxbeta(Num+10)-minbeta(Num+10);
        betaprior(Num+10)=-log(betarange(Num+10));
    end
end

%decrease numPCs if neccessary
if numPCs > length(betaprior)
    fprintf(1,'Error: numPCs may be at most 10\n');
    numPCs = length(betaprior);
    betas = zeros(maxits,numPCs);
    fprintf(1,'numPCs changed to %i\n',numPCs);
elseif numPCs > size(PCs_final,2)
    fprintf(1,'Error: numPCs must be less than model catalogue size\n');
    numPCs = size(PCs_final,2)-1;
    betas = zeros(maxits,numPCs);
    fprintf(1,'numPCs changed to %i\n',numPCs);
end

if(doTimeshift)
    maxT=0.01;%seconds
    minT=-0.01;%seconds
    Trange=maxT-minT;%size of time-shift-bin
    Tprior=-log(Trange);

    min_bin = -0.01;%seconds
    max_bin = 0.01;%seconds
    bin_sep = 0.0005;%seconds
    bin_size = Trange;
    mid_T0=[min_bin:bin_sep:max_bin];%positions of mid-points of bins
else
    maxT=0;
    minT=0;
    Trange=maxT-minT;
    Tprior=-log(1);
end

if(doDistance)
    maxd=1;
    %mind=0.05;
    mind=1e-10;
    drange=maxd-mind;
    dprior=-log(drange);
else
    maxd=1;
    mind=1;
    drange=maxd-mind;
    dprior=-log(1);
end


% set the number of active points in the Nested Sampling
numactive = 50;

% set the number of iterations in the MCMC for finding the next active
% point
nits = 10;

% This loads the waveform indicated in the input
wave=MDC_final(:,wv);

% XXX: HACK XXX
%wave=wave.*hann(length(wave));

% Sky position
% theta= pi/2;
% phi = pi/2;
% psi= 45;



%Preparations FFTs
len = length(wave);
duration = len*sample_deltaT;
deltaF = 1/duration;
N=len;

if ( mod(N,2)== 0 )
    numFreqs = N/2 - 1;
else
    numFreqs = (N-1)/2;
end

f = deltaF*[1:1:numFreqs]';

%Generate sigma according to a certain PSD file
noise_struct.deltaT = sample_deltaT;
noise_struct.data = zeros(length(wave),1);

% sig=0;
% for cntsigma=1:500;
% [noise_struct] = detectorNoise(noise_struct,noise_curve);
% sig1 = noise_struct.data(1:numFreqs+1).*conj(noise_struct.data(1:numFreqs+1));
% sig = sig+sig1;
% end
% sigma = sqrt(real(sig) / (cntsigma));
[noise_struct,sigma] = detectorNoise(noise_struct,noise_curve);
for z=1:detno;
    sigdet(:,z)=sigma;
    
    %sigdet(1:numFreqs,z)=sigi(1+1:numFreqs+1,z);
end
for i=1:detno
    %Generate noise to be added to wave
    noise_struct1.deltaT = sample_deltaT;
    noise_struct1.data = zeros(length(wave),1);
    [noise_struct1] = detectorNoise(noise_struct,noise_curve);
    noise = (noise_struct1.data);
    
%     size(noise)
%     length(wave)
    
    if ~isempty(find(strcmpi(typeofscaling,{'SNR'})))
    SNR=computeSNR_colourednoise(wave, sigma, lowfreq, highfreq);
    scale_factor = scaling/SNR;
    fprintf(1, 'Scale Factor: %f\n', scale_factor);
    wave0 = scale_factor*wave;
    effective_distance = 10/scale_factor;
    fprintf(1, 'Effective Distance: %f kpc\n', effective_distance);
    elseif ~isempty(find(strcmpi(typeofscaling,{'distance' 'dis' 'dist' 'd'})))
    scale_factor = 10/scaling;
    wave0 = scale_factor*wave;
    effective_distance = scaling;
    fprintf(1, 'Effective Distance: %f kpc\n', effective_distance);
    else
    wave0 = wave;
    effective_distance = 10;
    scale_factor = 1;
    end
    %wave(:,i) = wave0;
    
    % Matrix of frequencies
    NFFT = length(wave0);

    if iseven(NFFT);
    f = fs/2*linspace(0,1,NFFT/2+1);
    else
    f = fs/2*linspace(0,1,NFFT/2);
    end
    
    % Find index of frequency cut off
    lowfreq_index = find(round(f)==lowfreq,1);
    highfreq_index = find(round(f)==highfreq,1);

    %function to compute SNR of the waveform
    SNR=computeSNR_colourednoise(wave0, sigma, lowfreq,highfreq);
    fprintf(1, 'Detector %i: SNR = %f\n', i, SNR);
    
    % inject waveform into two streams of gaussian noise and compute FFT
    % wave_ft_full(:,i) = noise;
     wave_ft_full(:,i) = fft(wave0)*sample_deltaT+(noise);
     % wave_ft_full(:,i) = fft(wave(:,i))+(noise);
    %since data is real, only want half of DFT
    wave_ft(:,i) = wave_ft_full(1+1:numFreqs+1,i);
    %sigma = sigma(1+1:numFreqs+1,i);
    
    % calculates the log of the evidence for noise only
   %evnoise(i) = -2*(sample_deltaT/len)*sum((abs((wave_ft(:,i))).^2)./(abs((sigma)).^2));
    evnoise(i) = -2*deltaF*sum((abs(wave_ft(lowfreq_index:highfreq_index,i)).^2)./((sigma(lowfreq_index:highfreq_index))));
    fprintf(1, 'Detector %i: log(Noise evidence) = %f\n', i, evnoise(i));
   
    %BEGIN save plots for testing purposes (to check against the Shoemaker aLIGO noise curve)
    %added by PK 12/13/11
%     [fpk, asdpk] = textread('ZERO_DET_high_P.txt', '%f %f');
%     figure(42); clf;
%     noisetest = noise(1+1:numFreqs+1, i);
%     loglog(f, sqrt(sigma));
%     hold on;
%     loglog(f, abs(noisetest), 'g');
%     loglog(fpk, asdpk, 'r');
%     hold off;
%     legend('sigma', 'noise', 'DCC');
%     saveplot('.', 'check_of_variable_noise', 'pdf');
%     xlim([100, 150]);
%     saveplot('.', 'check_of_variable_noise_zoom', 'pdf');
    %END save plots for testing purposes 

 
end


V=PCs_final * 1.1963542944976191e-20 * scale_factor;

%PCs in fourier domain
V_ft_full = fft(V)*sample_deltaT;
V_ft = V_ft_full(1+1:numFreqs+1,:);
%V_ft=real(V_ft);
%wave_ft=real(wave_ft);
% create active points initially drawn from a uniform distribution for
% each value in beta
for cnt=1:numPCs
    activebeta(:,cnt) = rand(numactive,1)*(maxbeta(cnt)-minbeta(cnt)) + minbeta(cnt);
end


actived = rand(1,numactive)*(maxd-mind) + mind;

if(doTimeshift)
    %find sum of likelihood for points in all the time shift bins, and then
    %only use range of bins with highest summed likelihood
    for k=1:length(mid_T0)
        %activeT_range(:,k) = rand(1,numactive)*(maxT-minT) - Trange/2 + mid_T0(k);
        activeT(:,k) = rand(1,numactive)*(maxT-minT) - Trange/2 + mid_T0(k);
        
        like_sum = -1e37;
        % calculate the likelihood*prior at these current active points
        likeactive = zeros(1,numactive) + sum(betaprior(1:numPCs)) ...
            + Tprior + dprior;
        %coherent case: L is product over individual detector likelihoods
        for j=1:numactive
            for i=1:detno
                lp = like_gauss_fspace_td(wave_ft(:,i), sigdet(:,i), ...
                    deltaF, len, f, lowfreq_index, highfreq_index, activeT(j), ...
                    @findy, V_ft, activebeta(j,:), actived(j), numPCs);
                
                likeactive(j) = likeactive(j) + lp;% + sum(betaprior);
                
            end
            like_sum = logplus(like_sum,likeactive(j)-numactive);
        end
        
        like_avg(k) = like_sum;
    end
    
    %pick the bins with highest sumed likelihood
    like_avg_max = max(like_avg);
    like_avg_min = min(like_avg);
    idcs = find(like_avg_max-like_avg < 0.1*(like_avg_max-like_avg_min));
    %minT = mid_T0(idcs(1))-bin_size/2;
    %maxT = mid_T0(idcs(end))+bin_size/2;
    Tprior = -log(maxT-minT);
end
%time shift values for active points
activeT = rand(1,numactive)*(maxT-minT) + minT;

% calculate the likelihood*prior at these current active points
likeactive = zeros(1,numactive) + sum(betaprior(1:numPCs)) ...
    + Tprior + dprior;

%coherent case: L is product over individual detector likelihoods
for j=1:numactive
    for i=1:detno
        lp = like_gauss_fspace_td(wave_ft(:,i), sigdet(:,i), ...
            deltaF, len, f, lowfreq_index, highfreq_index, activeT(j), ...
            @findy, V_ft, activebeta(j,:), actived(j), numPCs);
        
        likeactive(j) = likeactive(j) + lp;% + sum(betaprior);
        
    end
end

Z(1) = -1e37;

j = 2;

logw(1) = log(1 - exp(-1/numactive));

% start nested sampling loop - use John Veitch's critereon of continuing until
% i <= numactive*infosafe*H, where we steal infosafe = 1.5 from John's code
infosafe = 1.5;
H(1) = 0;

%doPlot = 1;
if(doPlot)
    figure('OuterPosition',[70,100,1200,600]);
end

% Define the output directory so that we can save intermittently
resultsdir=run_name;
mkdir(resultsdir)

verbose=1;

while j-2 <= numactive * infosafe * H(j-1)

    if verbose
        disp(sprintf('j-2: %.2f | stop: %.2f', j-2, numactive * infosafe * H(j-1)))
    end

    if(doPlot)
        activebeta_avg = mean(activebeta);
        activeD_avg = mean(actived);
        activeT_avg = mean(activeT);
        
    
        y_fft = findy(V_ft, activebeta_avg, activeD_avg, numPCs);
        y = ifft([0; y_fft; 0; flipud(conj(y_fft))]);
        y_sig_fft = fft(wave0)*sample_deltaT;
    
      %plot waveforms corresponding average of active points, in time domain   
       % subplot(1,3,1),plot((-100:100),wave0(4000:4200),(-100:100),y(4000:4200))
    
      %same, in frequency domain, plus signal received by detector    
%         subplot(1,2,1),plot((50:4000),abs(wave_ft(50:4000,1)),...
%             (50:4000),abs(y_sig_fft(50:4000)),...
%             (50:4000),abs(y_fft(50:4000)));
         h = subplot(1,2,1);
            plot(wave0,'r');
%             hold on
            plot(y)
    
            min1 = min(activebeta(:,1));
            max1 = max(activebeta(:,1));
            min2 = min(activebeta(:,2));
            max2 = max(activebeta(:,2));
            min3 = min(activebeta(:,3));
            max3 = max(activebeta(:,3));
    
            % creates a plot of points for betas 1, 2 and 3, zooms in on final outcomes
            subplot(1,2,2),plot3(activebeta(:,1), activebeta(:,2), activebeta(:,3), '.');
            subplot(1,2,2),xlim([min1 max1]);
            subplot(1,2,2),ylim([min2 max2]);
            subplot(1,2,2),zlim([min3 max3]);
            subplot(1,2,2),text(max1,max2,max3,['Iterations: ' num2str(j-2)]);
            
        drawnow;
%         Movie(j) = getframe(h);
    end
    
    % find minimum of likelihoods
    [Lmin, idx] = min(likeactive);
    
    % get the log weight
    logWt = Lmin + logw(j-1);
    
    Z(j) = logplus(Z(j-1), logWt);
    
    % get values for using in posterior calculation
    Lw(j-1) = logWt;
    betas(j-1,:) = activebeta(idx,:);
    Ts(j-1) = activeT(idx);
    Distance(j-1) = actived(idx);
    
    % calculate the information H
    H(j) = exp(logWt - Z(j))*Lmin + exp(Z(j-1) - Z(j)) *...
        (H(j-1) + Z(j-1)) - Z(j);
    
    Activebeta=activebeta;

    if (doDistance && doTimeshift)
        Activebeta(:,numPCs+1)=actived;
        Activebeta(:,numPCs+2)=activeT;
    elseif (doDistance)
        Activebeta(:,numPCs+1)=actived;
    elseif (doTimeshift)
        Activebeta(:,numPCs+1)=activeT;
    end
    
    % update the log width
    logw(j) = logw(j-1) - 1/numactive;
    
    % will need to change 2nd argument to number of PC's used
    if mod(j-2, numPCs + doDistance + doTimeshift)== 0
        cholmat = cholcov(cov(Activebeta));
    end
    
    good = 0;
    
    % incase no useable point is returned in the first MCMC keep repeating
    % until one is
    while good == 0
        % from the active points randomly pick a starting point for the chain
        randval = ceil(rand(1)*numactive);
        
        betachain = zeros(nits+1, numPCs);
        Tchain = zeros(nits+1, 1);
        dchain = zeros(nits+1, 1);
        
        betachain(1,:) = activebeta(randval,:);
        Tchain(1) = activeT(randval);
        dchain(1) = actived(randval);
        
        Lcur = likeactive(randval);
        
        acc = 0; % count accepted points
        
        % start MCMC chain
        %if verbose; disp('sampling new point'); end
        for k=1:nits
            % generate new position from proposal
            gasvals = randn(numPCs + doTimeshift + doDistance,1);
            newvals = cholmat*gasvals;
            
            betanew = betachain(k,:) + newvals(1:numPCs)';
            
            if(doTimeshift)
                Tnew = Tchain(k) + newvals(end);
            else
                Tnew = Tchain(k);
            end
            if(doDistance)
                dnew = dchain(k) + newvals(end-doTimeshift);
            else
                dnew = dchain(k);
            end
            
            % if new point is not within prior range then reject point
            if any((betanew < minbeta(1:numPCs)) | (betanew > maxbeta(1:numPCs))) || ...
                    dnew < mind || dnew > maxd || ...
                    Tnew < minT || Tnew > maxT
                %disp('outside prior')
                
                betachain(k+1,:) = betachain(k,:);
                Tchain(k+1) = Tchain(k);
                dchain(k+1) = dchain(k);
                continue;
            end
            
            % get likelihood for new position
            %coherent case: L is product over individual detector likelihoods
            Lnew = sum(betaprior(1:numPCs)) + Tprior + dprior;
            for i=1:detno
                %disp('computing likelihood')
                lp = like_gauss_fspace_td(wave_ft(:,i), sigma, ...
                    deltaF, len, f, lowfreq_index, highfreq_index, Tnew, ...
                    @findy, V_ft, betanew, dnew, numPCs);
                
                Lnew = Lnew + lp;% + sum(betaprior) ;
            end
            
            % check whether to accept or reject the new point
            Lrat = Lnew - Lcur;
            
            if Lrat - log(rand) >= 0 && Lnew > Lmin
                
                betachain(k+1,:) = betanew;
                Tchain(k+1) = Tnew;
                dchain(k+1) = dnew;
                Lcur = Lnew;
                acc = acc + 1;
            else
                
                betachain(k+1,:) = betachain(k,:);
                Tchain(k+1) = Tchain(k);
                dchain(k+1) = dchain(k);
            end
        end % MCMC sampling
        
        if acc > 0
            good = 1;
        end
    end 
    % substitute the final values of the MCMC chain into the position of
    % Lmin
    activebeta(idx,:) = betachain(end,:);
    activeT(idx) = Tchain(end);
    actived(idx) = dchain(end);
    
    likeactive(idx) = Lcur;
    
    j=j+1;
    
    %if rem(j,500) == 0
    %    disp('saving workspace')
    %    workspace_filename = strcat('workspace_iteration',num2str(j),'.mat');
    %    save([resultsdir '/' workspace_filename])
    %end
end

% we're no longer reducing the prior distribution, so weights stay the same
% as (1/N)*X_j
logwend = logw(j-1) -logw(1) - log(numactive);

% add the remaining points
for k=1:numactive
    logw(j+k-1) = logwend;
    Z(j+k-1) = logplus(Z(j+k-2), likeactive(k) + logw(j+k-1));
    Lw(k+j-2) = likeactive(k) + logw(j+k-1);
    betas(k+j-2,:) = activebeta(k,:);
    Ts(k+j-2) = activeT(k);
    Distance(k+j-2) = actived(k);
end

%reduce size of variables to be exported
Z = Z(1:j-1+numactive);
Lw = Lw(1:j-2+numactive);
betas = betas(1:j-2+numactive,:);
Ts = Ts(1:j-2+numactive);
Z_end= Z(j-1+numactive);
Bayes=Z(j-1+numactive)-sum(evnoise);
distance= Distance(1:j-2+numactive);

% calculates the posterior distribution on the beta distributions
rnums = rand(1,length(Lw));

wt = logw(1:length(Lw)) + Lw;
maxwt = max(wt);
idx = find(wt' > maxwt+log(rnums));

postbetas = betas(idx,:);
postT= Ts(idx,:);
postdis= distance(idx,:);

posterior_params_savefile = ['smee_output_' catalogue num2str(wv) '_model' model '_PCs' num2str(numPCs)...
    '_detno' num2str(detno) '_' typeofscaling strrep(num2str(scaling), '.', 'p') '_seed' num2str(seed)];
save([resultsdir '/' posterior_params_savefile],'catalogue','wv','model','betas','activebeta','detno',...
    'Ts','seed','distance','Lw', 'evnoise','Z_end','Bayes','SNR','postbetas','postT','postdis','effective_distance', 'numPCs');

% save a simple txt file with a line of output
outFile = sprintf('%s/smee_results.txt', resultsdir);
fid = fopen(outFile, 'a');
fprintf(fid, '%i %1.4e %1.4e %1.4e %i \n', wv, evnoise, Bayes, SNR, numPCs);
fclose(fid);

%print evidence for model
fprintf(1, 'Network: log(%s Evidence) = %f\n', model, Z_end);
fprintf(1, 'log(nested sampling Bayes factor vs Noise) = %f\n', Bayes);

%print extra output data
fprintf(1, 'No. of Iterations: %i\n', j);
fprintf(1,'H = %f\n', H(j-1));
for cnt=1:numPCs
    fprintf(1,'<B_%i> = %f\n',cnt,mean(activebeta(:,cnt)));
end
fprintf(1,'<T> = %f\n',mean(activeT));
fprintf(1,'<D> = %f\n',mean(actived));

% Save betas from nested sampling
%BETAS=activebeta;
%betas_savefile = ['betas_model' model wv '-cat' catalogue '-seed' seed];
%save([resultsdir betas_savefile],'BETAS');

clear H Lcur betas Lmin Lnew Lrat acc activebeta betachain cholmat betanew V fs...
    gasvals good idx infosafe likeactive logWt logdbeta Lw logwend logw maxwt...
    mlog2 newvals randval rnums wt Z lp l j k Activebeta ActiveT activeT Ts Tchain Tnew model
toc
exit
