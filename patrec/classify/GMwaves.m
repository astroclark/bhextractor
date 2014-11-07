% GMMwaves.m

%% Load NR Data
data = load('allWavestheta-90.mat');

NRsignals = real(data.complex_strain);
% signals = signals + 0.01*randn(size(signals));


for i=1:length(signorm)
    NRsignals(:,i) = NRsignals(:,i)/sqrt(dot(NRsignals(:,i),NRsignals(:,i)));
end

%% Create a population of 'glitch' waveofrms

% --- SineGaussians
time = -0.5:1/length(NRsignals):0.5 - 1/length(NRsignals);

sineGausses = zeros(size(NRsignals));

f0s = 100+25*randn(1,50);
bws = 0.1+0.01*randn(1,50);
ps = 2*pi*randn(1,50);

for i=1:size(sineGausses,2)
    sineGausses(:,i) = real(gauspuls(time, f0s(i), bws(i)));
    sineGausses(:,i) = sineGausses(:,i)/dot(sineGausses(:,i),sineGausses(:,i));
end

% --- Gaussian monopulses
gmono = zeros(size(NRsignals));

f0s = 100+25*randn(1,50);

for i=1:size(gmono,2)
    gmono(:,i) = gmonopuls(time, f0s(i));
    gmono(:,i) = gmono(:,i)/dot(gmono(:,i),gmono(:,i));
end

%% combine

signals = [NRsignals];%, sineGausses, gmono];
% signals = NRsignals;
% signals = [sineGausses, gmono];

noise   = 0.0001*randn(size(signals));
% signals = signals + noise;

 


%% Perform PCA
disp('doing pca')
[coeff, score, latent, tsquared, explained] = ...
    pca(signals);


%% GMM fits with BIC

cum_explained = cumsum(explained);
k90 = find(abs(90-cumsum(explained)) == min(abs(90-cumsum(explained))));

reduced_score = score;%score(:,1:k90);

% Estimate number of clusters with AIC

nmax=5;

BIC = zeros(1,nmax);
GMModels = cell(1,nmax);

options.MaxIter = 1000;
return
disp('fitting GMMs')
for k = 1:nmax
    k
    GMModels{k} = fitgmdist(reduced_score,k, 'options', options);%,'Regularize', 0.1);
    BIC(k)= GMModels{k}.BIC;
end

[minBIC,numComponents] = min(BIC);
numComponents

BestModel = GMModels{numComponents}




% %% wavelet pca
% 
% x_orig = signals;
% x = signals + 0.01*randn(size(signals));
% 
% kp = 0;
% figure();
% for i = 1:4
%     subplot(4,2,kp+1), plot(x_orig(:,i)); axis tight;
%     title(['Original signal ',num2str(i)])
%     subplot(4,2,kp+2), plot(x(:,i)); axis tight;
%     title(['Noisy signal ',num2str(i)])
%     kp = kp + 2;
% end
% 
% level = 5;
% wname = 'sym4';
% npc = 'heur';
% [x_sim, qual, npc] = wmspca(x,level,wname,npc);
% 
% qual
% kp = 0;
% figure();
% for i = 1:4
%     subplot(4,2,kp+1), plot(x(:,i)); axis tight;
%     title(['Noisy signal ',num2str(i)])
%     subplot(4,2,kp+2), plot(x_sim(:,i)); axis tight;
%     title(['First PCA ',num2str(i)])
%     kp = kp + 2;
% end
% 
% 
% npc(1:3) = zeros(1,3);
% [x_sim, qual, npc] = wmspca(x,level,wname,npc);
% 
% figure();
% kp = 0;
% for i = 1:4
%     subplot(4,2,kp+1), plot(x(:,i)); axis tight;
%     title(['Noisy signal ',num2str(i)])
%     subplot(4,2,kp+2), plot(x_sim(:,i)); axis tight;
%     title(['Second PCA ',num2str(i)])
%     kp = kp + 2;
% end

