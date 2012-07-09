% Sequential-source data reconstruction from randomized 'marine' acquisiton

% this is for the lab(project2) from eosc 454

% lina miao
% 74721119

% team member
% Ganason, Thuvendran
% Hess, Conrad

% The exercise deals with the reconstruction of a fully sampled data-volume
% from data that was subsampled by firing randomly dithered sources in marine.

%% Installation
% download and install 
clear; close all;
%% Data
% Number of time samples
nt = 1024;
% Number of sources
ns = 178;
% Number of receivers
nr = 178;

% Time sampling interval
dt = 0.004;

% Read data
D = ReadSuFast('GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);

% Select small subset
D = D(1:256,30,1:100);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)');
xlabel('Shot number'); ylabel('Time sample')
%% Set the patameters for randomized experiment
I  = eye(10);
RM1 = opSimSourceRandTimeDither([10 1 10], [5*10 1], 10);
RM2 = opSimSourcePeriodTimeDither([10 1 10], [5*10 1], 10);

% plot very long time series
figure;
plot(I(:),1:length(I(:)),'*');xlim([0.5 1.5]);ylabel('time samples');

% plot compressed series
figure;
plot(RM1*I(:),1:length(RM1*I(:)),'o');xlim([0.5 1.5]);ylabel('time samples');
figure;
plot(RM2*I(:),1:length(RM2*I(:)),'o');xlim([0.5 1.5]);ylabel('time samples');


%%
% Construct the sampling operator RM for p = 0.5 that works on the vectorized 
% version of the data using opSimSourceRandTimeDither.
p = .5;
D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
D_RM1 = opSimSourcePeriodTimeDither([nt,nr,ns],[p*nt*ns,1],ns);

% % Test the sampling operator with the dottest.
% x_test = rand(size(D_RM,2),1);
% y_test = rand(size(D_RM,1),1);
% left = y_test'*(D_RM*x_test);
% right = (D_RM'*y_test)'*x_test;
% error = norm(left-right);
% fprintf('In dottest error:%5.5e\n',error);

%%
% Generate simultaneous data simD and display the result.
simD1 = D_RM1*D;
simD2 = D_RM1*D;
figure;
imagesc(reshape(simD1,p*nt,ns)); colormap(gray); colorbar;
figure;
imagesc(reshape(simD2,p*nt,ns)); colormap(gray); colorbar;


%% sparsifying transform
% Use this to create a Curvelet SPOT operator:
C = opCurvelet(nt, ns);

% Transform the data into the Curvelet domain and plot the sorted coefficients 
C_D = C*D;
sort_CD = sort(abs(C_D),'descend');
figure;plot(sort_CD);title('sorted curvelet coefficients')
%% exercises
% Construct the measurement operator A. HINT: See 'Constructing a suitable 
% matrix' in Lab 7.
% Using spgl1, estimate the curvelet coefficients xest.

fid = fopen('log.txt', 'w'); 
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
p_list = [.3  .5 .7];
for i = 1:3
    p = p_list(i);
    D_RM1 = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
    simD1 = D_RM1*D;
    A = D_RM1*C';
    xest(:,i) = spgl1(A,simD1,0,1e-3,[],options);
    f(:,i) = C'*xest(:,i);
    snr(i) = snr(D,f(:,i)); 
end


figure; 
subplot(1,2,1);imagesc(reshape(f(:,1),nt,ns)); colormap(gray);
title(strcat(['p = .3, SNR=' num2str(snr(1)) 'dB']))
subplot(1,2,2);imagesc(reshape(f(:,1)-D,nt,ns)); colormap(gray);
title('difference')
    
figure; 
subplot(1,2,1);imagesc(reshape(f(:,2),nt,ns)); colormap(gray);
title(strcat(['p = .5, SNR=' num2str(snr(2)) 'dB']))
subplot(1,2,2);imagesc(reshape(f(:,2)-D,nt,ns)); colormap(gray);
title('difference')


figure; 
subplot(1,2,1);imagesc(reshape(f(:,3),nt,ns)); colormap(gray);
title(strcat(['p = .7, SNR=' num2str(snr(3)) 'dB']))
subplot(1,2,2);imagesc(reshape(f(:,3)-D,nt,ns)); colormap(gray);
title('difference')


%% This experiment I just can not make it work, hope you can continue 
% %% what can we do if we sample within the same time, say .25 of the original
% % time, what is the difference between periodic and dithering
% 
% I = eye(nt*ns/10,nt*ns);
% I_D = I*D;
% 
% % invert this peridoc sample using sparse optimization
% I_D_inv_spgl1 = spgl1(I*C',I_D,0,1e-3,[],options);
% figure; 
% subplot(1,2,1);imagesc(reshape(C'*I_D_inv_spgl1,nt,ns)); colormap(gray);
% title(strcat(['half time perodic spgl1, SNR=' num2str(snr(D,C'*I_D_inv_spgl1)) 'dB']))
% subplot(1,2,2);imagesc(reshape(C'*I_D_inv_spgl1-D,nt,ns)); colormap(gray);
% title('difference')
% 
% % dithering
% D_RM = opSimSourceRandTimeDither([nt,nr,ns],[.1*nt*ns,nr],ns);
% simD = D_RM*D;
% A = D_RM*C';
% xest = spgl1(A,simD,0,1e-3,[],options);
% figure; 
% subplot(1,2,1);imagesc(reshape(C'*xest,nt,ns)); colormap(gray);
% title(strcat(['dithering with half time, SNR=' num2str(snr(D,C'*xest)) 'dB']))
% subplot(1,2,2);imagesc(reshape(C'*xest-D,nt,ns)); colormap(gray);
% title('difference')

%% the bigger picture
clear;

% Data
nt = 1024;
ns = 178;
nr = 178;
dt = 0.004;
D = ReadSuFast('GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);
D = D(1:64,1:80,1:80);
[nt nr ns] = size(D);
D = D(:);

% measurement matrix
p = .5;
RM1 = opSimSourceRandTimeDither([nt nr ns], [p*nt*ns,nr], ns);
RM2 = opSimSourcePeriodTimeDither([nt nr ns], [p*nt*ns,nr], ns);

simD1 = RM1*D;
simD2 = RM2*D;

% wavelet operator along time axis
W = opSplineWavelet1(nt, 1, nt, 3, 5);

% curvelet operator along source-receiver oordinates
C = opCurvelet(nr,ns,6,16,1,'ME',0);

% oppKron2Lo : kronecker tensor product to act on a distributed vector
S = oppKron2Lo(C, W', 1);

% check the operator
adj = RM2'*simD1;
adj = reshape(adj,nt,nr,ns);
figure;imagesc(squeeze(adj(:,10,:))); colormap(gray)
figure;imagesc(squeeze(adj(:,:,10))); colormap(gray)
 

% l1 matrix
A = RM1*S';

% 
fid = fopen('log_log.txt', 'w'); 
options = spgSetParms('optTol', 1e-4, 'iterations', 10);%, 'fid', fid); 

xest = spgl1(A,simD1,0,1e-3,[],options);
f = S'*xest;
snr = snr(D,f); 




% figure; 
% subplot(1,2,1);imagesc(reshape(f,nt,nr,ns)); colormap(gray);
% title(strcat(['p = .5, SNR=' num2str(SNR) 'dB']))
% subplot(1,2,2);imagesc(reshape(f-D,nt,ns)); colormap(gray);
% title('difference')










%% 
