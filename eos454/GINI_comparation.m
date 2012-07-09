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
D = D(1:1024,30,1:128);
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% %% another data
% clear;
% D = ReadSuFast('data_ex3.su');
% D = D(1:128,1:1024);
% [nt, ns] = size(D);
% D = vec(D);

%% compare different measurements for different transforms 
% -- itself, wavelet, curvelet
W = opWavelet(nt,ns);
W_D = W*D;
C = opCurvelet(nt, ns);
C_D = C*D;
DCT = opDCT(nt*ns);
DCT_D = DCT*D;
FFT = opFFT2d(nt,ns);
FFT = opFunction(nt*ns,nt*ns,FFT);
FFT_D = FFT*D;
% 0 norm
D_0 = nnz(D);
WD_0 = nnz(W_D);
CD_0 = nnz(C_D);
DCTD_0 = nnz(DCT_D);
FFTD_0 = nnz(FFT_D);
% 1 norm
D_1 = norm(D,1);
WD_1 = norm(W_D,1);
CD_1 = norm(C_D,1);
DCTD_1 = norm(DCT_D,1);
FFTD_1 = norm(FFT_D,1);
% gini index
D_GI = gini(D);
WD_GI = gini(W_D);
CD_GI = gini(C_D);
DCTD_GI = gini(DCT_D);
FFTD_GI = gini(FFT_D);
fprintf('Transform: 0_norm     1_norm      Gini \n')
fprintf('Itself:    %3.2e, %3.2e, %1.6f\n',D_0, D_1,D_GI);
fprintf('Wavelet:   %3.2e, %3.2e, %1.6f\n',WD_0, WD_1,WD_GI);
fprintf('Curvelet:  %3.2e, %3.2e, %1.6f\n',CD_0, CD_1,CD_GI);
fprintf('DCT:       %3.2e, %3.2e, %1.6f\n',DCTD_0, DCTD_1,DCTD_GI);
fprintf('FFT:       %3.2e, %3.2e, %1.6f\n',FFTD_0, FFTD_1,FFTD_GI);

