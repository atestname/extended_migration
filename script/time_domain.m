%% 
% time domain.

%addpath of some functions and operators
addpath /Volumes/Users/linamiao/Documents/Tools/Matlabtools/tuning/
addpath(genpath('./simu_functions'))
addpath ../data
addpath ../functions

clear
% close all
%% read data
[data,SuTraceHeaders,SuHeader]=ReadSu('data_ex3.su');

% source and receiver coordinates per trace
for k=1:length(SuTraceHeaders)
    xs(k)= SuTraceHeaders(k).SourceX;
    xr(k)= SuTraceHeaders(k).GroupX;
end

% note that the sample interval is given in miliseconds!
t = [0:SuHeader.ns-1]'*SuHeader.dt*1e-6;

% source and receiver coordinate vectors
xs = unique(xs); 
xr = unique(xr); 

% reshape data into cube
data = reshape(data,length(t),length(xr),length(xs));

% window and subsample
data(1:100,:,:) = 0;
data = data(101:101+127,1:10:end,:);
data = data(:,1:end-1,1:end-1);%SR2MH cannot handle odd no.
xr = xr(1:10);nrec = length(xr);
xs = xr;nsrc = length(xs);
t = t(101:127+101);
nt = length(t);

%% simulate simutaneous experiment
% nt = length(t);
% ns = length(xs);
% nr = length(xr);
% 
% % measurement matrix
% p = .5;
% RM = opSimSourcePeriodTimeDither([nt nr ns], [floor(p*nt*ns),nr], ns);
% 
% % observed data from simutaneous experiment
% simD = RM*vec(data); 
% simd = reshape(simD,64,11,11);
% t = t(1:nt/2);

%% simutaneous experiment with randomized source superposition
% source superposition
 ncsrc = 5;
M   = opGaussian(nsrc,ncsrc);
RM  = opKron(opDirac(nt),opDirac(nrec),M');

simD = RM*vec(data);
simd = reshape(simD,nt,nrec,ncsrc);
%% sparsity transform
% wavelet operator along time axis
W = opSplineWavelet(nt, 1, nt, 3, 5);

% curvelet operator along source-receiver oordinates
C = opCurvelet(nr,ns,6,16,1,'ME',0);

% oppKron2Lo : kronecker tensor product to act on a distributed vector
S = opKron2Lo(C, W', 1); 

%% DSR extended migration
z = 10:10:1000;
v = 1000*ones(length(z),1);

dim = [nt,nrec,nsrc];
K = opDSR(dim,t,xr,xs,z,v);
M = K*vec(data); % DSR result if not sim
dim1 = [nt,nrec,ncsrc];
K1 = opDSR(dim1,t,xr,xs,z,v);
M1 = K1*simD; % DSR resulr if sim


%% inversion via joint sparsity
A = RM*K'*S';
B = reshape(simD,nt*ns,nr);
fid = fopen('log_log.txt', 'w'); 
opts = spgSetParms('optTol',1e-4);
sigma = 1e-3;
[X_hat,R,G,INFO] = spg_mmv(A,B,sigma,opts);


%% evaluate result
M1 = S'*X_hat;
SNR1 = snr(vec(M),vec(M1));

D1 = K'*M1;
SNR2 = snr(vec(D1),vec(data));

