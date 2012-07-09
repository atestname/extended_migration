%% 
% this is something in freq domain, forget it for a second, let us move to
% time domain.



% addpath of some functions and operators
addpath /Volumes/Users/linamiao/Documents/Tools/Matlabtools/tuning/

clear
%close all
%% generate data 
z = 0:10:1000;
x = 0:10:1000;
[zz,xx] = ndgrid(z,x);

% background velocity [m/s]
v0 = 1000 + 0.*zz;
figure;imagesc(v0);title('background velcocity')

% true velocity model
dv = 0*xx;
% dv(zz >= 0) = 1600; 
% v(zz >= 300) = 2000;
dv(zz >= 220) = 600;
dv(zz >= 250) = 0;
figure; imagesc(dv);

% perturbation 
%dv = v - v0;
figure;imagesc(dv);title('velocity perterbation')

% Modeling
model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [30 30 0];
  
% frequencies [Hz]
model.freq = 5:20:125; 
model.freq = 5:5:20; 
nfreq = length(model.freq);

% Ricker wavelet peak frequency and phase shift
model.f0 = 20;
model.t0 = -.08;

% source and receiver positions
model.zsrc = 10;
model.xsrc = 0:10:1000; nsrc = length(model.xsrc);
model.zrec = 10;
model.xrec = 0:10:1000; nrec = length(model.xrec);

% define point sources, each column of this matrix represents a source
% function defined on the grid {model.zsrc,model.xsrc}. A point source is
% represented as a spike on one of the gridp-points. If we take Q to be an
% identity matrix, each column represents a point-source on a different
% gridpoint.
Q = speye(nsrc);

% define model in [km^2/s^2]
m = 1e6./(v0(:) + dv(:)).^2; 

% create data
[DD,J] = F(m,Q,model,1);

% linearize data
D = J*m;

% reshape vectorized data into data-cube for plotting purposes
D = reshape(D,[nrec,nsrc,nfreq]); 

% plot frequency slices
figure;imagesc(real(D(:,:,1)));title('frequency slice of observed data')

% re-organize data into frequency-receiver-source order
D = permute(D,[3,1,2]);

% save modelled data
save model_data DD D model x z dv v0 ;



  
%% simulate simutaneous experiment
load model_data
nf = length(model.freq);
addpath(genpath('./simu_functions'))
% marine acquisition
p = .5;
RM = opSimSourceRandTimeDither([nf,nr,ns],[p*nf*ns,1],ns);
% simutaneous data
simD = RM*vec(D);

% sparsity representation
C = opCurvelet(nf, ns);

% suppose we have the expended migration operator opDSRmig
K = opDSRmig([nf,nr,ns],model.freq,x,z,v);
% operator
A = RM*K*C';

%% recover original data via joint sparsity
opts = spgSetParms('optTol',1e-4);
sigma = 1e-3;
[X_hat,R,G,INFO] = spg_mmv(A,B,sigma,opts);
f = C'*X_hat;
SNR = snr(vec(D),f);
