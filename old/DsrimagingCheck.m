% this script is to test whether DSRImaing is working properly

% addpath of some functions and operators
addpath /Volumes/Users/linamiao/Documents/Tools/Matlabtools/tuning/

%% use a point reflector
%% generate data
z = 0:10:500;
x = 0:10:500;
[zz,xx] = ndgrid(z,x);

% background velocity [m/s]
v0 = 2500 + 0*xx;

% perturbation
epsilon = .1;
dv = 0*xx; dv((zz-250).^2 + (xx-250).^2 <= 100^2) = epsilon*2500;
v = v0+dv;
figure;imagesc(v);xlabel('x [m]');ylabel('z [m]');zlabel('velocity [m/s]');zlim([2500 3000]);

% Modeling
% grid, z = o(1) + [0:n(1)-1]*d(1), z = o(1) + [0:n(2)-1]*d(2);
model.o = [0 0];
model.d = [10 10];
model.n = [51 51];
model.nb = [10 10 0];
 
% frequencies [Hz]
model.freq = linspace(1,125,2^3); nfreq = length(model.freq);

% Ricker wavelet peak frequency and phase shift
model.f0 = 10;
model.t0 = 0;

% source and receiver positions
model.zsrc = 10;
model.xsrc = 0:10:500; nsrc = length(model.xsrc);
model.zrec = 10;
model.xrec = 0:10:500; nrec = length(model.xrec);

% define point sources, each column of this matrix represents a source
% function defined on the grid {model.zsrc,model.xsrc}. A point source is
% represented as a spike on one of the gridp-points. If we take Q to be an
% identity matrix, each column represents a point-source on a different
% gridpoint.

Q = speye(nsrc);

% define model in [km^2/s^2]
m = 1e6./(v0(:) + dv(:)).^2; 

% create data
D = F(m,Q,model,1);

% reshape vectorized data into data-cube for plotting purposes
D = reshape(D,[nrec,nsrc,nfreq]); 

% plot frequency slices
figure;imagesc(real(D(:,:,1)))



%% migration
data = permute(D,[3,1,2]);
image = DSRImagingDFT(data,v(1:length(z)),z,x,model.freq,1); 

figure; imagesc(real(image))





