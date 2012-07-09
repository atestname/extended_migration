% try it in only two dimension, frequency and receiver

% % addpath of some functions and operators
% addpath /Volumes/Users/linamiao/Documents/Tools/Matlabtools/tuning/
% 
% clear
% close all
% %% generate data 
% z = 0:10:630;
% x = 0:10:630;
% [zz,xx] = ndgrid(z,x);
% 
% % background velocity [m/s]
% v0 = 1100 + 0*xx;
% 
% % perturbation 
% epsilon = .1;
% dv = 0*xx; 
% dv((zz-10)>= 200) = epsilon*1100;
% dv((zz-100)>= 200) = epsilon*2000;
% dv((zz-200)>= 300) = epsilon*4000;
% 
% v = v0+dv;
% figure;imagesc(v);xlabel('x [m]');ylabel('z [m]');zlabel('velocity [m/s]');zlim([2100 3000]);
% 
% % Modeling
% % grid, z = o(1) + [0:n(1)-1]*d(1), z = o(1) + [0:n(2)-1]*d(2);
% model.o = [0 0];
% model.d = [10 10];
% model.n = [64 64];
% model.nb = [3 3 0];
%   
% % frequencies [Hz]
% model.freq = linspace(-250,250,64); nfreq = length(model.freq);
% 
% % Ricker wavelet peak frequency and phase shift
% model.f0 = 10;
% model.t0 = 0;
% 
% % source and receiver positions
% model.zsrc = 10;
% model.xsrc = 0:10:630; nsrc = length(model.xsrc);
% model.zrec = 10;
% model.xrec = 0:10:630; nrec = length(model.xrec);
% 
% % define point sources, each column of this matrix represents a source
% % function defined on the grid {model.zsrc,model.xsrc}. A point source is
% % represented as a spike on one of the gridp-points. If we take Q to be an
% % identity matrix, each column represents a point-source on a different
% % gridpoint.
% 
% Q = speye(nsrc);
% 
% % define model in [km^2/s^2]
% m = 1e6./(v0(:) + dv(:)).^2; 
% 
% % create data
% D = F(m,Q,model);
% 
% % reshape vectorized data into data-cube for plotting purposes
% D = reshape(D,[nrec,nsrc,nfreq]); 
% 
% D = permute(D,[3,1,2]);
% 
% save DD D model x z v ;
% 
% % plot frequency slices
% figure;imagesc(real(D(:,:,1)))
% 
% Fx = opDirac(64);
% Ft = opDFT(64);
% F = opKron(Fx,Fx,Ft);
% dd = F'*vec(D);
% d = reshape(dd,64,64,64);
% figure;imagesc(real(d(:,:,1)));colormap(gray);
% 
%   
% 
% 
% 
% % keyboard;
%% migration
clear;close all;
load DD;
figure;imagesc(v);colorbar;
data = D;
% image = DSRImagingDFT2D(data,v(1:length(z)),z,x,model.freq,1); 
%[image,uu] = DSRImagingDFT3D(data,v(1:length(z)),z,x,model.freq,1); 
image = my_mig(data,model.freq,x,z,v(1:length(z)),1);
figure; imagesc(real(image));colorbar;











