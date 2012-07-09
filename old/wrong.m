% something wrong
clear;
close all;

% time coordinate in seconds as column vector
t = [0:.004:1.020]';

% spatial coordinate in meters as row vector
x = 0:10:1270;
z = 0:10:1270;

% generate 2D grid:
[tt,xx] = ndgrid(t,x);

% velocity
% background velocity [m/s]
v = 3000;

figure;imagesc(v); colormap(gray);
% Source wavefield is impulsive source at t=0.1 and x=500,
source = (tt-.1).*exp(-1e-3*(xx-500).^2 - 1e3*(tt-.1).^2);
figure;imagesc(real(source));colormap(gray)

%%  if source is given in frequency domain
Ft = opDFT(size(u,1));
Fr = opDirac(size(u,2));
F = opKron(Fr,Ft);
%%
u = source;

Ft = opDFT(size(u,1));
Fr = opDFT(size(u,2));
F = opKron(Fr,Ft);
U = F*vec(u);
fnyq = 1. / (2*(t(2)-t(1)));
nf = size(u,1);
df = 2*fnyq/nf;
f = [0:df:fnyq-df -fnyq:df:-df]; 
%f = 2*pi*f;

dx = (x(2)-x(1));
fnyq = 1. / (2*(dx));
nf = size(u,2);
dff = 2*fnyq/nf;
fr = [0:dff:fnyq-dff -fnyq:dff:-dff];
k = fr;
U = reshape(U,size(u,1),size(u,2));
figure;imagesc(real(U));colormap(gray)

%% how about this
% [U,f] = fftrl(u,t);
% U = fft(U,[],2);



%% if fktran everything is fine
  [U1,f1,k1] = fktran(u,t,x);
%  source_fk1 = U1;
% 
%% well, then I will center 0
% U = fftshift(U,2);
% k = sort(k);
% f = sort(f);
source_fk = U;
figure;imagesc(real(source_fk));colormap(gray)
% and remove negative frequency part
% UU = U(1:length(t)/2+1,:);
% source_fk = UU;
% f = f(length(t)/2:end);
% figure;imagesc(real(source_fk));colormap(gray)

%% 
[ff,kk] = ndgrid(f,k);
dz = 10;
kz = 2*pi*sqrt((ff/v).^2 - kk.^2);
% fix sign of kz
kz = -real(kz)+1i*abs(imag(kz));
W  = exp(500*1i*kz);



%% use ifktran
% source_extrap = ifktran(W.*source_fk,f,k);

%% fft
% ext = W.*source_fk;
% ext = ifftshift(ext);
% ee = ifftrl(ext,f);
% ee = ifft(ee,[],2);
% ee = reshape(ee,length(t)+1,length(x));
% source_extrap = ee;

%% not using ifktran
ext = W.*source_fk;
% ext = ifftshift(ext);
ee = F'*vec(ext);
ee = reshape(ee,length(t),length(x));
source_extrap = ee;

%%  trim back to orinigal size
source_extrap = source_extrap(1:length(t),1:length(x));


% plot
figure;imagesc(x,t,real(source_extrap));colormap(gray);
xlabel('x [m]');ylabel('t [s]');title('wavefield at z=500 m.');




%% kz = 2*pi*sqrt((ff/v).^2 - kk.^2);

% wrong sign on imaginary part:
W1 = exp(500*1i*(-real(kz) - 1i*abs(imag(kz))));

% wrong sing on real part:
W2 = exp(500*1i*( real(kz) + 1i*abs(imag(kz))));

% plot
% extrapolate and transform back to t-x
source_extrap1 = ifktran(W1.*source_fk,f,k);
source_extrap2 = ifktran(W2.*source_fk,f,k);

% trim back to orinigal size
source_extrap1 = source_extrap1(1:length(t),1:length(x));
source_extrap2 = source_extrap2(1:length(t),1:length(x));

% plot
figure;imagesc(x,t,source_extrap1);colormap(gray);
xlabel('x [m]');ylabel('t [s]');title('wavefield at z=500 m, wrong sign imag. part');
figure;imagesc(x,t,source_extrap2);colormap(gray);
xlabel('x [m]');ylabel('t [s]');title('wavefield at z=500 m, wrong sign real part.');

% transform W to the t-x domain
Wtx = ifktran(W,f,k);

% it turns out that |Wtx| is `flipped' in the |x| direction and of the
% wrong size. We fix it here:
Wtx = fftshift(Wtx,2);
Wtx = Wtx(1:251,14:end-14);

% plot
figure;imagesc(x,t,Wtx,[-1 1]*1e1);colormap(gray)
xlabel('x [m]');ylabel('t [s]');title('extrapolator in the t-x domain');








