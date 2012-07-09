function [D,J] = F(m,Q,model)
% Frequency domain FD modeling operator
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use: 
%   [D,J] = F(m,Q,model)
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition. size(Q,2) determines the number of
%                       sources
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%   model.pml         - [n, beta, l] (optional)
%   model.order       - order of stencil {2,4} (optional).
%
% output:
%   D  - Data cube (nrec x nsrc x nfreq) as vector. nsrc  = size(Q,2);
%                                                   nrec  = length(zrec)*length(xrec) 
%                                                   nfreq = length(freq)
%   J  - Jacobian as SPOT operator

% optional arguments
nb    = [1 1]*floor(max(model.n)/2);
beta  = 100;
l     = 2;
order = 2;
if isfield(model,'pml')
    nb    = [1 1]*model.pml(1);
    beta  = model.pml(2);
    l     = model.pml(3);
end
if isfield(model,'order')
    order = model.order;
end

% physical grid
N = prod(model.n);

% comp. grid
ot = model.o-nb.*model.d;
dt = model.d;
nt = model.n+2*nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),nb(2)),opExtension(model.n(1),nb(1)));

% model parameter: squared slowness [s^2/m^2] on computational grid.
nu = 1e-6*Px*m;

% loop over frequencies
D = zeros(nsrc*nrec,nfreq);
for k = 1:nfreq
    Hk = Helm2D(model.freq(k),nu,ot,dt,nt,nb,beta,l,order);
    Uk = Hk\(w(k)*Ps'*Q);
    D(:,k) = vec(Pr*Uk);
end

% output
D = vec(D);
J = opDF(m,Q,model);
