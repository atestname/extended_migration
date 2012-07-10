function [Imag] = cig_imaging(F,t,ky,kh,h,z)
% input 
%     F : wavefield vector(z,h,y)
%     t : time axis
%     kh: offset wavenumber
%     ky: midpoint wavenumber
%
% output
%     Imag : common p image gather(z_hat,p,m)
%     z_hat: vertical coordinate after radon transform
%     p    : wave parameter

% determine dimentions
nt = length(t);
ny = length(ky);
nh = length(kh);
nz = length(z);

% compute f
dt   = t(2)-t(1);
fny  = .5/dt;
if mod(nt,2) == 1
    df   = 2*fny/(nt-1);
    f    = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nt;
    f    = [0:df:fny -fny+df:df:-df];
end

% compute p
[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
p = unique(p);

% form Radon operator
R = opRadon(z,h,p,1);

% deal with the midpoint dimension
Im = opDirac(ny);

% the imaging operator
CIG = opKron(R,Im); 

% apply the operator
%FF = reshape(F,nh*nz,ny);
Imag = CIG*vec(F);
