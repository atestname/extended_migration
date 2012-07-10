function wave = DSR_inv(eximage,t,x,y,z,v)

nz = length(z);
nx = length(x); 
ny = length(y);
nt = length(t);

if size(size(eximage)) ~= 3
    eximage = reshape(eximage,nz,nx,ny);
end

data = zeros(nz,nt,nx,ny);

% basic computation of frequencies and wave numbers
dt   = t(2)-t(1);
nt   = length(t);
dx   = x(2) - x(1);
nx   = length(x);
dy   = y(2) - y(1);
ny   = length(y);
fny  = .5/dt;
if mod(nt,2) == 1
    df   = 2*fny/(nt-1);
    f    = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nt;
    f    = [0:df:fny -fny+df:df:-df];
end
fny  = .5/dx;
if mod(nx,2) == 1
    df   = 2*fny/(nx-1);
    kx   = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nx;
    kx   = [0:df:fny -fny+df:df:-df];
end
fny  = .5/dy;
if mod(ny,2) == 1
    df   = 2*fny/(ny-1);
    ky   = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/ny;
    ky   = [0:df:fny -fny+df:df:-df];
end

% fkk grid
[ff kkx kky] = ndgrid(f,kx,ky);
nf = length(f);

% depth interval
dz = z(2)-z(1);

% loop through all depth levels
for iz = 1:nz
    SE = eximage(iz,:,:);
    % scattering zero time to all frequencies component
    S = opKron(opDirac(ny),opDirac(nx),ones(1,nf));
    E = S'*vec(SE);
    data(iz,:,:,:) = DSR_inv_step(E,ff,kkx,kky,iz*dz,v(iz));
end

wave = squeeze(mean(data,1));
%wave = squeeze(data(2,:,:,:));

function u = DSR_inv_step(v,ff,kkx,kky,dz,vel)
% v in frequency space wavenumber domain
% fkk transform
if size(size(v)) ~= 3
    v = reshape(v,size(ff));
end
spec = fft(v,[],2);
spec = fft(spec,[],3);

% DSR operators
Px = 2*pi*sqrt((ff/vel).^2-kkx.^2);
Py = 2*pi*sqrt((ff/vel).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

% inverse apply phase shift
E = exp(1i*abs(dz)*(Px + Py));
spec = reshape(spec,size(E));
spec = spec./E;
% spec = diag(vec(E)'*vec(spec));

% inverse f-k-k transform
u = ifft(spec,[],1);
u = ifft(u,[],2);
u = ifft(u,[],3);

