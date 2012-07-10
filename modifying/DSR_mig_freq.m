function image = DSR_mig_freq(data,f,x,y,z,v)

% DSR migration if data is given in frequency(space) domain

% x is xr -- receiver grid
% y is xs -- source grid
% f       -- time axis
% z       -- vertical axis
% v       -- velocity as a function of z, must be the same size of z


if norm(x-y)
    error('source and receiver grids must be the same!');
end

if length(v) - length(z)
    error('velocity must be the same size as z')
end

% initialize image
image = zeros(length(z),length(x));

% depth step 
dz = z(2) - z(1);

% compute first
dx   = x(2) - x(1);
nx   = length(x);
dy   = y(2) - y(1);
ny   = length(y);


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

% f-k-k grid
[ff,kkx,kky] = ndgrid(f,kx,ky);

% loop over depth levels
for iz = 1:length(z)
    % survey sinking
    E = DSR_step(data,ff,kkx,kky,iz*dz,v(iz));
    
    % imaging condition
    image(iz,:) = diag(squeeze(E(1,:,:)));
end
% end loop over depth levels


function v = DSR_step(u,ff,kkx,kky,dz,v)

% f-k-k transform
spec = u;
spec = fft(spec,[],2);
spec = fft(spec,[],3);

% DSR operators
Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

% apply phase shift
spec = exp(1i*abs(dz)*(Px + Py)).*spec;

% pad zero to other frequencies part
% frequency axis
% n = length(f);
 n = size(ff,1);
% faxis = 0:5:125;
% nf = length(faxis);
% temp = zeros(nf,length(x),length(x));
% for i = 1:n
%     idx = find(faxis==f(i));
%     temp(idx,:,:) = spec(i,:,:);
% end

% spec = temp;

% compensate negative frequency part 
if mod(n,2) == 1
   tmp = [spec ; conj(spec(end:-1:2,:,:))];
else
   tmp = [spec ; conj(spec(end-1:-1:2,:,:))];
end

% inverse fkk transform
v = ifft(tmp,[],1);%*sqrt(n);
v = ifft(v,[],2);
v = ifft(v,[],3);

% take real part
v = real(v);

