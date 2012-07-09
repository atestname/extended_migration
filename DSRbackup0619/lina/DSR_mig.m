function image = DSR_mig(data,t,xr,xs,z,v)

if norm(xr-xs)
    error('source and receiver grids must be the same!');
end

% initialize image
image = zeros(length(z),length(xr));

% depth step 
dz = z(2) - z(1);

% loop over depth levels
for iz = 1:length(z)
    % survey sinking
    E = DSR_step(data,t,xr,xs,iz*dz,v);
    
    % imaging condition
    image(iz,:) = diag(squeeze(E(1,:,:)));
end
% end loop over depth levels


function v = DSR_step(u,t,x,y,dz,v)

tmax = t(end);
dt   = t(2)-t(1);
xmax = x(end);
dx   = x(2) - x(1);
ymax = y(end);
dy   = y(2) - y(1);
f    = [0:1/tmax:.5/dt -.5/dt:1/tmax:-1/tmax];
kx   = [0:1/xmax:.5/dx -.5/dx:1/xmax:-1/xmax];
ky   = [0:1/ymax:.5/dy -.5/dy:1/ymax:-1/ymax];

% f-k-k transform
spec = fft(u,[],1);
spec = fft(spec,[],2);
spec = fft(spec,[],3);


% f-k-k grid
[ff,kkx,kky] = ndgrid(f,kx,ky);

% DSR operators
Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

% apply phase shift
spec = exp(1i*abs(dz)*(Px + Py)).*spec;

% inverse fkk transform
v = ifft(spec,[],1);
v = ifft(v,[],2);
v = ifft(v,[],3);

% take real part
v = real(v);

