% this is the back up for DSR_mig in frequency domain at Jun 20, 
% able to get migrated image, 
% without odd even check
% without velocity changed in horizontal axis
% without vector initial velocity

function image = DSR_mig_freq(data,f,xr,xs,z,v)

% DSR migration if data is given in frequency(space) domain

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
    E = DSR_step(data,f,xr,xs,iz*dz,v);
    
    % imaging condition
    image(iz,:) = diag(squeeze(E(1,:,:)));
end
% end loop over depth levels


function v = DSR_step(u,f,x,y,dz,v)

xmax = x(end);
dx   = x(2) - x(1);
ymax = y(end);
dy   = y(2) - y(1);
kx   = [0:1/xmax:.5/dx -.5/dx:1/xmax:-1/xmax];
ky   = [0:1/ymax:.5/dy -.5/dy:1/ymax:-1/ymax];

% f-k-k transform
spec = u;
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

% pad zero to other frequencies part
% frequency axis
n = length(f);
faxis = 0:5:125;
nf = length(faxis);
temp = zeros(nf,length(x),length(x));
for i = 1:n
    idx = find(faxis==f(i));
    temp(idx,:,:) = spec(i,:,:);
end

spec = temp;

% compensate negative frequency part 
if mod(n,2) == 0
   tmp = [spec ; conj(spec(end:-1:2,:,:))];
else
   tmp = [spec ; conj(spec(end-1:-1:2,:,:))];
end

% inverse fkk transform
v = ifft(tmp,[],1)*sqrt(n);
v = ifft(v,[],2);
v = ifft(v,[],3);

% take real part
v = real(v);

