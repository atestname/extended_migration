function image = DSR_mig_spot_2(data,t,x,y,z,v)

% x is xr -- receiver grid
% y is xs -- source grid
% t       -- time axis
% z       -- vertical axis
% v       -- velocity as a function of z, must be the same size of z

if norm(x-y)
    error('source and receiver grids must be the same!');
end

if length(v) - length(z)
    error('velocity must be the same size as z')
end

% initialize image
image = zeros(length(z),length(x));%,length(x));

% depth step 
dz = z(2) - z(1);

% basic computation of frequencies and wave numbers
[f kx ky] =  fkk(t,x,y);

% f-k-k grid
[ff,kkx,kky] = ndgrid(f,kx,ky);

% loop over depth levels
% initialize
% data  = opFFT1(size(ff))*vec(data);
for iz = 1:length(z)
    
    % survey sinking
    DSR = opDSR_step_spot(ff,kkx,kky,dz,v(iz));
    E = DSR*vec(data); % sink for dz, and E is in f-k-k domain
    data = opFFT1(size(ff))'*E; % updata time domain data for next layer
    S = opKron(opDirac(length(ky)),opDirac(length(kx)),ones(1,length(f)));
    SE = S*E;% sum over all the freqs 
    SE = reshape(SE,length(ky),length(kx));
    image(iz,:) = diag((SE));
    %image(iz,:,:) = ((SE));
end












