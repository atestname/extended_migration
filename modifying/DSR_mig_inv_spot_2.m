function data = DSR_mig_inv_spot_2(image,t,x,y,z,v)

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
    input = image(iz,:,:);
    % form all the operators
    if iz == 1
        DSR = opDSR_step_spot(ff,kkx,kky,dz,v(iz));
    else 
        DSR = DSR*F*
    end
    F = opFFT1(size(ff));
    S = opKron(opDirac(length(ky)),opDirac(length(kx)),ones(1,length(f)));
    A = S*D;
    output = output + A'*vec(input);
end












