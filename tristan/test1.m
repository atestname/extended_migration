function output = test1(data,t,x,y,z,v)

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


% depth step 
dz = z(2) - z(1);

% basic computation of frequencies and wave numbers
[f kx ky] =  fkk(t,x,y);

% f-k-k grid
[ff,kkx,kky] = ndgrid(f,kx,ky);

% initial data
%image = zeros(size(ff));
image = zeros(length(z),length(x));

% operator
S = opKron(opDirac(length(ky)),opDirac(length(kx)),ones(1,length(f)));
F = opFFT1(size(ff))';

% copy data
% R = repmat(opDirac(numel(data)),length(z),1);
% data = R*vec(data);

% loop over depth levels
for iz = 1:length(z)
    D = opDSR_step_spot(ff,kkx,kky,dz,v(iz));
    if iz == 1
        A = D;
    else
        A = D*F*A;
    end
    op = S*A;
    output = op*vec(data);
    output = reshape(output,length(x),length(y));
    %image(iz,:,:) = ((output(:,:)));
    image(iz,:) = diag((output(:,:)));
end

output = image;




% for iz = 1:length(z)
%     image(iz,:) = diag(squeeze(output(iz,:,:)));
% end
    











