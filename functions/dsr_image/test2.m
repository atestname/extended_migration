function data = test2(input,t,x,y,z,v)

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

if size(size(input)) ~= 3
    input = reshape(input,length(z),length(x),length(y));
end

% operator
S = opKron(opDirac(length(ky)),opDirac(length(kx)),ones(1,length(f)));
F = opFFT1(size(ff))';


% copy data operator
% R = repmat(opDirac(numel(input)),length(z),1);
% initialize output 
data = zeros(length(f),length(x),length(y));
data = vec(data);

% loop over depth levels
for iz = 1:length(z)
    D = opDSR_step_spot(ff,kkx,kky,dz,v(iz));
    if iz == 1
        A = D;
    else
        A = D*F*A;
    end
    
    op = S*A;
    
    data = op'*vec(input(iz,:,:)) + data;
end




% for iz = 1:length(z)
%     image(iz,:) = diag(squeeze(output(iz,:,:)));
% end
    











