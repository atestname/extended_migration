function op = opSimSourceRandTimeDither(dim,dims,nss)

% opSimSourceRandTimeDither creates random time-shifts between sequentially deployed sources.
% Input parameters:  
%        dim  : Dimensions of the input data
%        dims : Dimensions of randomized/simultaneous data
%        nss  : Number of sources


% Determine sizes
m = prod(dims);
n = prod(dim);


% Possible random time indices
ind = randperm(dims(1) - dim(1));


% Time index for each source 
shifts = sort(ind(1:nss));


% Shot positions
shotIND = linspace(1,nss,nss);


% Construct the operator
fh = @(x,mode) opSimSourceRandTimeDither_intrnl(m, n, x, mode, dim, dims, shifts, shotIND, nss);
op = opFunction(m, n, fh);


function y = opSimSourceRandTimeDither_intrnl(m, n, x, mode, dim, dims, shifts, shotIND, nss)

if (mode == 0)
    
    y = {m,n,[0,0,0,0],{'opSimSourceRandTimeDither'}};
    
elseif (mode == 1)
    
    x = reshape(x,dim);
    y = zeros(dims);
    
    for j = 1 : nss;
        shot = x(:,:,shotIND(j));
        y(shifts(j):shifts(j) + dim(1) - 1,:) = y(shifts(j):shifts(j) + dim(1) - 1,:) + shot;
    end
    
    y = y(:);
    
else
    
    x = reshape(x,dims);
    y = zeros(dim);
    
    for j = 1 : nss;
        y(:,:,shotIND(j)) = x(shifts(j):shifts(j) + dim(1) - 1,:);
    end
    
    y = y(:);
    
end


