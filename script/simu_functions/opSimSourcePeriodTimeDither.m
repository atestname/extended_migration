function op = opSimSourcePeriodTimeDither(dim,dims,nss)

% opSimSourceRandTimeDither creates peroid time-shifts between sequentially deployed sources.
% Input parameters:  
%        dim  : Dimensions of the input data
%        dims : Dimensions of simultaneous data
%        nss  : Number of sources


% Determine sizes
m = prod(dims);
n = prod(dim);


% Time index for each source 
shifts = linspace(1,dims(1)-dim(1),nss);
shifts = floor(shifts);

% Shot positions
shotIND = linspace(1,nss,nss);


% Construct the operator
fh = @(x,mode) opSimSourcePeriodTimeDither_intrnl(m, n, x, mode, dim, dims, shifts, shotIND, nss);
op = opFunction(m, n, fh);


function y = opSimSourcePeriodTimeDither_intrnl(m, n, x, mode, dim, dims, shifts, shotIND, nss)

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


