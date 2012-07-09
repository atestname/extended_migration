function op = opDSRmig(dim,t,x,y,z,v)

% Determine sizes
m = prod(dim);
n = length(z)*dim(2)*dim(3);

% Construct the operator
fh = @(data,mode) opDSRmig_intrnl(data, t, x, y, z, v ,mode);
op = opFunction(n, m, fh);


function y = opDSRmig_intrnl(data,t,x,y,z,v,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDSRmig'}};
        
    elseif (mode == 1)
        
        y = DSR_mig(data,t,x,y,z,v);
        
    else % mode = -1
        
        y = DSR_inv(data,t,x,y,z,v);
        
    end % end of if mode == 0




end

end