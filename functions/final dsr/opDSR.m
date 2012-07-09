function op = opDSR(dim,t,x,y,z,v)

% Determine sizes
m = prod(dim);

[f, ky, kh, ~] = kykh(t,x,y);
[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
p = unique(p);
n = length(z)*length(p)*length(ky);

% Construct the operator
fh = @(data,mode) opDSRmig_intrnl(data, t, x, y, z, v ,mode);
op = opFunction(n, m, fh);


function y = opDSRmig_intrnl(data,t,x,y,z,v,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDSR'}};
        
    elseif (mode == 1)
        
        y = DSR(data,t,x,y,z,v);
        y = vec(y);
        
    else % mode = -1
        
        y = inv_DSR(data,t,x,y,z,v);
        y= vec(y);
        
    end % end of if mode == 0




end

end