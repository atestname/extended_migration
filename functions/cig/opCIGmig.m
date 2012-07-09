function op = opCIGmig(f,ky,kh,h,z)
% Determine sizes
% nf = length(f);
nz = length(z);
nh = length(h);
ny = length(ky);

% compute p
[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
p = unique(p);
np = length(p);

% clear f p;

m = nz*nh*ny;
n = nz*np*ny;

% Construct the operator
fh = @(data,mode) opCIGmig_intrnl(data,f,ky,kh,h,z,mode);
op = opFunction(n, m, fh);


function y = opCIGmig_intrnl(data,f,ky,kh,h,z,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDSRmig'}};
        
    elseif (mode == 1)
        
        y = cig_imaging(data,f,ky,kh,h,z);
        y = vec(y);
        
    else % mode = -1
        
        y = cig_inv_imaging(data,f,ky,kh,h,z);
        y = vec(y);
        
    end % end of if mode == 0




end

end