function op = opCIGmig(t,ky,kh,h,z)
% Determine sizes
nz = length(z);
nh = length(h);
ny = length(ky);

% compute f
dt   = t(2)-t(1);
fny  = .5/dt;
nt = length(t);
if mod(nt,2) == 1
    df   = 2*fny/(nt-1);
    f    = [0:df:fny -fny:df:-df];
else
    df   = 2*fny/nt;
    f    = [0:df:fny -fny+df:df:-df];
end

% compute p
[ff,kkh] = ndgrid(f,kh);
p = kkh./ff;
p = unique(p);
np = length(p);

clear f p;

m = nz*nh*ny;
n = nz*np*ny;

% Construct the operator
fh = @(data,mode) opCIGmig_intrnl(data,t,ky,kh,h,z,mode);
op = opFunction(n, m, fh);


function y = opCIGmig_intrnl(data,t,ky,kh,h,z,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDSRmig'}};
        
    elseif (mode == 1)
        
        y = cig_imaging(data,t,ky,kh,h,z);
        y = vec(y);
        
    else % mode = -1
        
        y = cig_inv_imaging(data,t,ky,kh,h,z);
        y = vec(y);
        
    end % end of if mode == 0




end

end