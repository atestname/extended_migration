function op = opDSRmig(dim,tf,x,z,v)

% Determine sizes
m = prod(dim);

% Construct the operator
fh = @(data,mode) opDSRmig_intrnl(data, tf, x, z, v);
op = opFunction(m, m, fh);


function y = opDSRmig_intrnl(data, tf, x, z, v, mode)

if (mode == 0)
    
    y = {m,m,[0,0,0,0],{'opDSRmig'}};
    
elseif (mode == 1)
    
    % initialize image
    image = zeros(length(z),length(x));
    
    % depth step
    dz = z(2) - z(1);
    
    % transform operator
    Ft = opDFT(size(data,1));
    Ir = opDirac(size(data,2));
    Is = opDirac(size(data,3));
    F = opKron(Is,Ir,Ft);
    
    ppi = data; % initialize
    
    for iz = 1:length(z)
        ppi = my_step(ppi,tf,x,v(iz,:),dz,flag);
        if flag == 1
            % transform back to  t-s-r domain
            ppishift = ifftshift(ppi,1);
            ppitemp = F'*vec(ppishift);
        else
            ppitemp = ppi;
        end
        
        % update image
        ppitemp = reshape(ppitemp,length(tf),length(x),length(x));
        %clear ppitemp
        for ix = 1:length(x)
            image(iz,ix) = ppitemp(1,ix,ix);
        end
        
    end
    
else
    
    %
end


