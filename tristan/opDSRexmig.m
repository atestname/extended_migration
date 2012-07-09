function op = opDSRexmig(dim,t,x,y,z,v)

% return extended image

% Determine sizes
m = prod(dim);

% Construct the operator
fh = @(data,mode) opDSRmig_intrnl(data, t, x, y, z, v ,mode);
op = opFunction(m, m, fh);


function y = opDSRmig_intrnl(data,t,x,y,z,v,mode)

if (mode == 0)
    
    y = {m,m,[0,0,0,0],{'opDSRmig'}};
    
elseif (mode == 1)
    
    if norm(x-y)
        error('source and receiver grids must be the same!');
    end
    
    if length(v) - length(z)
        error('velocity must be the same size as z')
    end
    
    % initialize image
    image = zeros(length(z),length(x));
    if nargout > 1
        eximage = zeros(length(z),length(x),length(x));
    end
    
    % depth step
    dz = z(2) - z(1);
    
    % compute first
    dt   = t(2)-t(1);
    nt   = length(t);
    dx   = x(2) - x(1);
    nx   = length(x);
    dy   = y(2) - y(1);
    ny   = length(y);
    
    fny  = .5/dt;
    if mod(nt,2) == 1
        % df   = 2*fny/(nt-1);
        df   = 2*fny/nt;
        f    = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/nt;
        f    = [0:df:fny -fny+df:df:-df];
    end
    fny  = .5/dx;
    if mod(nx,2) == 1
        %df   = 2*fny/(nx-1);
        df   = 2*fny/nx;
        kx   = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/nx;
        kx   = [0:df:fny -fny+df:df:-df];
    end
    fny  = .5/dy;
    if mod(ny,2) == 1
        %df   = 2*fny/(ny-1);
        df   = 2*fny/ny;
        ky   = [0:df:fny -fny:df:-df];
    else
        df   = 2*fny/ny;
        ky   = [0:df:fny -fny+df:df:-df];
    end
    
    % f-k-k grid
    [ff,kkx,kky] = ndgrid(f,kx,ky);
    
    % loop over depth levels
    for iz = 1:length(z)
        
        % survey sinking
        data = reshape(data,dim);
        E = DSR_step(data,ff,kkx,kky,iz*dz,v(iz));
        
        % imaging condition
        %image(iz,:) = diag(squeeze(E(1,:,:)));
        eximage(iz,:,:) = squeeze(E(1,:,:));

    end % end loop over depth levels
    y = eximage;
    
else % mode = -1
    
    %
end % end of if mode == 0



    function v = DSR_step(u,ff,kkx,kky,dz,v)
        
        % f-k-k transform
        spec = fft(u,[],1);
        spec = fft(spec,[],2);
        spec = fft(spec,[],3);
        
        % DSR operators
        Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
        Py = 2*pi*sqrt((ff/v).^2-kky.^2);
        
        Px = real(Px)+1i*abs(imag(Px));
        Py = real(Py)+1i*abs(imag(Py));
        
        % apply phase shift
        spec = exp(1i*abs(dz)*(Px + Py)).*spec;
        
        % inverse fkk transform
        v = ifft(spec,[],1);
        v = ifft(v,[],2);
        v = ifft(v,[],3);
        
        % take real part
        v = real(v);
    end

end

end
