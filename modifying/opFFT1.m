function op = opFFT1(dim)

% fft and inverse fft on the first dimension of input

m = prod(dim);

% Construct the operator
fh = @(data,mode) opFFT3_intrnl(data,dim,mode);
op = opFunction(m, m, fh);


function y = opFFT3_intrnl(u,dim,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        if size(size(u))~=3
            u = reshape(u,dim);
        end
        spec = fft(u,[],1)/sqrt(size(u,1));        
        y = vec(spec);
        
    else % mode = -1
        
        if size(size(u))~=3
            u = reshape(u,dim);
        end
        v = ifft(u,[],1)*sqrt(size(u,1));
        y = vec(v);
        
    end % end of if mode == 0




end

end