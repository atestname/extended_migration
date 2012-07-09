function op = opDFT(dim)

% Construct the operator
fh = @(data,mode) opDFT_intrnl(data,mode);
op = opFunction(dim, dim, fh);


function y = opDFT_intrnl(data,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        y = fft(data);
        y = vec(y);
        
    else % mode = -1
        
        y = ifft(data);
        y= vec(y);
        
    end % end of if mode == 0




end

end