function op = opTEST(dim,t,x,y,z,v)

m = prod(dim);
n = length(z)*length(x)*length(y);

% Construct the operator
fh = @(data,mode) opTEST_intrnl(data,t,x,y,z,v,mode);
op = opFunction(n, m, fh);


function y = opTEST_intrnl(data,t,x,y,z,v,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        y = test1(data,t,x,y,z,v);
        
    else % mode = -1
        
        y = test2(data,t,x,y,z,v);
    end % end of if mode == 0

    y = vec(y);


end

end
