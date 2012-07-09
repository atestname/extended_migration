function op = opDSR_step_spot(ff,kkx,kky,dz,v)

m = numel(ff);

% Construct the operator
fh = @(u,mode) opDSR_step_spot_intrnl(u,ff,kkx,kky,dz,v,mode);
op = opFunction(m, m, fh);


function y = opDSR_step_spot_intrnl(u,ff,kkx,kky,dz,v,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        y = DSR_step(u,ff,kkx,kky,dz,v);
        y = vec(y);
        
    else % mode = -1
        
        y = DSR_inv_step(u,ff,kkx,kky,dz,v);
        y = vec(y);
        
    end % end of if mode == 0




end

end