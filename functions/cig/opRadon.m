classdef opRadon < opSpot
% SPOT operator for NMO correction. 
% 
% use:
%   op = opRadon(t,h,q,power);

% input:
%   t     - time vector in seconds 
%   h     - offset vecror in meters
%   q     - radon parameter
%   power - 1: linear radon, 2: parabolic radon
%
% output:
%   op    - spot operator of size [length(t)*length(q) length(t)*length(h)]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        t,h,q,power;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opRadon(t,h,q,power)
            op = op@opSpot('opRadon', length(t)*length(q), length(t)*length(h));
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = []; 
            op.sweepflag = false;
            op.t         = t;   
            op.h         = h;                      
            op.q         = q;
            op.power     = power;
        end %constructor
    end
    
    methods ( Access = protected )
        function out = multiply(op,in,mode)
            
            if mode==1
                in  = reshape(in,length(op.t),length(op.h));
                out = lpradon(in,op.t,op.h,op.q,op.power,1);
            else
                in  = reshape(in,length(op.t),length(op.q));
                out = lpradon(in,op.t,op.h,op.q,op.power,-1);
            end
            out = out(:);
        end %multiply
    end
end


