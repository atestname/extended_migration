classdef opNMO < opSpot
% SPOT operator for NMO correction. Uses nmo.m
% 
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca

% use:
%   op = opNMO(t,h,v);
%
% input:
%   t    - time vector [s]
%   h    - offset vector [m]
%   v    - NMO velocity [m/s] as vector of size [length(t) x 1].
%
% output
%   out  - spot operator of size [length(t)*length(h) x length(t)*length(h)]

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        t,h,v;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opNMO(t,h,v)
            N  = length(t)*length(h);
            op = op@opSpot('opNMO', N, N);
            op.cflag     = 1;  
            op.linear    = 1; 
            op.children  = []; 
            op.sweepflag = false;
            op.t         = t;   
            op.h         = h;                      
            op.v         = v;
        end %constructor
    end
    
    methods ( Access = protected )
        function out = multiply(op,in,mode)
            in = reshape(in,length(op.t),length(op.h));
            if mode==1
                out = nmo(in,op.t,op.h,op.v,1);
            else
                out = nmo(in,op.t,op.h,op.v,-1);
            end
            out = out(:);
        end %multiply
    end
end


