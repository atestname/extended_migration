classdef opDF < opSpot
% SPOT wrapper for DF.m
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   J = opDF(m,Q,model)
%
% see DF.m for further documentation
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        %put any variables that you need to store here
        % for instance, any variables you pass to the function handle
        % should be stored here instead, and then they will be available in
        % the multiply function
        % m, n, linear, cflag, children - are already provided by
        % superclass no need to redefine them here.
        mt,Q,model,nfreq,nt;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opDF(mt,Q,model)
           %put all preprocessing thats done in the sparco operator before
           %passing back the function handle.
            
            nsrc  = size(Q,2);
            nrec  = length(model.xrec)*length(model.zrec);
            nfreq = length(model.freq);
            m = nsrc*nrec*nfreq;
            n = length(mt);
           % preprocessing ends here
           
           % now we instantiate the operator and set the property values
           op = op@opSpot('opDF', m, n);
           op.cflag     = 1;  % is your operator complex?
           op.linear    = 1; % is your operator linear?
           op.children  = []; % optional - these are other operators your 
           op.sweepflag = 1;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;
           op.nfreq     = nfreq;
           op.nt        = nsrc*nrec;
       end %constructor
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode == 1
                y = DF(op.mt,op.Q,x, 1,op.model);
           else %adjoint
                y = DF(op.mt,op.Q,x,-1,op.model);  
           end
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
