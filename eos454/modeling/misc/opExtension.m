classdef opExtension < opSpot
% Extension operator. Pads input with constant values on either side
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   op = opExtension(n,nb,flag)
%
% input:
%   n    - length of input
%   nb   - number of points to add on both sides
%   flag - 0: padd with zeros, 1: padd with boundary value.
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       nb;
       a;
    end % Properties


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opExtension(n,nb,a)
          m = n + 2*nb;
          if nargin<3
              a = 1;
          end
          % Construct operator
          op = op@opSpot('Extension', m, n);
          op.nb = nb;
          op.sweepflag = 1;
          op.a = a;
       end % Constructor
       
       function out = test(op)
       
           x  = randn(op.n,1);
           y  = opExtension_intrnl(op.n,op.nb,op.a,x,1);
           xt = opExtension_intrnl(op.n,op.nb,op.a,y,-1);
           
           e1 = abs((y'*y)/(x'*xt)-1);
           
           if~(e1<1e-10); fprintf(2,'opExtension: adjoint test failed, error = %g\n',e1); end
           
           out = (e1<1e-10);
       end
       

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = opExtension_intrnl(op.n,op.nb,op.a,x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opExtension_intrnl(n,nb,a,x,mode)
    nx = size(x,2);
    if mode == 1
        if a
            y = [repmat(x(1,:),nb,1);x;repmat(x(end,:),nb,1)];
        else
            y = [zeros(nb,nx);x;zeros(nb,nx)];
        end
    else
       y = x(nb+1:end-nb,:);
       if a
            y(1,:)   = y(1,:) + sum(x(1:nb,:),1);
            y(end,:) = y(end,:) + sum(x(end-nb+1:end,:),1);
       end
    end
end

