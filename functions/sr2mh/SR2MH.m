% SR2MH(D, MODE) reads a seismic line D in source-receiver (s-r)
% coordinates (or midpoint-offset (m-h) coordinates) and returns it in
% midpoint-offset coordinates (or source-receiver coordinates). MODE
% determines the type of conversion.
% MODE = 1 ---> s-r to m-h
% MODE = 2 ---> m-h to s-r

% Dependencies: opSR2MH.m
%
% Hassan Mansour, 2011
%
% Haneet Wason
% April, 2012 : Wrapped opSR2MH into a SPOT operator.
% May, 2012   : Combined conversion from s-r to m-h, and m-h to s-r.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Dout = SR2MH(D, MODE)

if MODE == 1
    
    % Data dimensions
    dim = size(D);
    
    % Operator that converts source-receiver slices to midpoint-offset slices
    SR = opSR2MH(dim(2));
    
%     % Convert the subsampled measurements from s-r to m-h
%     Dout = zeros(dim(1), dim(2), 2*dim(2));
%     
%     for t = 1:dim(1)
%         
%         % Do conversion of full data to midpoint offset
%         Dout(t,:,:) = reshape(SR*(real(D(t,:)))',1,dim(2),2*dim(2));
%         
%     end
    
    opMH = opKron(SR, opEye(dim(1)));
    temp = opMH*D(:);
    Dout = reshape(temp, dim(1), dim(2), 2*dim(2));
    
elseif MODE == 2
    
    % Data dimensions
    dim = size(D);
    
    % Operator that converts midpoint-offset slices to source-receiver slices
    SR = opSR2MH(dim(2));
   
%     % Conversion of full data from m-h to s-r
%     Dout = zeros(dim(1), dim(2), dim(2));
%     
%     for t = 1:dim(1)
%         
%         Dout(t,:,:) = reshape(SR'*(real(D(t,:)))',1,dim(2),dim(2));
%         
%     end
    
    opSR = opKron(SR, opEye(dim(1)));
    temp = opSR'*D(:);
    Dout = reshape(temp, dim(1), dim(2), dim(2));
    
    
end

