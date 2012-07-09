function out = nmo(in,t,h,v,flag)
% NMO correction and adjoint 
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   out = nmo(in,t,h,v,flag)
%
% input:
%   in   - data matrix of size [length(t) x length(h)], each column is a trace 
%   t    - time vector [s]
%   h    - offset vector [m]
%   v    - NMO velocity [m/s] as vector of size [length(t) x 1].
%   flag - 1:forward, -1:adjoint
%
% output
%   out  - data matrix of size [length(t) x length(h)], each column is a trace 


[nt,nh] = size(in);
t = reshape(t,nt,1);
v = reshape(v,nt,1);

out = zeros(nt,nh);

for i = 1:nh
    tau = sqrt(t.^2 + h(i).^2./v.^2);
    A   = getLA(t,tau);
    if flag > 0
        out(:,i) = A*in(:,i);
    else
        out(:,i) = A'*in(:,i);
    end
end


end


function A = getLA(x1,x2)
    % find interior points
    ik = find(x2 >= x1(1) & x2 <= x1(end));
    % sizes
    n1=length(x1);n2=length(x2);nk=length(ik);
    % check
    if ~nk
        A = [];
        return;
    end
    
    % initialize stuff
    I = zeros(4*nk); J = I; S = I;
    a=1;b=2;c=3;d=4;
    l = 1;
    
    % loop
    for i = 1:nk
        k = ik(i);
        if x2(k)<x1(b)
            while (x2(k)<x1(b))&&(b-1>1)
                b=b-1;
            end
            a=b-1;c=b+1;d=c+1;
        elseif x2(k)>x1(c)
            while (x2(k)>x1(c))&&(c+1<n1)
                c=c+1;
            end
            a=c-2;b=c-1;d=c+1;
        end
        
        I(l)   = k; J(l)   = a; S(l)   = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(a)-x1(b))*(x1(a)-x1(c))*(x1(a)-x1(d)));
        I(l+1) = k; J(l+1) = b; S(l+1) = ((x2(k)-x1(a))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(b)-x1(a))*(x1(b)-x1(c))*(x1(b)-x1(d)));
        I(l+2) = k; J(l+2) = c; S(l+2) = ((x2(k)-x1(b))*(x2(k)-x1(a))*(x2(k)-x1(d)))/((x1(c)-x1(b))*(x1(c)-x1(a))*(x1(c)-x1(d)));
        I(l+3) = k; J(l+3) = d; S(l+3) = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(a)))/((x1(d)-x1(b))*(x1(d)-x1(c))*(x1(d)-x1(a)));
        l = l + 4;
    end
    % construct sparse matrix
    A = sparse(I(1:l-1),J(1:l-1),S(1:l-1),n2,n1); 
end
