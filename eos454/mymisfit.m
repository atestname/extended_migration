function [f,g] = mymisfit(m,Q,D,model)

% Write a matlab function mymisfit(m,Q,D,model) that takes in a model, 
% a source function, observed data and a structure containing the modeling 
% parameters and returns the value of the misfit and the gradient for the 
% given model.

% compute predicted data
[D_pre,J] = F(m,Q,model);

% calculate misfit
f = .5 * norm(D_pre - D(:),2).^2;

% calculate gradient
g = J'*(D_pre - D(:));
