function idx = jitter1d(n, ns, gap)
% create roughly ns jitter-sampled indices from 1:n.
%
% idx = jitter1d(n, ns, gap)


% average spacing
d = floor(n/ns);

% jittered sampling
jittered_samples = zeros(1, n);

for i = 1:d:n
    
    
    perturb = sign(randn)*round(.5*rand*d);
    
    jittered_samples(min(max(1,i + perturb),n)) = 1;
    
end

% Collect indices of jittered locations
idx = find(jittered_samples);

end
