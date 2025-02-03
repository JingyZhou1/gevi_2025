function tc = convCut(stim, irf)

% for temporal convolution

n_terms = max(length(stim), length(irf));

%n_terms = length(stim);

tc = conv(stim, irf, 'full');
tc = tc(1 : n_terms);

end