function prd = normModel2(t, prm, stim, basis_s, basis_f)

% prm(1 : 2) : filter for the first lope
% prm(3 : 4) : second lope of the filter

%% useful funcitons

norm_sum = @(x) x ./sum(x(:)); 
%norm_max = @(x) x ./max(abs(x(:)));

 
%% make linear predictions

nFast = size(basis_f, 1);

num_rsp  = prm(1 : nFast) * basis_f;
%num_rsp  = convCut(stim, num_rsp);

%% make denominator response

%prm(nFast + 1) = prm(nFast + 1);
%prm(nFast + 2) = prm(nFast + 2);
%prm(nFast + 5) = prm(nFast + 5);
%prm(nFast + 6) = prm(nFast + 6);

d_filter = gammaPDF(t, prm(nFast + 1), prm(nFast + 2)) - prm(nFast + 4) * gammaPDF(t, prm(nFast + 5), prm(nFast + 6));
d_filter = [0, d_filter(1 : end - 1)];

dem_rsp  = convCut(stim, d_filter);

%% make final response

prd = num_rsp./max(1e-9, (prm(nFast + 3) + dem_rsp));


end