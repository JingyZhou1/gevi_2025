function prd = normModel3(t, prm, stim, basis_f, ctr_range)

% normalization with a scaling contrast factor
% different from normModel2, here the fast basis is not convolved with the
% stimulus yet



% prm(1 : 2) : filter for the first lope
% prm(3 : 4) : second lope of the filter

%% useful funcitons

norm_sum = @(x) x ./sum(x(:)); 
%norm_max = @(x) x ./max(abs(x(:)));


tf_range = 1 : ctr_range(1) - 1;
 
%% make linear predictions

nFast = size(basis_f, 1);

% for k = 1 : size(basis_f, 1)
%     f_prd(k, :) = convCut(basis_f(k, :), stim);
% end

num_rsp  = prm(1 : nFast) * basis_f;

%num_rsp(ctr_range) = num_rsp(ctr_range) * prm(end);

%% make denominator response

%prm(nFast + 1) = prm(nFast + 1);
%prm(nFast + 2) = prm(nFast + 2);
%prm(nFast + 5) = prm(nFast + 5);
%prm(nFast + 6) = prm(nFast + 6);

d_filter = [0, gammaPDF(t, prm(nFast + 1), prm(nFast + 2)) - prm(nFast + 4) * gammaPDF(t, prm(nFast + 5), prm(nFast + 6))];
d_filter = norm_sum(d_filter(1 : end - 1));

dem_rsp  = convCut(stim, d_filter);

%dem_rsp(ctr_range) = dem_rsp(ctr_range);

%% make final response

prd = num_rsp./max(1e-9, (prm(nFast + 3) + dem_rsp));

%prd(ctr_range) = prd(ctr_range) * prm(end);

prd(tf_range) = prd(tf_range) * prm(nFast + 7);

end