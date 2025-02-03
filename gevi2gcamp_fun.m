function [prd, trans_weights] = gevi2gcamp_fun(gevi, prm, transbasis, gcamp, rng)

% assuming the first 3 parameters are for exponentiation

%mk_n = @(x, prm) prm(3) ./(1 + exp(prm(2) +  x * prm(1))) + prm(4);


%gevi_exp = exp(gevi .* prm(1) ) .* prm(2)  ;
gevi_exp = max(0, gevi) .^prm(2) * prm(1);
%gevi_exp = max(0, gevi).^ 1.7282* prm(1) ;
%gevi_exp = exp(max(0, gevi .* prm(2))) * prm(1);

%gevi_exp = mk_n(gevi, prm);

% make filters
pre = [];
for k = 1 : size(transbasis, 1)
    pre(k, :) = convCut(gevi_exp, transbasis(k, :));
end

% final prediction
trans_weights = least_square(pre', gcamp');

prd = trans_weights' * pre;

prd(rng) = prd(rng) * prm(3);


end