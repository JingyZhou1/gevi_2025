function prd = filterModel(t, prm)

%% make the first gamma function

bump1 = norm_sum(gammaPDF(t, prm(1), prm(2)));
%bump1 = exp(-t * prm(1)) * prm(2);
%bump1 = t.^prm(1) * prm(2); % exp(-t * prm(1)) * prm(2);
gamma2 = norm_sum(gammaPDF(t, prm(3), prm(4)));

prd = bump1 * prm(5) + prm(6) * gamma2;
prd = norm_max([0, prd(1 : end - 1)]);

end