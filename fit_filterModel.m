function mse = fit_filterModel(t, prm, filter, f_lth)

prd = filterModel(t, prm);

mse = sum((prd(1 : f_lth) - filter).^2);

end