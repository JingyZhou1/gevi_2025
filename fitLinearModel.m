function mse = fitLinearModel(t, prm, stim, fBasis, dt, scale_range)

prd = linearModel(t, prm, stim, fBasis, scale_range);


mse = sum((prd - dt).^2);


end