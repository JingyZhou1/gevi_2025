function prd = linearModel(t, prm, stim, fBasis, scale_range)


%% make linear model prediction without scaling

nFast = size(fBasis, 1);

prd = prm(1 : nFast) * fBasis;

%% scale the linear model prediction

prd(scale_range) = prd(scale_range) * prm(end);



end