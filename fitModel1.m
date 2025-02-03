function r2 = fitModel1(t, prm, stim, basis_f, dt, ctr_range)

prd =  normModel3(t, prm, stim, basis_f, ctr_range);

r2 = sum((prd - dt).^2);

%% visualize

% figure (100), clf
% plot(prd, 'r-'), hold on
% plot(dt, 'k-')


end