function r2 = fitModel(t, prm, stim, basis_s, basis_f, dt)

prd =  normModel2(t, prm, stim, basis_s, basis_f);

%%
r2 = sum((prd - dt).^2);

%% visualize

% figure (100), clf
% plot(prd), hold on
% plot(dt)


end