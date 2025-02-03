function error = fitBoxCart_exponential_filtered(stim, prm, dt, t)


% forward model
filter = gammaPDF(t, prm(5), 2);
filter = [0, filter(1 : length(filter) - 1)];
prd = convCut(stim, filter) * prm(1) + prm(2) * exp(-t * prm(3)) + prm(4);

% compute error
error = sum((prd - dt).^2);


%% visualize

% figure (100), clf
% plot(prd), hold on
% plot(dt)







end