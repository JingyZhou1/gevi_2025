function error = fitBoxCart_exponential(stim, prm, dt, t)


% forward model
prd = stim * prm(1) + prm(2) * exp(-t * prm(3)) + prm(4);

% compute error
error = sum((prd - dt).^2);


%% visualize

% figure (100), clf
% plot(prd), hold on
% plot(dt)







end