function error = fitBoxCart(stim, prm, dt, t)


% forward model
prd = stim * prm(1) + prm(2) * t;

% compute error
error = sum((prd - dt).^2);


%% visualize

% figure (100), clf
% plot(prd), hold on
% plot(dt)







end