% vsd_analysis

%% load vsd TF data

vsd = [];

[vsd.tf, vsd.mtf, vsd.lintf, vsd.normtf, vsd.bs_normtf] = gevi_loadData('vsd_tf', 140, 7);

%% visualize vsd data

figure (1), clf
plot(squeeze(vsd.normtf(:, :, 3))')

%%

n_boots = 50;

vsd.tofit = (vsd.bs_normtf(:, :, 1 : 7) + vsd.bs_normtf(:, :, 8 : end))/2;
vsd.tofit = reshape(vsd.tofit, n_boots, []);

%% visualize vsd data

figure (1), clf
subplot(311), 
plot(vsd.mtf'), box off
set(gca, 'xtick', 0 : 140 : 140 * 7)

subplot(312)
%tmp = (vsd.normtf(:, :, 1 : 7) + vsd.normtf(:, :, 8 : end)
plot(reshape(vsd.normtf(:, :, 1 : 7), 6, [])'), box off
set(gca, 'xtick', 0 : 140 : 140 * 7)

subplot(313),
plot(vsd.tofit'), box off
set(gca, 'xtick', 0 : 140 : 140 * 7)

%% make stimulus and basis functions

stim = mkStim('TF', 140);
stim = stim(:, 3 : 9);
stim = reshape(stim, [], 1);

t_lth = 140;

% % make basis functions
t = linspace(0.01, t_lth/100, t_lth);
%
nFast = 6;
nSlow = 5;

f_lth = 50;

fBasis  = mkBasis(t(1 : f_lth), nFast, 'fast_longdt');
sBasis  = mkBasis(t, nSlow, 'slow');

% make basis functions
basis = concatenateBasisAcrossConditions(fBasis, sBasis, stim, t);

%% linear model fit

vsd.lin_w = []; vsd.lin_prd = [];

fBasis_long = basis(end - nFast + 1 : end, :);

for k = 1 : n_boots
    vsd.lin_w(k, :)   = least_square(fBasis_long', vsd.tofit(k, :)');
    vsd.lin_prd(k, :) = vsd.lin_w(k, :) * fBasis_long;
    vsd.lin_r2(k)     = compute_r2(vsd.tofit(k, :), vsd.lin_prd(k, :));
end

figure (2), clf
subplot(211), 
shadedErrorBar(1 : length(vsd.tofit), mean(vsd.tofit), std(vsd.tofit)), hold on
plot(mean(vsd.lin_prd), 'b-'), hold on

%% normalization model fit

vsd.norm_prm = []; vsd.norm_prd = [];

for k = 1 : n_boots
    init = [vsd.lin_w(k, :), 0.03, 4, .1, .6, 0.007, 20, 1];
    vsd.norm_prm(k, :) = fminsearch(@(x) fitModel1(t, x, stim',  fBasis_long, vsd.tofit(k, :), 140 * 0 + 1 : 140 * 7), init);
    vsd.norm_prd(k, :) = real(normModel3(t, vsd.norm_prm(k, :), stim', fBasis_long, 140 * 0 + 1 : 140 * 7));
    % compute r2
    vsd.norm_r2(k) = compute_r2(vsd.tofit(k, :), vsd.norm_prd(k, :));
end

subplot(212), 
shadedErrorBar(1 : length(vsd.tofit), mean(vsd.tofit), std(vsd.tofit)), hold on
plot(mean(vsd.norm_prd), 'r-'), hold on

%% exporting parameters

fileLoc = '/Users/jyzhou/Desktop/GEVI/';

% save file

save(fullfile(fileLoc, 'vsd_parameters'), 'vsd');
