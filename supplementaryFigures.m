% supplementary figures 

%% load data
fileLoc = '/Users/jyzhou/Desktop/GEVI/';
dataLoc = '/Users/jyzhou/Desktop/GEVI/code/';

% figure location
figureLoc = '/Users/jyzhou/Desktop/GEVI/figures/';

a = load(fullfile(dataLoc, 'preprocessed_data.mat'));

gcamp = a.gcamp;

n_boots = 50;

%% GCaMP preprocessing

% TEMPORAL FREQUENCY ------------------------------------------------------
% plot temporal frequency raw data
gcamp.tf_toplot = reshape(gcamp.tf(:, :, 3 : 9), 11, []);

t_plot_tf = 0.01 : 0.01 : 0.01 * 140 * 7;

c = parula(49);

figure (1), clf
% plot raw data
subplot(611),  hold on, colormap(c),
plot(t_plot_tf, gcamp.tf_toplot'), box off, 
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'), 

% plot detrended data
gcamp.tf_toplot_norm = reshape(gcamp.normtf(:, :, 3 : 9), 11, []);

subplot(612), 
plot(t_plot_tf, gcamp.tf_toplot_norm'), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:')

% plot bootstrapped data
gcamp.tf_toplot_bs = reshape(gcamp.bs_normtf(:, :, 3 : 9), n_boots, []);

subplot(613), 
plot(t_plot_tf, gcamp.tf_toplot_bs'), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:')

% CONTRAST ----------------------------------------------------------------
% plot contrast raw data
gcamp.ctr_toplot = reshape(gcamp.ctr(:, :, 3 : 8), 8, []);

t_plot_ctr = 0.01 : 0.01 : 0.01 * 140 * 6;

figure (1), 
% plot raw data
subplot(614), 
plot(t_plot_ctr, gcamp.ctr_toplot'), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 6, 'xticklabel', 1.4 * ones(1, 6))
plot([0, 0.01 * 140 * 6], [0, 0], 'k:')

% plot detrended data
gcamp.ctr_toplot_norm = reshape(gcamp.normctr(:, :, 1 : 6), 8, []);

subplot(615), 
plot(t_plot_ctr, gcamp.ctr_toplot_norm'), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 6, 'xticklabel', 1.4 * ones(1, 6))
plot([0, 0.01 * 140 * 6], [0, 0], 'k:')

% plot bootstrapped data
gcamp.ctr_toplot_bs = reshape(gcamp.bs_normctr(:, :, 1 : 6), n_boots, []);

subplot(616), 
plot(t_plot_ctr, gcamp.ctr_toplot_bs'), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 6, 'xticklabel', 1.4 * ones(1, 6))
plot([0, 0.01 * 140 * 6], [0, 0], 'k:')

%% VSD preprocessing

% Load VSD data
b = load(fullfile(fileLoc, 'vsd_parameters.mat'));
vsd = b.vsd;

% plot temporal frequency raw data
vsd.tf_toplot = reshape(vsd.tf(:, :, 3 : 9), 6, []);

figure (2), clf
subplot(311)
plot(t_plot_tf, vsd.tf_toplot), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'), 

% plot detrended temporal frequency data
vsd.tf_toplot_norm = reshape(vsd.normtf(:, :, 3 : 9), 6, []);

subplot(312), 
plot(t_plot_tf, vsd.tf_toplot_norm), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'),

% plot bootstrapped temporal frequency data
vsd.tf_toplot_bs = reshape(vsd.bs_normtf(:, :, 3 : 9), n_boots, []);

subplot(313), 
plot(t_plot_tf, vsd.tf_toplot_bs), box off, hold on
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'),

%% plot vsd model fits

vsd_toplot = (vsd.bs_normtf(:, :, 1 : 7) + vsd.bs_normtf(:, :, 8 : 14))/2;
vsd_toplot = reshape(vsd_toplot, 50, []);

figure (3), clf
subplot(211), 
shadedErrorBar(t_plot_tf, mean(vsd_toplot), std(vsd_toplot)), box off, hold on
plot(t_plot_tf, mean(vsd.lin_prd), 'b-')
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'),

subplot(212), 
shadedErrorBar(t_plot_tf, mean(vsd_toplot), std(vsd_toplot)), box off, hold on
plot(t_plot_tf, mean(vsd.norm_prd), 'r-')
set(gca, 'xtick', 0 : 1.4 : 1.4 * 7, 'xticklabel', 1.4 * ones(1, 7))
plot([0, 0.01 * 140 * 7], [0, 0], 'k:'),