% Main analysis

% In this file, we generate Figure x - Figure x in the main text.

%% locations

% data location
dataLoc = '/Users/jyzhou/Desktop/GEVI/GEVI_paper/';

% figure location   
figureLoc = '/Users/jyzhou/Desktop/GEVI/figures/';

%% load both GEVI and GCaMP data

a = load(fullfile(dataLoc, 'preprocessed_data.mat'));

gevi = a.gevi;
gcamp = a.gcamp;

% pre-defined parameters
n_boots = 50;

%% concatenate GEVI and GCaMP data

% make gevi fitting data
gevi.tofit = cat(3, (gevi.bs_normtf(:, :, 1 : 7) + gevi.bs_normtf(:, :, 8 : 14))/2, ...
    (gevi.bs_normctr(:, :, 1 : 6) + gevi.bs_normctr(:, :, 7 : 12))/2); % averaging for the temporal frequency, and for the contrast data respectively


gevi.tofit = reshape(gevi.tofit, n_boots, []); % 50 bootstraps

% make gcamp fitting data
gcamp.tofit = cat(3, (gcamp.bs_normtf(:, :, 1 : 7) + gcamp.bs_normtf(:, :, 8 : 14))/2, ...
    (gcamp.bs_normctr(:, :, 1 : 6) + gcamp.bs_normctr(:, :, 7 : 12))/2);

gcamp.tofit = reshape(gcamp.tofit, n_boots, []);

figure (1), clf
subplot(211), plot(gevi.tofit'), title('GEVI bootstrapped data')
subplot(212), plot(gcamp.tofit'), title('GCaMP bootstrapped data')

%% Preparations for model fitting

t_lth = 140; % number of time points sampled

% make stimulus -------------------------------------------
stim1 = mkStim('TF', t_lth);
stim1 = stim1(:, 3 : 9); % the first two are blank stimulus conditions

stim2 = mkStim('contrast', t_lth);
stim2 = stim2(:, 3 : 8);

stim = reshape(cat(2, stim1, stim2), [], 1);

% make basis functions ------------------------------------
t = linspace(0.01, t_lth/100, t_lth);

nFast = 6;
nSlow = 5;

f_lth = 50; % length of the fast filter

fBasis  = mkBasis(t(1 : f_lth), nFast, 'fast_longdt');
sBasis  = mkBasis(t, nSlow, 'slow');

% make basis functions
basis = concatenateBasisAcrossConditions(fBasis, sBasis, stim, t);

%% GEVI AND GCAMP MODEL FIT

% Here, we only fit the fast component for both linear and normalization
% model

%% Fit linear model

gevi.lin_w = []; gevi.lin_prm = []; gevi.lin_prd = []; gcamp.lin_w = []; gcamp.lin_prd = []; gcamp.lin_prm = [];

fBasis_long = basis(end - nFast + 1 : end, :);

for k = 1 : n_boots
    % GEVI ----------------------------------------------------------
    % gevi linear weight
    gevi.lin_w(k, :) = least_square(fBasis_long', gevi.tofit(k, :)');

    init = [gevi.lin_w(k, :), 1];

    % fit scaled linear model
    gevi.lin_prm(k, :) = fminsearch(@(x) fitLinearModel(t, x, stim', fBasis_long, gevi.tofit(k, :), 140 * 7 + 1 : 140 * 13), init);
    


    % gevi linear predictions
    gevi.lin_prd(k, :) = linearModel(t, gevi.lin_prm(k, :), stim', fBasis_long, 140 * 7 + 1 : 140 * 13);
    % compute r2
    gevi.lin_r2(k) = compute_r2(gevi.tofit(k, :), gevi.lin_prd(k, :));

    % GCaMP ---------------------------------------------------------
    % gevi linear weight
    gcamp.lin_w(k, :) = least_square(fBasis_long', gcamp.tofit(k, :)');
    init = [gcamp.lin_w(k, :), 1];

    gcamp.lin_prm(k, :) = fminsearch(@(x) fitLinearModel(t, x, stim', fBasis_long, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13), init);

    % gevi linear predictions
    gcamp.lin_prd(k, :) = linearModel(t, gcamp.lin_prm(k, :), stim', fBasis_long, 140 * 7 + 1 : 140 * 13);

    % compute r2
    gcamp.lin_r2(k) = compute_r2(gcamp.tofit(k, :), gcamp.lin_prd(k, :));
end

% visualize linear model predictions
t_plot = [1 : length(gevi.tofit)]./100;

% visualize GEVI linear fit
figure (2), clf
subplot(2, 1, 1)
shadedErrorBar(t_plot, mean(gevi.tofit), std(gevi.tofit)), hold on
plot(t_plot, mean(gevi.lin_prd), 'b-'), box off
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

% visualize GCaMP linear fit
figure (3), clf
subplot(2, 1, 1)
shadedErrorBar(t_plot, mean(gcamp.tofit), std(gcamp.tofit)), hold on
plot(t_plot, mean(gcamp.lin_prd), 'b-'), box off
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

%% fit normalization model

gevi.norm_prm = []; gevi.norm_prd = []; gcamp.norm_prm = []; gcamp.norm_prd = [];

for k = 1 : n_boots
    % fit normalization model to GEVI -------------------------
    init = [gevi.lin_w(k, :), 0.03, 4, .1, .6, 0.007, 20, 1];
    gevi.norm_prm(k, :) = fminsearch(@(x) fitModel1(t, x, stim',  fBasis_long, gevi.tofit(k, :), 140 * 7 + 1 : 140 * 13), init);
    gevi.norm_prd(k, :) = real(normModel3(t, gevi.norm_prm(k, :), stim', fBasis_long, 140 * 7 + 1 : 140 * 13));
    % compute r2
    gevi.norm_r2(k) = compute_r2(gevi.tofit(k, :), gevi.norm_prd(k, :));

    % fit normalization model to GCaMP -------------------------
    init1 = [gcamp.lin_w(k, :), 0.03, 4, .1, .6, 0.05, 5, 1];
    gcamp.norm_prm(k, :) = fminsearch(@(x) fitModel1(t, x, stim',  fBasis_long, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13), init);
    gcamp.norm_prd(k, :) = real(normModel3(t, gcamp.norm_prm(k, :), stim', fBasis_long, 140 * 7 + 1 : 140 * 13));
    % compute r2
    gcamp.norm_r2(k) = compute_r2(gcamp.tofit(k, :), gcamp.norm_prd(k, :));

end

%% visualize

% visualizing gevi normalization fit
figure (2),
subplot(2, 1, 2)
shadedErrorBar(t_plot, mean(gevi.tofit), std(gevi.tofit)), hold on
plot(t_plot, mean(gevi.norm_prd), 'r-'), box off
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

% visualizing gcamp normalization fit
figure (3),
subplot(2, 1, 2)
shadedErrorBar(t_plot, mean(gcamp.tofit), std(gcamp.tofit)), hold on
plot(t_plot, mean(gcamp.norm_prd), 'r-'), box off
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

%% GEVI and GCaMP filters

% load vsd data
fileLoc = '/Users/jyzhou/Desktop/GEVI/';
b = load(fullfile(fileLoc, 'vsd_parameters.mat'));
vsd = b.vsd;

% make vsd linear filters --------------------------------
vsd.lin_filters = vsd.lin_w * fBasis;
for k = 1 : n_boots
    vsd.lin_filters(k, :) = norm_max(vsd.lin_filters(k, :));
end


% Linear filters (GEVI and VSD) --------------------------
gevi.lin_filters   = gevi.lin_w * fBasis;
% compute the numerator filter for the normalization model
gevi.norm_filters1 = gevi.norm_prm(:, 1 : nFast) * fBasis;
% normalize the numerator filter to the maximum value.
for k = 1 : n_boots
    gevi.lin_filters(k, :) = norm_max(gevi.lin_filters(k, :));
    gevi.norm_filters1(k, :) = norm_max(gevi.norm_filters1(k, :));
end

% compute the denominator filter for the normalization model
for k = 1 : n_boots
    gevi.norm_filters2(k, :) = gammaPDF(t, gevi.norm_prm(k, nFast + 1), gevi.norm_prm(k, nFast + 2)) - ...
        gevi.norm_prm(k, nFast + 4) * gammaPDF(t, gevi.norm_prm(k, nFast + 5), gevi.norm_prm(k, nFast + 6));
    gevi.norm_filters2(k, :) = norm_max([0, gevi.norm_filters2(k, 1 : end - 1)]);
end

% exclude an outlier
gevi.norm_filters2(42, :) = nan(1, length(gevi.norm_filters2));

figure (4), clf
subplot(131), shadedErrorBar(t(1 : f_lth), mean(gevi.lin_filters), std(gevi.lin_filters)), axis square, box off, hold on
shadedErrorBar(t(1 : f_lth), mean(vsd.lin_filters), std(vsd.lin_filters), 'b')
plot([t(1), t(f_lth)], [0, 0], 'k:'), title('GEVI filters')

% normalization filters (GEVI) ---------------------------
subplot(132), shadedErrorBar(t(1 : f_lth), mean(gevi.norm_filters1), std(gevi.norm_filters1)), axis square, box off, hold on
plot([t(1), t(f_lth)], [0, 0], 'k:')
subplot(133), %plot(t, gevi.norm_filters2), axis square, box off,
shadedErrorBar(t, nanmean(gevi.norm_filters2), nanstd(gevi.norm_filters2)), axis square, box off, hold on
xlim([t(1), t(f_lth)]), plot([t(1), t(f_lth)], [0, 0], 'k:')

% Linear filters (GEVI and GCaMP) ------------------------
gcamp.lin_filters = gcamp.lin_w * fBasis;
% compute the numerator filter for the normalization model
gcamp.norm_filters1 = gcamp.norm_prm(:, 1 : nFast) * fBasis;

for k = 1 : n_boots
    gcamp.lin_filters(k, :) = norm_max(gcamp.lin_filters(k, :));
    gcamp.norm_filters1(k, :) = norm_max(gcamp.norm_filters1(k, :));
end

% compute the denominator filter for the normalization model
for k = 1 : n_boots
    gcamp.norm_filters2(k, :) = real(gammaPDF(t, gcamp.norm_prm(k, nFast + 1), gcamp.norm_prm(k, nFast + 2)) - ...
        gcamp.norm_prm(k, nFast + 4) * gammaPDF(t, gcamp.norm_prm(k, nFast + 5), gcamp.norm_prm(k, nFast + 6)));
    gcamp.norm_filters2(k, :) = norm_max([0, gcamp.norm_filters2(k, 1 : end - 1)]);
end

figure (5), clf
subplot(131), shadedErrorBar(t(1 : f_lth), mean(gevi.lin_filters), std(gevi.lin_filters)), axis square, box off, hold on
shadedErrorBar(t(1 : f_lth), mean(gcamp.lin_filters), std(gcamp.lin_filters)), axis square, box off, hold on
plot([t(1), t(f_lth)], [0, 0], 'k:'), title('GCaMP filters')

% normalization filters (GCaMP) --------------------------
subplot(132), shadedErrorBar(t(1 : f_lth), mean(gcamp.norm_filters1), std(gcamp.norm_filters1)), hold on, axis square, box off
plot([t(1), t(f_lth)], [0, 0], 'k:')
subplot(133), %plot(t, gcamp.norm_filters2), axis square, box off,
shadedErrorBar(t, mean(gcamp.norm_filters2), std(gcamp.norm_filters2)), axis square, box off, hold on
xlim([t(1), t(f_lth)]), plot([t(1), t(f_lth)], [0, 0], 'k:')

%% Temporal frequency comparison

% Here, I want to compare temporal frequencies for GEVI, VSD and GCaMP

gevi.frq_amp = []; vsd.frq_amp = []; gcamp.frq_amp = [];

freq = [2, 4, 8.3, 10, 16.7, 20, 33.3];

fs = 100; % sampling frequency

for iboot = 1 : n_boots
    for k = 1 : length(freq)
        rng = (k - 1) * 140 + 1 : k * 140;
        % vsd -------------------------
        x1 = vsd.tofit(iboot, rng);
        N1 = length(x1);
        y1 = fft(x1);

        f1 = (0:N1-1)*(fs/N1);
        p  = abs(y1/N1);

        % find frequency
        sampling_frq = f1(1 : N1/2);
        % find the most
        [~, idx] = min(abs(sampling_frq - freq(k)));
        vsd.frq_amp(iboot, k) = p(idx)/p(1);

        % gevi ------------------------------------
        x2 = gevi.tofit(iboot, rng);
        N2 = length(x2);
        y2 = fft(x2);

        p2 = abs(y2 / N2);
        gevi.frq_amp(iboot, k) = p2(idx)/p2(1);

        % gcamp ------------------------------------
        x3 = gcamp.tofit(iboot, rng);
        N3 = length(x3);
        y3 = fft(x3);

        p3 = abs(y3 / N3);
        gcamp.frq_amp(iboot, k) = p3(idx)/p3(1);
    end
    vsd.frq_amp(iboot, :) = norm_max(vsd.frq_amp(iboot, :));
    gevi.frq_amp(iboot, :) = norm_max(gevi.frq_amp(iboot, :));
    gcamp.frq_amp(iboot, :) = norm_max(gcamp.frq_amp(iboot, :));
end

figure (11), 

shadedErrorBar(freq, mean(vsd.frq_amp), std(vsd.frq_amp), 'k'), hold on
shadedErrorBar(freq, mean(gevi.frq_amp), std(gevi.frq_amp), 'r')
shadedErrorBar(freq, mean(gcamp.frq_amp), std(gcamp.frq_amp), 'b'), 
axis square, set(gca, 'xtick', freq), box off

%% Temporal frequency comparison (using model fits)


gevi.frq_amp = []; vsd.frq_amp = []; gcamp.frq_amp = [];

freq = [2, 4, 8.3, 10, 16.7, 20, 33.3];

fs = 100; % sampling frequency

for iboot = 1 : n_boots
    for k = 1 : length(freq)
        rng = (k - 1) * 140 + 1 : k * 140;
        % vsd -------------------------
        x1 = vsd.norm_prd(iboot, rng);
        N1 = length(x1);
        y1 = fft(x1);

        f1 = (0:N1-1)*(fs/N1);
        p  = abs(y1/N1);

        % find frequency
        sampling_frq = f1(1 : N1/2);
        % find the most
        [~, idx] = min(abs(sampling_frq - freq(k)));
        vsd.frq_amp(iboot, k) = p(idx)/p(1);

        % gevi ------------------------------------
        x2 = gevi.norm_prd(iboot, rng);
        N2 = length(x2);
        y2 = fft(x2);

        p2 = abs(y2 / N2);
        gevi.frq_amp(iboot, k) = p2(idx)/p2(1);

        % gcamp ------------------------------------
        x3 = gcamp.norm_prd(iboot, rng);
        N3 = length(x3);
        y3 = fft(x3);

        p3 = abs(y3 / N3);
        gcamp.frq_amp(iboot, k) = p3(idx)/p3(1);
    end
    vsd.frq_amp(iboot, :) = norm_max(vsd.frq_amp(iboot, :));
    gevi.frq_amp(iboot, :) = norm_max(gevi.frq_amp(iboot, :));
    gcamp.frq_amp(iboot, :) = norm_max(gcamp.frq_amp(iboot, :));
end

figure (12), 

shadedErrorBar(freq, mean(vsd.frq_amp), std(vsd.frq_amp), 'k'), hold on
shadedErrorBar(freq, mean(gevi.frq_amp), std(gevi.frq_amp), 'r')
shadedErrorBar(freq, mean(gcamp.frq_amp), std(gcamp.frq_amp), 'b'), 
axis square, set(gca, 'xtick', freq), box off


%% contrast response functions

ctr_levels = [3.125, 6.25, 12.5, 25, 50, 100];

n_ctrs = 6;

% extracting contrast responses from gevi and gamp data
gcamp_ctr = []; gevi_ctr = [];
for k = 1 : n_ctrs
    idx = 140 * (n_ctrs + k) + 1 : 140 * (n_ctrs + k + 1);
    gcamp_ctr(k, :, :) = gcamp.tofit(:, idx);
    gevi_ctr(k, :, :) = gevi.tofit(:, idx);
end
s_gcamp_ctr = flipud(sum(gcamp_ctr, 3));
s_gevi_ctr  = flipud(sum(gevi_ctr, 3));

for k = 1 : 50
     s_gcamp_ctr(:, k) = norm_max(s_gcamp_ctr(:, k));
     s_gevi_ctr(:, k) = norm_max(s_gevi_ctr(:, k));
end

% visualize contrast response functions
figure (6), clf
subplot(131), plot(ctr_levels,  median(s_gevi_ctr, 2), 'ko:', 'markerfacecolor', 'k'), axis square, hold on
% plot error bar
for k = 1 : 6
    plot([ctr_levels(k), ctr_levels(k)], [prctile(s_gevi_ctr(k, :), 25), prctile(s_gevi_ctr(k, :), 75)], 'k-', 'linewidth', 1)
end
set(gca, 'xtick', ctr_levels), box off, set(gca, 'fontsize', 14)
title('GEVI contrast responses')

subplot(132), plot(ctr_levels, median(s_gcamp_ctr, 2), 'ko:', 'markerfacecolor', 'k'), axis square, hold on
for k = 1 : 6
    plot([ctr_levels(k), ctr_levels(k)], [prctile(s_gcamp_ctr(k, :), 25), prctile(s_gcamp_ctr(k, :), 75)], 'k-', 'linewidth', 1)
end
set(gca, 'xtick', ctr_levels), box off, set(gca, 'fontsize', 14)
title('GCaMP contrast responses')

%% transform between contrast response functions

ctr_prm = [];

s_trans = @(dt, prm) prm(2).* dt.^prm(1);
ls = @(dt, prm, dt2) sum((dt2 - s_trans(dt, prm)).^2);

init = [1, 1];

for k = 1 : 50
    ctr_prm(k, :) = fminsearch(@(x) ls(s_gevi_ctr(:, k), x, s_gcamp_ctr(:, k)), init);
    ctr_power(k, :) = norm_max(ctr_prm(k, 2) * s_gevi_ctr(:, k).^ctr_prm(k, 1));

end

% show the predicted power transform
figure (6), subplot(133), cla
for k = 1 : 50
    %x = [0 : 0.01: 1];
    %plot(x, x.^ctr_prm(k, 1)), hold on, axis square
    plot(s_gevi_ctr(:, k), ctr_power(k, :), 'r'), hold on, axis square
end

plot(mean(s_gevi_ctr, 2), mean(ctr_power), 'b-', 'linewidth', 2)

plot(mean(s_gevi_ctr, 2), mean(s_gcamp_ctr, 2), 'ko', 'markerfacecolor', 'k')
box off
xlabel('normalized GEVI contrast response')
ylabel('predicted gcamp contrast response')
set(gca, 'fontsize', 14)

%% do the overall transforms

trans_prm = []; trans_prd = []; trans_weights = [];

nTrans = 5;

tBasis3 = mkBasis(t(1 : 40), nTrans, 'fast'); %32

%prm = fminsearch(@(x) fit_gevi2gcamp(m_gevi_dt, x, tBasis3, m_gcamp_dt), [1, 1]); %[3, 0, 1, 0]

% make prediction
%[trans_prd, trans_weights] = gevi2gcamp_fun(m_gevi_dt, prm, tBasis3, m_gcamp_dt);

for k = 1 : n_boots
    trans_prm(k, :) = fminsearch(@(x) fit_gevi2gcamp(gevi.tofit(k, :), x, tBasis3, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13), [1, 1, 1]);
    [trans_prd(k, :), trans_weights(:, k)] = gevi2gcamp_fun(gevi.tofit(k, :), trans_prm(k, :) , tBasis3, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13);
    trans_r2(k) = compute_r2(gcamp.tofit(k, :), trans_prd(k, :));

end

% compute r2 separately for temporal frequency and contrast
% for k = 1 : n_boots
%     transr2_tf(k) = compute_r2(gcamp.tofit(k, 1 : 140 * 7), trans_prd(k, 1 : 140 * 7));
%     transr2_ctr(k) = compute_r2(gcamp.tofit(k, 140 * 7 + 1 : end), trans_prd(k, 140 * 7 + 1 : end));
% end

% do the overall transform, but separately for temporal frequency and for
% contrast
for k = 1 : n_boots
    trans_prm_tf(k, :) = fminsearch(@(x) fit_gevi2gcamp(gevi.tofit(k, 1 : 140 * 7), x, tBasis3, gcamp.tofit(k, 1 : 140 * 7), 1 : 140 * 7), [1, 1, 1]);
    trans_prm_ctr(k, :) = fminsearch(@(x) fit_gevi2gcamp(gevi.tofit(k, 140 * 7 + 1 : end), x, tBasis3, gcamp.tofit(k, 140 * 7 + 1 : end), 1 : 140 * 6), [1, 1, 1]);

end

% visualize

figure (7), clf
shadedErrorBar(t_plot, mean(gcamp.tofit), std(gcamp.tofit)), hold on
shadedErrorBar(t_plot, mean(trans_prd), std(trans_prd), 'r-'), box off
%set(gca, 'xtick', 0 : 140 : 140 * 13)
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

figure (8), clf % filters
for k = 1 : 50
    trans_filter(k, :) = norm_max(trans_weights(:, k)' * tBasis3);
end

trans_filter(isinf(trans_filter)) = nan;

shadedErrorBar(t_plot(1 : size(trans_filter, 2)), nanmean(trans_filter), nanstd(trans_filter))

% for k = 1 : 50
%     plot(0.01 : 0.01 : 0.01 * 40, norm_max(trans_weights(:, k)' * tBasis3)), hold on
% end
axis square, box off

%% do the overall transform (from model instead of from data)

trans_prm = []; trans_prd = []; trans_weights = [];

nTrans = 5;

tBasis3 = mkBasis(t(1 : 40), nTrans, 'fast'); %32

%prm = fminsearch(@(x) fit_gevi2gcamp(m_gevi_dt, x, tBasis3, m_gcamp_dt), [1, 1]); %[3, 0, 1, 0]

% make prediction
%[trans_prd, trans_weights] = gevi2gcamp_fun(m_gevi_dt, prm, tBasis3, m_gcamp_dt);

for k = 1 : n_boots
    trans_prm(k, :) = fminsearch(@(x) fit_gevi2gcamp(gevi.norm_prd(k, :), x, tBasis3, gcamp.norm_prd(k, :), 140 * 7 + 1 : 140 * 13), [1, 1, 1]);
    [trans_prd(k, :), trans_weights(:, k)] = gevi2gcamp_fun(gevi.norm_prd(k, :), trans_prm(k, :) , tBasis3, gcamp.norm_prd(k, :), 140 * 7 + 1 : 140 * 13);
    trans_r2(k) = compute_r2(gcamp.tofit(k, :), trans_prd(k, :));

end

% visualize

figure (7), clf
shadedErrorBar(t_plot, mean(gcamp.norm_prd), std(gcamp.norm_prd)), hold on
shadedErrorBar(t_plot, mean(trans_prd), std(trans_prd), 'r-'), box off
%set(gca, 'xtick', 0 : 140 : 140 * 13)
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

figure (8), clf % filters
for k = 1 : 50
    trans_filter(k, :) = norm_max(trans_weights(:, k)' * tBasis3);
end

trans_filter(isinf(trans_filter)) = nan;

shadedErrorBar(t_plot(1 : size(trans_filter, 2)), nanmean(trans_filter), nanstd(trans_filter))

% for k = 1 : 50
%     plot(0.01 : 0.01 : 0.01 * 40, norm_max(trans_weights(:, k)' * tBasis3)), hold on
% end
axis square, box off

%% do the overall transform with fixed power

trans_prm1 = []; trans_prd1 = []; trans_weights1 = [];


for k = 1 : 50
    trans_prm1(k, :) = fminsearch(@(x) fit_gevi2gcamp1(gevi.tofit(k, :), x, ctr_prm(k, 1), tBasis3, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13), [1, 1]);
    [trans_prd1(k, :), trans_weights1(:, k)] = gevi2gcamp_fun1(gevi.tofit(k, :), trans_prm1(k, :) , ctr_prm(k, 1), tBasis3, gcamp.tofit(k, :), 140 * 7 + 1 : 140 * 13);
    trans_r21(k) = compute_r2(gcamp.tofit(k, :), trans_prd1(k, :));
end


% transform 


figure (9), clf
shadedErrorBar(t_plot, mean(gcamp.tofit), std(gcamp.tofit)), hold on
shadedErrorBar(t_plot, mean(trans_prd1), std(trans_prd1), 'r-'), box off
xlim([1.4 * 7,  1.4 * 13])
%set(gca, 'xtick', 0 : 140 : 140 * 13)
%set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

figure (10), clf % filters

for k = 1 : 50
    trans_filter1(k, :) = norm_max(trans_weights1(:, k)' * tBasis3);
end

trans_filter1(isinf(trans_filter1)) = nan;

shadedErrorBar(t_plot(1 : size(trans_filter1, 2)), nanmean(trans_filter1), nanstd(trans_filter1)), hold on

%shadedErrorBar(t_plot(1 : size(trans_filter, 2)), nanmean(trans_filter), nanstd(trans_filter), 'r-')

axis square, box off

figure (11), clf
plot(1, mean(trans_r2), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), hold on
plot(2, mean(trans_r21), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), 
plot([1, 1], [mean(trans_r2) - std(trans_r2), mean(trans_r2) + std(trans_r2)], 'linewidth', 1), 
plot([2, 2], [mean(trans_r21) - std(trans_r21), mean(trans_r21) + std(trans_r21)], 'linewidth', 1)
xlim([0.5, 2.5]), ylim([0.75, 1]), grid on, box off, axis square

%% transform between contrast part of the data with fixed power


trans_prm1 = []; trans_prd1 = []; trans_weights1 = [];


for k = 1 : 50
    trans_prm1(k, :) = fminsearch(@(x) fit_gevi2gcamp1(gevi.tofit(k, 140 *7 + 1 : 140 * 13), x, ctr_prm(k, 1), tBasis3, gcamp.tofit(k, 140 * 7 +1: 140 * 13), 140 * 1 : 140 * 6), [1, 1]);
    [trans_prd1(k, :), trans_weights1(:, k)] = gevi2gcamp_fun1(gevi.tofit(k, 140 *7 + 1 : 140 * 13), trans_prm1(k, :) , ctr_prm(k, 1), tBasis3, gcamp.tofit(k, 140 *7 + 1 : 140 * 13), 140 * 1 : 140 * 6);
    trans_r21(k) = compute_r2(gcamp.tofit(k, 140 *7 + 1 : 140 * 13), trans_prd1(k, :));
end


% transform 


figure (9), clf
shadedErrorBar(t_plot(140 * 7 + 1 : 140 * 13), mean(gcamp.tofit(:, 140 *7 + 1 : 140 * 13)), std(gcamp.tofit(:, 140 *7 + 1 : 140 * 13))), hold on
shadedErrorBar(t_plot(140 * 7 + 1 : 140 * 13), mean(trans_prd1), std(trans_prd1), 'r-'), box off
%xlim([1.4 * 0,  1.4 * 6])
set(gca, 'xtick', 0 : 140 : 140 * 13)
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

figure (10), clf % filters

for k = 1 : 50
    trans_filter1(k, :) = norm_max(trans_weights1(:, k)' * tBasis3);
end

trans_filter1(isinf(trans_filter1)) = nan;

shadedErrorBar(t_plot(1 : size(trans_filter1, 2)), nanmean(trans_filter1), nanstd(trans_filter1)), hold on

%shadedErrorBar(t_plot(1 : size(trans_filter, 2)), nanmean(trans_filter), nanstd(trans_filter), 'r-')

axis square, box off

figure (11), clf
plot(1, mean(trans_r2), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), hold on
plot(2, mean(trans_r21), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), 
plot([1, 1], [mean(trans_r2) - std(trans_r2), mean(trans_r2) + std(trans_r2)], 'linewidth', 1), 
plot([2, 2], [mean(trans_r21) - std(trans_r21), mean(trans_r21) + std(trans_r21)], 'linewidth', 1)
xlim([0.5, 2.5]), ylim([0.75, 1]), grid on, box off, axis square


%% do the overall transform with fixed power (between model and model)

trans_prm1 = []; trans_prd1 = []; trans_weights1 = [];


for k = 1 : 50
    trans_prm1(k, :) = fminsearch(@(x) fit_gevi2gcamp1(gevi.norm_prd(k, :), x, ctr_prm(k, 1), tBasis3, gcamp.norm_prd(k, :), 140 * 7 + 1 : 140 * 13), [1, 1]);
    [trans_prd1(k, :), trans_weights1(:, k)] = gevi2gcamp_fun1(gevi.norm_prd(k, :), trans_prm1(k, :) , ctr_prm(k, 1), tBasis3, gcamp.norm_prd(k, :), 140 * 7 + 1 : 140 * 13);
    trans_r21(k) = compute_r2(gcamp.tofit(k, :), trans_prd1(k, :));
end

figure (9), clf
shadedErrorBar(t_plot, mean(gcamp.norm_prd), std(gcamp.norm_prd)), hold on
shadedErrorBar(t_plot, mean(trans_prd1), std(trans_prd1), 'r-'), box off
%set(gca, 'xtick', 0 : 140 : 140 * 13)
set(gca, 'xtick', t_plot(140 : 140 : end), 'XTickLabel', 1.4 * ones(1, 13))

figure (10), clf % filters

for k = 1 : 50
    trans_filter1(k, :) = norm_max(trans_weights1(:, k)' * tBasis3);
end

trans_filter1(isinf(trans_filter1)) = nan;

shadedErrorBar(t_plot(1 : size(trans_filter1, 2)), nanmean(trans_filter1), nanstd(trans_filter1)), hold on

%shadedErrorBar(t_plot(1 : size(trans_filter, 2)), nanmean(trans_filter), nanstd(trans_filter), 'r-')

axis square, box off

figure (11), clf
plot(1, mean(trans_r2), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), hold on
plot(2, mean(trans_r21), 'ko', 'markersize', 4, 'markerfacecolor', 'k'), 
plot([1, 1], [mean(trans_r2) - std(trans_r2), mean(trans_r2) + std(trans_r2)], 'linewidth', 1), 
plot([2, 2], [mean(trans_r21) - std(trans_r21), mean(trans_r21) + std(trans_r21)], 'linewidth', 1)
xlim([0.5, 2.5]), ylim([0.75, 1]), grid on, box off, axis square

%% Additional analysis to understand the transform

%% understand the filter

% first, fit two gamma functions to the filters
mfilter = nanmean(trans_filter1);

%init = [0.0014, 2, 0.0034, 70, 1, .1];
init = [0.0013, 2, 0.0034, 70, 1, .1];
filter_prm = fminsearch(@(x) fit_filterModel(t, x, mfilter, length(mfilter)), init);
filter_prd = filterModel(t, filter_prm);

figure (13), clf
plot(mfilter, 'k-'), hold on
plot(filter_prd(1 : 40), 'r-')


%% make simulations

% make simulation in steps:

k = 2;
% First the biphasic filter
f_biphasic = trans_weights(1 : 2, k)' * tBasis3(1 : 2, :);
mgevi = mean(gevi.tofit);

nl_mgevi = real(mgevi.^2.231);

figure (14), clf

subplot(4, 3, 4), plot(linspace(0, 1, 100), linspace(0, 1, 100).^2.231, 'k-'), axis square, box off
subplot(4, 3, 5 : 6), plot(nl_mgevi), box off
set(gca, 'xtick', 0 : 140 : 140 * 13)

subplot(437), plot(t(1 : 40), f_biphasic(1 : 40), 'k-'), box off, axis square, %axis square %xlim([0, 0.4])


subplot(4, 3, 2 : 3 ), plot(mgevi), box off, 
set(gca, 'xtick', 0 : 140 : 140 * 13)
subplot(4, 3, 8 : 9), plot(convCut(nl_mgevi, filter1)), box off, 


% Second the slow filter
f_slow = trans_weights(:, k)' * tBasis3;

subplot(4, 3, 10), plot(t(1 : length(f_slow)), f_slow, 'r-'), box off, axis square

hold on, plot(t(1 : length(f_slow)), trans_weights(3 : 5, k)' * tBasis3(3 : 5, 1 : 40), 'k--')
subplot(4, 3, 11 : 12), plot(convCut(nl_mgevi, f_slow)), box off

% then the exponential component
figure (15), clf

x = convCut(nl_mgevi, f_slow);
y = mean(gcamp.tofit);

scale = least_square(x', y');

plot(scale * convCut(nl_mgevi, f_slow)), hold on, 
plot(mean(gcamp.tofit), 'r-')



% % first make just the first bump's prediction
% filter1 = gammaPDF(t, filter_prm(1), filter_prm(2));
% filter1 = [0, filter1];
% figure (14), clf
% subplot(334), plot(t(1 : 40), filter1(1 : 40)), box off, %axis square %xlim([0, 0.4])
% 
% mgevi = mean(gevi.tofit);
% subplot(3, 3, 2 : 3), plot(mgevi), box off, 
% subplot(3, 3, 5 : 6), plot(convCut(mgevi, filter1)), box off, 
% 
% subplot(3, 3, 7), plot(t(1 : 40), filter_prd(1 : 40)), box off, %axis square
% 
% subplot(3, 3, 8 : 9), plot(convCut(mgevi, filter_prd)), box off

% now make the combined prediction







