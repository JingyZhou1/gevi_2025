

figureLoc = '/Users/jyzhou/Desktop/GEVI/figures/';

%% load data

a = load('preprocessed_data.mat');

gevi = a.gevi;
gcamp = a.gcamp;

%% concatenate data

% make gevi fitting data
gevi.tofit = cat(3, (gevi.bs_normtf(:, :, 1 : 7) + gevi.bs_normtf(:, :, 8 : 14))/2, (gevi.bs_normctr(:, :, 1 : 6) + gevi.bs_normctr(:, :, 7 : 12))/2);

gevi.tofit = reshape(gevi.tofit, 50, []);

% make gcamp fitting data
gcamp.tofit = cat(3, (gcamp.bs_normtf(:, :, 1 : 7) + gcamp.bs_normtf(:, :, 8 : 14))/2, (gcamp.bs_normctr(:, :, 1 : 6) + gcamp.bs_normctr(:, :, 7 : 12))/2);

gcamp.tofit = reshape(gcamp.tofit, 50, []);

figure (1), clf
subplot(211), plot(gevi.tofit')
subplot(212), plot(gcamp.tofit')

%% prepping for model fit

% make stimulus
stim1 = mkStim('TF', 140);
stim1 = stim1(:, 3 : 9);

stim2 = mkStim('contrast', 140);
stim2 = stim2(:, 3 : 8);

stim = reshape(cat(2, stim1, stim2), [], 1);

t_lth = 140;

% % make basis functions
t = linspace(0.01, t_lth/100, t_lth);
%
nFast = 6;
nSlow = 5;

f_lth = 50;

fBasis  = mkBasis(t(1 : f_lth), nFast, 'fast_longdt');
sBasis  = mkBasis(t, nSlow, 'slow');

% make basis function -------------------------------------
basis = concatenateBasisAcrossConditions(fBasis, sBasis, stim, t);

%% fit linear model

gevi.lin_w = []; gcamp.lin_w = []; gevi.norm_prm = []; gcamp.norm_prm = [];

for k = 1 : 50
    k
    % fit GEVI MODEL
    % ------------------------------------------------------------
    % fit linear model
    gevi.lin_w(k, :) = least_square(basis', gevi.tofit(k, :)');
    gevi_fast_left(k, :) = gevi.tofit(k, :)- gevi.lin_w(k, 1 : nSlow * 13) * basis(1 : nSlow * 13, :);

    % fit normalization model
    %init = [gevi.lin_w(k, 13 * nSlow+ 1 : end), 0.03, 4, .1, .6, 0.05, 5, 1];
    init = [gevi.lin_w(k, 13 * nSlow + 1 : end), 0.03, 4, .1, .6, 0.007, 20, 1];%

    for iter = 1 : 30
        if iter == 1
            gevi.norm_prm(k, :) = init;
            fast_left = gevi_fast_left(k, :);
        end
        gevi.norm_prm(k, :) = fminsearch(@(x) fitModel1(t, x, stim',  basis(13 * nSlow+ 1 : end, :), fast_left, 140 * 7 + 1 : 140 * 13), gevi.norm_prm(k, :));
        % fit slow component
        fastprd   = real(normModel3(t, gevi.norm_prm(k, :), stim', basis(13 * nSlow  + 1 : end, :), 140 * 7 + 1 : 140 * 13));
        slow_w    = least_square(basis(1 : 13 * nSlow, :)', (gevi.tofit(k, :) - fastprd)');
        slowprd   = slow_w' * basis(1 : 13 * nSlow, :);
        fast_left = gevi.tofit(k, :) - slowprd;
    end
    gevi.norm_fastleft(k, :) = fast_left';
    gevi.norm_fastprd(k, :) = fastprd';

    % FIT GCAMP MODEL ----------------------------------------
    % fit linear model
    gcamp.lin_w(k, :) = least_square(basis', gcamp.tofit(k, :)');
    gcamp_fast_left(k, :) = gcamp.tofit(k, :) - gcamp.lin_w(k, 1 : nSlow * 13) * basis(1 : nSlow * 13, :);

    % fit normalization model
    init = [gcamp.lin_w(k, 13 * nSlow+ 1 : end), 0.03, 4, .1, .6, 0.05, 5, 1];

    for iter = 1 : 30
        if iter == 1
            gcamp.norm_prm(k, :) = init;
            fast_left = gcamp_fast_left(k, :);
        end
        gcamp.norm_prm(k, :) = fminsearch(@(x) fitModel1(t, x, stim',  basis(13 * nSlow+ 1 : end, :), fast_left, 140 * 7 + 1 : 140 * 13), gcamp.norm_prm(k, :));
        % fit slow component
        fastprd   = real(normModel3(t, gcamp.norm_prm(k, :), stim', basis(13 * nSlow  + 1 : end, :), 140 * 7 + 1 : 140 * 13));
        slow_w    = least_square(basis(1 : 13 * nSlow, :)', (gcamp.tofit(k, :) - fastprd)');
        slowprd   = slow_w' * basis(1 : 13 * nSlow, :);
        fast_left = gcamp.tofit(k, :) - slowprd;
    end
    gcamp.norm_fastleft(k, :) = fast_left';
    gcamp.norm_fastprd(k, :) = fastprd';
end

%% visualization gevi


m_gevi_dt = mean(gevi.norm_fastleft); s_gevi_dt = std(gevi.norm_fastleft);

% linear model fit
gevi_lin_fast_left = gevi.tofit- gevi.lin_w(:, 1 : nSlow * 13) * basis(1 : nSlow * 13, :);

m_gevi_lindt = mean(gevi_lin_fast_left); s_gevi_lindt = std(gevi_lin_fast_left);

gevi_linprd = gevi.lin_w(:, 13 * nSlow + 1 : end) * basis(nSlow * 13 + 1 : end, :);

t_plot = [1 : length(m_gevi_dt)]./100;

figure (2), clf
subplot(211), 
shadedErrorBar(t_plot, m_gevi_lindt, s_gevi_lindt), hold on
plot(t_plot, gevi_linprd, 'r-')

subplot(212)
shadedErrorBar(t_plot, m_gevi_dt, s_gevi_dt), hold on
plot(t_plot, mean(gevi.norm_fastprd), 'r-')

%% visualization gcamp


m_gcamp_dt = mean(gcamp.norm_fastleft); s_gcamp_dt = std(gcamp.norm_fastleft);

% linear model fit
gcamp_lin_fast_left = gcamp.tofit- gcamp.lin_w(:, 1 : nSlow * 13) * basis(1 : nSlow * 13, :);

m_gcamp_lindt = mean(gcamp_lin_fast_left); s_gcamp_lindt = std(gcamp_lin_fast_left);

gcamp_linprd = gcamp.lin_w(:, 13 * nSlow + 1 : end) * basis(nSlow * 13 + 1 : end, :);

t_plot = [1 : length(m_gcamp_dt)]./100;

figure (3), clf
subplot(211), 
shadedErrorBar(t_plot, m_gcamp_lindt, s_gcamp_lindt), hold on, box off
set(gca, 'xtick', 0 : 1.4 : 1.4 * 13)
plot(t_plot, gcamp_linprd, 'r-')
set(gca, 'xticklabel', ones(1, 13) *1.4)

subplot(212)
shadedErrorBar(t_plot, m_gcamp_dt, s_gcamp_dt), hold on, box off
plot(t_plot, mean(gcamp.norm_fastprd), 'r-')
set(gca, 'xtick', 0 : 1.4 : 1.4 * 13)
set(gca, 'xticklabel', ones(1, 13) *1.4)

%% do the transform

nTrans = 5;

tBasis3 = mkBasis(t(1 : 40), nTrans, 'fast'); %32

% prm = fminsearch(@(x) fit_gevi2gcamp(m_gevi_dt, x, tBasis3, m_gcamp_dt), [1, 1]); %[3, 0, 1, 0]
% 
% % make prediction
% [trans_prd, trans_weights] = gevi2gcamp_fun(m_gevi_dt, prm, tBasis3, m_gcamp_dt);
% 
% 
% figure (4), clf
% plot(m_gcamp_dt, 'k-'), hold on
% plot(trans_prd, 'r-')
% 
% % filter shape
% figure (5), clf
% 
% plot(trans_weights' * tBasis3)

%% contrast response function bootstrap

ctr_levels = [3.125, 6.25, 12.5, 25, 50, 100]

% extracting contrast responses from gevi and gamp data
gcamp_ctr = []; gevi_ctr = [];
for k = 1 : 6
    idx = 140 * (7 + k - 1) + 1 : 140 * (7 + k); 
    gcamp_ctr(k, :, :) = gcamp.norm_fastleft(:, idx);
    gevi_ctr(k, :, :) = gevi.norm_fastleft(:, idx);
end
s_gcamp_ctr = flipud(sum(gcamp_ctr, 3));
s_gevi_ctr  = flipud(sum(gevi_ctr, 3));

% plotting contrast response functions
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

%% contrast transform

s_trans = @(dt, prm) prm(2).* dt.^prm(1);
ls = @(dt, prm, dt2) sum((dt2 - s_trans(dt, prm)).^2);

init = [1, 1];

for k = 1 : 50
    ctr_prm(k, :) = fminsearch(@(x) ls(s_gevi_ctr(:, k), x, s_gcamp_ctr(:, k)), init);
end

% show the predicted power transform
figure (6), subplot(133), 
for k = 1 : 50
    x = [0 : 0.01: 1];
    plot(x, x.^ctr_prm(k, 1)), hold on, axis square
end
box off
xlabel('normalized GEVI contrast response')
ylabel('predicted gcamp contrast response')
set(gca, 'fontsize', 14)

%% do the overall transforms

trans_prm = []; trans_prd = []; trans_weights = [];

%prm = fminsearch(@(x) fit_gevi2gcamp(m_gevi_dt, x, tBasis3, m_gcamp_dt), [1, 1]); %[3, 0, 1, 0]

% make prediction
%[trans_prd, trans_weights] = gevi2gcamp_fun(m_gevi_dt, prm, tBasis3, m_gcamp_dt);

for k = 1 : 50
    trans_prm(k, :) = fminsearch(@(x) fit_gevi2gcamp(gevi.norm_fastleft(k, :), x, tBasis3, gcamp.norm_fastleft(k, :)), [1, 1]);
    [trans_prd(k, :), trans_weights(:, k)] = gevi2gcamp_fun(gevi.norm_fastleft(k, :), trans_prm(k, :) , tBasis3, gcamp.norm_fastleft(k, :));
end

%% visualize

figure (7), clf
shadedErrorBar(1 : 140 * 13, mean(gcamp.norm_fastleft), std(gcamp.norm_fastleft)), hold on
plot(1 : 140 * 13, mean(trans_prd), 'r-'), box off
set(gca, 'xtick', 0 : 140 : 140 * 13)

figure (8), clf % filters
for k = 1 : 50
    plot(0.01 : 0.01 : 0.01 * 40, norm_max(trans_weights(:, k)' * tBasis3)), hold on
end
axis square, box off

%% do the overall transform with fixed power

trans_prm1 = []; trans_prd1 = []; trans_weights1 = [];


for k = 1 : 50
    trans_prm1(k, :) = fminsearch(@(x) fit_gevi2gcamp1(gevi.norm_fastleft(k, :), x, ctr_prm(k, 1), tBasis3, gcamp.norm_fastleft(k, :)), [1, 1]);
    [trans_prd1(k, :), trans_weights1(:, k)] = gevi2gcamp_fun1(gevi.norm_fastleft(k, :), trans_prm(k, :) , ctr_prm(k, 1), tBasis3, gcamp.norm_fastleft(k, :));
end

%% visualize

figure (9), clf
shadedErrorBar(1 : 140 * 13, mean(gcamp.norm_fastleft), std(gcamp.norm_fastleft)), hold on
plot(1 : 140 * 13, mean(trans_prd1), 'r-'), box off
set(gca, 'xtick', 0 : 140 : 140 * 13)

figure (10), clf % filters
for k = 1 : 50
    plot(0.01 : 0.01 : 0.01 * 40, norm_max(trans_weights1(:, k)' * tBasis3)), hold on
end
axis square, box off








