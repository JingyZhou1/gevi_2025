% gevi_fitModel

% load data
% fit linear and normalization model for pre-processing

% extract the fast component for each data type, then take an average and
% bootstrap.

figLoc = '/Users/jyzhou/Library/CloudStorage/GoogleDrive-jyz205@nyu.edu/My Drive/gd_Projects/GEVI/figures';

%% load gevi tf data

gevi = [];

gevi.tf = []; gevi.mtf = [];, gevi.lintf = []; gevi.normtf = []; gevi.bs_normtf = [];

[gevi.tf, gevi.mtf, gevi.lintf, gevi.normtf, gevi.bs_normtf] = gevi_loadData('gevi_tf', 140, 7);

%% visualize
figure (1), clf
new_colors = bone(13);
%set(gca, 'colororder', new_colors), hold on
for k = 3 : 9
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 10
    plot(t, squeeze(gevi.tf(k1, :, k))', 'linewidth', 1), hold on
    end
    plot([t(end), t(end)], [-0.03, 0.03], 'k--')
end
set(gca, 'xtick', 1.4 * 3 : 1.4 : 1.4 * 9, 'XTickLabel', ones(1, 7) * 1.4), box off
plot([1.4 * 2, 1.4 * 9], [0, 0], 'k:'), ylim([-0.03, 0.03])

figure (2), clf
for k = 1 : 14
    subplot(4, 4, k)
    plot(squeeze(gevi.lintf(:, :, k))'), hold on
end



figure (3), clf
% for k = 3 : 9
%     subplot(4, 4, k)
%     plot(squeeze(gevi.normtf(:, :, k))'), hold on
% end
for k = 1 : 7
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 10
    plot(t, squeeze(gevi.normtf(k1, :, k))', 'linewidth', 1), hold on
    end
    plot([t(end), t(end)], [-0.01, 0.02], 'k--')
end
set(gca, 'xtick', 1.4 * 1 : 1.4 : 1.4 * 7, 'XTickLabel', ones(1, 7) * 1.4), box off
plot([1.4 * 0, 1.4 * 7], [0, 0], 'k:'), ylim([-0.005, 0.015])

figure (4), clf
% for k = 1 : 14
%     subplot(4, 4, k)
%     plot(squeeze(gevi.bs_normtf(:, :, k))'), hold on
% end
new_colors = flipud(pink(50));
for k = 1 : 7
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 50
    plot(t, squeeze(gevi.bs_normtf(k1, :, k))', 'color', new_colors(k1, :)), hold on
    end
    plot([t(end), t(end)], [-0.01, 0.02], 'k--')
end
set(gca, 'xtick', 1.4 * 1 : 1.4 : 1.4 * 7, 'XTickLabel', ones(1, 7) * 1.4), box off
plot([1.4 * 0, 1.4 * 7], [0, 0], 'k:'), ylim([-0.005, 0.015])

%% load gevi contrast data

gevi.ctr = []; gevi.mctr = []; gevi.linctr = []; gevi.normctr = []; gevi.bs_normctr = [];

[gevi.ctr, gevi.mctr, gevi.linctr, gevi.normctr, gevi.bs_normctr] = gevi_loadData('gevi_ctr', 140, 6);

%% visualize

figure (1), clf
%new_colors = bone(10);
%set(gca, 'colororder', new_colors), hold on
for k = 3 : 8
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 7
    plot(t, squeeze(gevi.ctr(k1, :, k))', 'linewidth', 1), hold on
    end
    plot([t(end), t(end)], [-0.05, 0.03], 'k--')
end
set(gca, 'xtick', 1.4 * 3 : 1.4 : 1.4 * 8, 'XTickLabel', ones(1, 6) * 1.4), box off
plot([1.4 * 2, 1.4 * 8], [0, 0], 'k:'), ylim([-0.02, 0.03])


figure (2), clf
for k = 1 : 12
    subplot(4, 4, k)
    plot(squeeze(gevi.linctr(:, :, k))'), hold on, ylim([-5, 15] * 10^(-3))
end

figure (3), clf
for k = 1 : 6
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 7
    plot(t, squeeze(gevi.normctr(k1, :, k))', 'linewidth', 1), hold on
    end
    plot([t(end), t(end)], [-0.01, 0.02], 'k--')
end
set(gca, 'xtick', 1.4 * 1 : 1.4 : 1.4 * 6, 'XTickLabel', ones(1, 6) * 1.4), box off
plot([1.4 * 0, 1.4 * 6], [0, 0], 'k:'), ylim([-0.005, 0.015])

figure (4), clf

new_colors = flipud(pink(50));
for k = 1 : 6
    idx = 140 * (k -1) + 1 : 140 * k;
    t = idx ./100;
    for k1 = 1 : 50
    plot(t, squeeze(gevi.bs_normctr(k1, :, k))'), hold on
    end
    plot([t(end), t(end)], [-0.01, 0.02], 'k--')
end
set(gca, 'xtick', 1.4 * 1 : 1.4 : 1.4 * 6, 'XTickLabel', ones(1, 6) * 1.4), box off
plot([1.4 * 0, 1.4 * 6], [0, 0], 'k:'), ylim([-0.005, 0.015])


%% load gcamp tf data

gcamp = [];

gcamp.tf = []; gcamp.mtf = []; gcamp.lintf = []; gcamp.normtf = []; gcamp.bs_normtf = [];

[gcamp.tf, gcamp.mtf, gcamp.lintf, gcamp.normtf, gcamp.bs_normtf] = gevi_loadData('gcamp_tf', 140, 7);

figure (1), clf
for k = 1 : 16
    subplot(4, 4, k)
    plot(squeeze(gcamp.tf(:, :, k))'), hold on
end

figure (2), clf
for k = 1 : 14
    subplot(4, 4, k)
    plot(squeeze(gcamp.lintf(:, :, k))'), hold on
end

figure (3), clf
for k = 1 : 14
    subplot(4, 4, k)
    plot(squeeze(gcamp.normtf(:, :, k))'), hold on, %ylim([-5, 15] * 10^(-3))
end

figure (4), clf
for k = 1 : 14
    subplot(4, 4, k)
    plot(squeeze(gcamp.bs_normtf(:, :, k))'), hold on, %ylim([-5, 15] * 10^(-3))
end

%% load gcamp contrast data

gcamp.ctr = []; gcamp.mctr = []; gcamp.linctr = []; gcamp.normctr = []; gcamp.bs_normctr = [];

[gcamp.ctr, gcamp.mctr, gcamp.linctr, gcamp.normctr, gcamp.bs_normctr] = gevi_loadData('gcamp_ctr', 140, 6);


figure (1), clf
for k = 1 : 14
    subplot(4, 4, k)
    plot(squeeze(gcamp.ctr(:, :, k))'), hold on, ylim([-15, 45] * 10^(-3))
end

figure (2), clf
for k = 1 : 12
    subplot(4, 4, k)
    plot(squeeze(gcamp.linctr(:, :, k))'), hold on, ylim([-15, 45] * 10^(-3))
end

figure (3), clf
for k = 1 : 12
    subplot(4, 4, k)
    plot(squeeze(gcamp.normctr(:, :, k))'), hold on, ylim([-5, 45] * 10^(-3))
end

figure (4), clf
for k = 1 : 12
    subplot(4, 4, k)
    plot(squeeze(gcamp.bs_normctr(:, :, k))'), hold on, ylim([-5, 45] * 10^(-3))
end

%% save data

save('preprocessed_data.mat', 'gevi', 'gcamp')

