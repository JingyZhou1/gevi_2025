%   % a1_load data

function [data, mdata, lin_data, norm_data, bsdata] = gevi_loadData(condNm, dtLength, nConds)
%
% INPUTS ------------------------------------------
% condNm: which datasets to load. Options are "gevi_tf", "gevi_ctr",
%         "gcamp_tf", "gcamp_ctr", "vsd_tf".
% dtLength: length of a trial. A typical number would be 140. 
% nConds: total number of conditions. 

% OUTPUTS -----------------------------------------
% data: all data
% mdata: the averaged dataset

visualize = 1;

%% dataset names

switch condNm
    case 'gevi_tf'

        fileNms = {'M29D20230503R1ROITC.mat', 'M28D20230508R0ROITC.mat', 'M28D20230426R1ROITC.mat', 'M28D20230502R0ROITC.mat',  ...
    'M29D20230510R0ROITC.mat', 'M29D20230510R1ROITC.mat', 'M29D20230510R2ROITC.mat', 'M30D20240119R0ROITC.mat', 'M30D20240126R2ROITC.mat', 'M30D20240201R0ROITC.mat'};
        
        dtType = 'gevi';
        tpLoc  = 'tf_alldata';
        stim   = mkStim('TF', 200);
        f_lth = 45; % fast filter length

    case 'gevi_ctr'
        
        fileNms = {'M28D20230424R0ROITC.mat', 'M28D20230426R0ROITC.mat', 'M28D20230502R1ROITC.mat', 'M28D20230505R1ROITC.mat', 'M30D20240119R1ROITC.mat', ...
    'M30D20240124R2ROITC.mat', 'M30D20240126R0ROITC.mat'};

        % excluded 'M30D20240126R1ROITC.mat' for some large artifact that I
        % cannot get rid of through pre-processing model fitting

        dtType = 'gevi';
        tpLoc  = 'ctr_alldata';
        stim   = mkStim('contrast', 200);
        f_lth = 45; % fast filter length

    case 'gcamp_tf'
    
        fileNms = {'M30D20240119R3ROITC.mat', 'M30D20240124R0ROITC.mat', 'M30D20240209R1ROITC.mat', 'M31D20230509R2ROITC.mat', 'M31D20230509R3ROITC.mat', ...
    'M31D20230516R1ROITC.mat', 'M31D20230516R2ROITC.mat', 'M31D20230516R3ROITC.mat', 'M31D20230523R2ROITC.mat', 'M31D20230523R3ROITC.mat', 'M31D20230523R4ROITC.mat'};

        dtType = 'gcamp';
        tpLoc = 'tf_alldata';
        stim   = mkStim('TF', 200);
        f_lth = 50; % fast filter length

    case 'gcamp_ctr'

        fileNms = {'M19D20180301R1ROITCExcNuc.mat', 'M19D20180312R0ROITCExcNuc.mat', 'M19D20180320R0ROITCExcNuc.mat', ...
    'M19D20180320R1ROITCExcNuc.mat', 'M19D20180405R2ROITCExcNuc.mat', 'M30D20240119R2ROITC.mat', 'M30D20240124R1ROITC.mat', 'M30D20240209R0ROITC.mat'};

        % excluded this dataset: 'M19D20180209R3ROITCExcNuc.mat',

        dtType = 'gcamp';
        tpLoc  = 'ctr_alldata';
        stim   = mkStim('contrast', 200);
        f_lth  = 50; % fast filter length

    case 'vsd_tf'
       
        fileNms = {'M28D20230526R6ROITC.mat', 'M29D20230505R2ROITC.mat', 'M29D20230505R3ROITC.mat', 'M29D20230511R0ROITC.mat', 'M29D20230511R1ROITC.mat', 'M29D20230511R2ROITC.mat'};

        dtType = 'vsd';
        tpLoc  = 'tf_alldata';
        stim   = mkStim('TF', 200); %??????? double check
        f_lth  = 45; % fast filter length
      

    otherwise

        error('condNm does not exist.')

end

%% prepping to fit data

% linear filter length
nFast = 6;
nSlow = 6;


%% load data

data = []; mdata = []; lin_w = []; lin_data = []; norm_data = [];

nFiles = length(fileNms);

totalNconds = nConds * 2 + 2; % assuming that we have 2 repeats in each dataset, and 2 blank trials. DOUBLE CHECK!!

dtLoc = '/Users/jyzhou/Library/CloudStorage/GoogleDrive-jyz205@nyu.edu/My Drive/gd_Projects/GEVI/';


for k = 1 : nFiles
    a = load(fullfile(dtLoc, tpLoc, dtType, fileNms{k}));
    % deal with the gcamp contrast conditions (different data order) -----
    if strcmp(condNm, 'gcamp_ctr') & (k < 6)
        ctr_idx1 = [1, 2, 13, 11, 9, 7, 5, 3, 14, 12, 10, 8, 6, 4];
        a.TC = a.TC(:, ctr_idx1);
    end

    this_dtLength = size(a.TC, 1)

    % cut stimulus --------------------------------------------
    stim1 = stim(1 : this_dtLength, :);

    % make fast filter ----------------------------------------
    t_trial = [1 : this_dtLength] * 0.01;
    fBasis  = mkBasis(t_trial(1 : f_lth), nFast, 'fast_longdt');
   

    % make slow filter ----------------------------------------
    if this_dtLength < 140,
        sBasis  = mkBasis(t_trial, nSlow, 'slow_gevi');
    elseif this_dtLength == 140,
        sBasis  = mkBasis(t_trial, nSlow, 'slow');
    else
        sBasis  = mkBasis(t_trial, nSlow, 'slow_longdt');
    end
    sBasis = [zeros(size(sBasis, 1), 1), sBasis(:, 1 : end - 1)];

    % make the initial point 0
    for k1 = 1 : size(a.TC, 2)
        a.TC(:, k1) = a.TC(:, k1) - a.TC(1, k1); 
    end

    % concatenate data and stimulus ---------------------------
    dt_concat   = reshape(a.TC(:, 3 : end), [], 1);
    stim_concat = reshape(stim1(:, 3 : end), [], 1);

    % make basis function -------------------------------------
    basis = concatenateBasisAcrossConditions(fBasis, sBasis, stim_concat, t_trial);


    % fit linear model ----------------------------------------
    % ---------------------------------------------------------
    lin_w(k, :) = least_square(basis', dt_concat);

    % get de-trended data:
    fast_left = dt_concat' - lin_w(k, 1 : nConds * 2 * nSlow) * basis(1 : nConds * 2 * nSlow, :);

    % reshape the de-trended data:
    lin_fastleft{k} = reshape(fast_left, [], nConds * 2);


    % fit normalization model ---------------------------------
    % ---------------------------------------------------------
    init = [lin_w(k, nConds * 2 * nSlow + 1 : end), 0.03, 4, .1, .6, 0.05, 5];

    for iter = 1 : 10
        if iter == 1
            norm_prm = init;
        end
            norm_prm = fminsearch(@(x) fitModel(t_trial, x, stim_concat',  basis(1 : nConds * 2 * nSlow, :), ...
                basis(nConds * 2 * nSlow + 1 : end, :), fast_left), norm_prm);
        % fit slow component
        fastprd   = real(normModel2(t_trial, norm_prm, stim_concat', basis(1 : nConds * 2 * nSlow , :), basis(nConds * 2 * nSlow  + 1 : end, :)));
        slow_w    = least_square(basis(1 : nConds * 2 * nSlow, :)', (dt_concat' - fastprd)');
        slowprd   = slow_w' * basis(1 : nConds * 2 * nSlow, :);
        fast_left = dt_concat' - slowprd;
    end
     norm_fastleft{k} = reshape(fast_left, [], nConds * 2);

     figure (200), clf
     plot(dt_concat'), hold on
     plot(slowprd + fastprd)


    % cut or pad data (the original data) ---------------------
    if size(a.TC, 1) < dtLength % when the data length is shorter, we pad the data
        pad_lth = dtLength - size(a.TC, 1);
        data(k, :, :)      = cat(1, a.TC, nan(pad_lth, totalNconds));
        lin_data(k, :, :)  = cat(1, lin_fastleft{k}, nan(pad_lth, nConds * 2));
        norm_data(k, :, :) = cat(1, norm_fastleft{k}, nan(pad_lth, nConds * 2));
    else % when the data length is longer than the desirable length, we cut the data
        data(k, :, :)      = a.TC(1 : dtLength, :);
        lin_data(k, :, :)  = lin_fastleft{k}(1 : dtLength, :);
        norm_data(k, :, :) = norm_fastleft{k}(1 : dtLength, :);
    end
end

%% fix a couple inconsistencies in contrast data

% if strcmp(condNm, 'gcamp_ctr')
%     ctr_idx1 = [1, 2, 13, 11, 9, 7, 5, 3, 14, 12, 10, 8, 6, 4];
% 
%     for k = 1  : 5
%         data(k, :, :) = data(k, :, ctr_idx1);
%     end
% end

%% average across data

mdata = (data(:, :, 3 : 3 + nConds - 1) + data(:, :, 3 + nConds : 2 + nConds * 2))/2;

%% concatenate data

mdata = reshape(mdata, nFiles, dtLength * nConds);

%% create the bootstrapped mean data

% create 50 sets of mean data

nSamples1 = 50; % how many 
nSamples2 = 50; % how many synthesized subjects to average across

bsdata = [];

for iSample = 1 : nSamples1
    for iCond = 1 : nConds * 2
        idx = randi([1, nFiles], 1, nSamples2);
        bsdata(iSample, :, iCond) = nanmean(norm_data(idx, :, iCond));
    end
end


%% visualize data

if visualize

    figure (100), clf
    for k1 = 1 : nConds * 2 + 2
        for k = 1 : nFiles
            subplot(4, 5, k1),
            plot(squeeze(data(k, :, k1))), hold on
        end
        set(gca, 'ylim', [-0.05, 0.05]), axis square
    end
%       figure (100), clf
%       for k = 1 : nConds
%           idx = 140 * (k -1) + 1 : 140 * k;
%           t = idx ./100;
%     
%           plot(t, squeeze(data(k, :, :))), hold on
%       end

end

end