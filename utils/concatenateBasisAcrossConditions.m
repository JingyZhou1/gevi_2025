function finalBasis = concatenateBasisAcrossConditions(fBasis, sBasis, stim, t)

% fBasis: fast basis, n basis x time course
% sBasis: slow basis, n basis x time course
% stim  : all stimulus conditions concatenated, a vector
% t     : time course for a single trial

%% set up

nFast = min(size(fBasis));
nSlow = min(size(sBasis));

nt    = length(t);
nCond = length(stim)/nt; % assuming each trial share the same time course

%% concatenate slow basis

sBasis_allConditions = zeros(nSlow * nCond, nt * nCond); % number of basis x basis time courses

for istim = 1 : nCond
    range1 = nSlow * (istim - 1) + 1 : nSlow * istim;
    range2 = nt * (istim - 1) + 1 : nt * istim;
    sBasis_allConditions(range1, range2) = sBasis;
end

%% convolve fast basis with the stimulus vector

fBasis_allConditions = zeros(nFast, size(sBasis_allConditions, 2));

for ibasis = 1 : nFast
    fBasis_allConditions(ibasis, :) = convCut(stim, fBasis(ibasis, :));
end

%% combine fast and slow bases

finalBasis = [sBasis_allConditions; fBasis_allConditions];

end