function d = tau_diagnostics(mc, which_est, early_H)
% PURPOSE
% -------
% Compact NUMERICAL diagnostics for the adaptive block scales tau, for
% ONE scale estimator ('block' = tau independent by horizon, 'pooled' =
% log tau smoothed across horizons).  These are the numbers the thesis
% needs in order to argue that the adaptive layer is a usable
% SPECIFICATION DIAGNOSTIC, separately from whether it lowers RMSE.
%
% WHAT tau MEASURES -- AND WHAT IT DOES NOT
% -----------------------------------------
% tau_{i,g,h} is the posterior scale of the deviation of block g of
% equation i's horizon-h LP coefficients FROM THE VAR-IMPLIED CENTRE.
% A large tau is evidence of PRIOR-DATA DISAGREEMENT at that
% (equation, block, horizon).  Disagreement has several possible
% sources: genuine dynamic misspecification of the fitted VAR, but also
% sampling noise in the prior centre itself, a badly chosen lambda_h,
% or simply a horizon where the LP is weakly identified.  Nothing here
% licenses reading a large tau as "the VAR is structurally wrong at
% this block"; the correct-DGP false-positive rates below exist
% precisely to calibrate that reading.
%
% METRICS
% -------
% Posterior summaries, averaged over replications:
%   .tau_bar, .tau_med_bar, .p_gt1_bar   (K x G x H)
%       MC averages of the posterior mean, posterior median and
%       posterior P(tau > 1) of each tau_{i,g,h}.
%   .tau_q_bar (K x G x H x nQ) at .tau_q_probs
%       MC average of the posterior QUANTILES of each tau_{i,g,h}.  A
%       large posterior mean can come from a shifted posterior or from a
%       heavy right tail; only the quantiles tell them apart, and the
%       half-Cauchy prior makes the second possibility a real one.
%   .tau_rep_q (G x 5)
%       quantiles ACROSS REPLICATIONS (0.05/0.25/0.5/0.75/0.95) of the
%       replication-level aggregate tau_rep(g,r) = mean over (i,h).
%
% Localisation (which block does the diagnostic point at?):
%   .argmax_freq (G x 1)          aggregate winner distribution
%   .argmax_freq_by_eq (K x G)    per-equation winner distribution
%   .argmax_freq_by_eq_early      the same on horizons h <= early_H
%   .detect_prob / _by_eq / _by_eq_early
%       share of replications whose winner is the truly misspecified
%       block g* (NaN when the DGP has no unique g*).
%   .mean_rank_true / .mean_rank_true_by_eq
%       average RANK of the true block (1 = largest tau).  A rank of
%       (G+1)/2 is what pure noise delivers, so this says how much
%       better than chance the ranking is even when it is not first.
%
% CONCENTRATION (the statistic that actually discriminates):
%   .max_argmax_freq / .max_argmax_freq_early
%       max over (equation, block) of the share of replications in which
%       that block wins in that equation, over all horizons and over the
%       early window.  Under the correct DGP the winner is uniform over
%       the G blocks in every equation, so this statistic has a
%       CALIBRATED NULL: it is the maximum of K*G frequencies each with
%       mean 1/G, and on the correct DGP it lands a little above 1/G by
%       chance alone.  On a DGP with a genuine localised conflict it goes
%       to 1 in the affected equation.  The aggregate flag rules below,
%       by contrast, average tau over equations BEFORE comparing blocks,
%       which dilutes a signal confined to one equation almost to
%       nothing -- on this project's designs they do not discriminate at
%       all, and that is reported rather than hidden.
%   .argmax_cell / .argmax_cell_early
%       which (equation, block) attains that maximum.
%
% Flag rules and FALSE POSITIVES (the same numbers on the correct DGP):
%   .flag_ratio(thr)  = share of reps with max_g tau_rep / median_g
%                       tau_rep > thr, for thr in .flag_thresholds.
%   .flag_prob(q)     = share of reps in which some (i,g) has
%                       mean_h P(tau_{i,g,h} > 1) > q, for q in
%                       .prob_thresholds.
%   .joint_detect_ratio / .joint_detect_prob
%                       flag fired AND pointed at the true block.
% On the 'correct' DGP every flag is by construction a FALSE POSITIVE,
% and .false_positive = true marks the struct so a caller cannot
% misread the same field on a misspecified DGP.
%
% Chain diagnostics (whatever the Monte Carlo stored):
%   .ess_beta_mean, .ess_logtau_mean, .lag1_beta_mean, .acc_rate_mean,
%   .n_tau_clip_mean
%
% INPUTS
% ------
% mc        : output struct of run_montecarlo.
% which_est : 'block' (default) or 'pooled'.
% early_H   : end of the early-horizon window (default min(6, H)).
%
% OUTPUTS
% -------
% d : struct documented above, plus .est, .dgp_name, .misspec_block,
%     .R, .K, .G, .H, .early_H.
%
% DIMENSIONS
% ----------
% K equations, G blocks, H horizons, R replications.
%
% NOTES
% -----
% Deterministic given mc.  Results saved before mc.tau existed are
% handled: with which_est = 'block' the legacy mc.tau_mean is used and
% the probability-based fields come back as NaN.

if nargin < 2 || isempty(which_est), which_est = 'block'; end

% --- locate the scale arrays -------------------------------------------
have_prob = false;
if isfield(mc, 'tau') && isstruct(mc.tau) && isfield(mc.tau, which_est)
    S = mc.tau.(which_est);
    tau_mean = S.mean;
    if isfield(S, 'med'), tau_med = S.med; else, tau_med = nan(size(tau_mean)); end
    if isfield(S, 'p_gt1') && ~isempty(S.p_gt1)
        p_gt1 = S.p_gt1;  have_prob = any(p_gt1(:) ~= 0);
    else
        p_gt1 = nan(size(tau_mean));
    end
elseif strcmp(which_est, 'block') && isfield(mc, 'tau_mean')
    tau_mean = mc.tau_mean;                      % legacy result file
    tau_med  = nan(size(tau_mean));
    p_gt1    = nan(size(tau_mean));
else
    error('tau_diagnostics: no scales stored for "%s".', which_est);
end

[K, G, H, R] = size(tau_mean);
if nargin < 3 || isempty(early_H), early_H = min(6, H); end
early_H = min(early_H, H);

d.est = which_est;
d.dgp_name = mc.dgp_name;
d.misspec_block = mc.misspec_block;
d.R = R;  d.K = K;  d.G = G;  d.H = H;  d.early_H = early_H;
d.false_positive = strcmp(mc.dgp_name, 'correct') || isempty(mc.misspec_block);

% --- posterior summaries averaged over replications --------------------
d.tau_bar     = mean(tau_mean, 4);
d.tau_med_bar = mean(tau_med, 4);
d.p_gt1_bar   = mean(p_gt1, 4);
d.have_prob   = have_prob;
% Posterior quantiles of each tau_{i,g,h}, averaged over replications.
% Reported alongside the mean because a large mean can come either from
% a shifted posterior or from a heavy right tail, and only the quantiles
% distinguish the two.
d.tau_q_bar = [];  d.tau_q_probs = [];
if isfield(mc, 'tau') && isstruct(mc.tau) && isfield(mc.tau, which_est)
    Sq = mc.tau.(which_est);
    if isfield(Sq, 'q') && ~isempty(Sq.q)
        d.tau_q_bar   = mean(Sq.q, 5);          % (K x G x H x nQ)
        d.tau_q_probs = Sq.probs;
    end
end

% --- replication-level aggregates --------------------------------------
tau_rep = zeros(G, R);
for r = 1:R
    tau_rep(:, r) = squeeze(mean(mean(tau_mean(:, :, :, r), 3), 1));
end
probs = [0.05 0.25 0.50 0.75 0.95];
d.tau_rep_probs = probs;
d.tau_rep_q = zeros(G, numel(probs));
for g = 1:G
    d.tau_rep_q(g, :) = empirical_quantile(tau_rep(g, :)', probs);
end

[~, winner] = max(tau_rep, [], 1);
d.argmax_freq = zeros(G, 1);
for g = 1:G, d.argmax_freq(g) = mean(winner == g); end

% --- per-equation localisation ------------------------------------------
win_eq = zeros(K, R);  win_eqe = zeros(K, R);
rank_true_eq = nan(K, R);
gstar = mc.misspec_block;
for r = 1:R
    for i = 1:K
        t_all  = squeeze(mean(tau_mean(i, :, :, r), 3));
        t_earl = squeeze(mean(tau_mean(i, :, 1:early_H, r), 3));
        [~, win_eq(i, r)]  = max(t_all);
        [~, win_eqe(i, r)] = max(t_earl);
        if ~isempty(gstar)
            [~, ord] = sort(t_all, 'descend');
            rank_true_eq(i, r) = find(ord == gstar, 1);
        end
    end
end
d.argmax_freq_by_eq       = zeros(K, G);
d.argmax_freq_by_eq_early = zeros(K, G);
for i = 1:K
    for g = 1:G
        d.argmax_freq_by_eq(i, g)       = mean(win_eq(i, :)  == g);
        d.argmax_freq_by_eq_early(i, g) = mean(win_eqe(i, :) == g);
    end
end

if ~isempty(gstar)
    d.detect_prob              = mean(winner == gstar);
    d.detect_prob_by_eq        = mean(win_eq  == gstar, 2);
    d.detect_prob_by_eq_early  = mean(win_eqe == gstar, 2);
    [~, ord_agg] = sort(tau_rep, 1, 'descend');
    rank_agg = zeros(1, R);
    for r = 1:R, rank_agg(r) = find(ord_agg(:, r) == gstar, 1); end
    d.mean_rank_true       = mean(rank_agg);
    d.mean_rank_true_by_eq = mean(rank_true_eq, 2);
    d.rank_chance          = (G + 1) / 2;
else
    d.detect_prob = NaN;
    d.detect_prob_by_eq       = nan(K, 1);
    d.detect_prob_by_eq_early = nan(K, 1);
    d.mean_rank_true = NaN;  d.mean_rank_true_by_eq = nan(K, 1);
    d.rank_chance = (G + 1) / 2;
end

% --- concentration statistic (calibrated null: 1/G per cell) -----------
[mx, lin] = max(d.argmax_freq_by_eq(:));
[i_mx, g_mx] = ind2sub([K, G], lin);
d.max_argmax_freq = mx;
d.argmax_cell = [i_mx, g_mx];
[mxe, line_] = max(d.argmax_freq_by_eq_early(:));
[i_me, g_me] = ind2sub([K, G], line_);
d.max_argmax_freq_early = mxe;
d.argmax_cell_early = [i_me, g_me];
d.argmax_chance = 1 / G;

% --- flag rules ---------------------------------------------------------
ratio = max(tau_rep, [], 1) ./ (median(tau_rep, 1) + eps);
d.flag_thresholds = [1.25 1.5 2.0];
d.flag_ratio = zeros(1, numel(d.flag_thresholds));
d.joint_detect_ratio = nan(1, numel(d.flag_thresholds));
for k = 1:numel(d.flag_thresholds)
    fl = ratio > d.flag_thresholds(k);
    d.flag_ratio(k) = mean(fl);
    if ~isempty(gstar)
        d.joint_detect_ratio(k) = mean(fl & (winner == gstar));
    end
end
d.flag_rate = d.flag_ratio(2);          % thr = 1.5, the legacy number

d.prob_thresholds = [0.50 0.75 0.90];
d.flag_prob = nan(1, numel(d.prob_thresholds));
d.joint_detect_prob = nan(1, numel(d.prob_thresholds));
if have_prob
    pmax = zeros(1, R);  pwin = zeros(1, R);
    for r = 1:R
        Pig = mean(p_gt1(:, :, 1:early_H, r), 3);    % (K x G)
        [pmax(r), lin] = max(Pig(:));
        [~, gcol] = ind2sub([K, G], lin);
        pwin(r) = gcol;
    end
    for k = 1:numel(d.prob_thresholds)
        fl = pmax > d.prob_thresholds(k);
        d.flag_prob(k) = mean(fl);
        if ~isempty(gstar)
            d.joint_detect_prob(k) = mean(fl & (pwin == gstar));
        end
    end
end

% --- chain diagnostics ---------------------------------------------------
d.ess_beta_mean = NaN;  d.ess_logtau_mean = NaN;
d.lag1_beta_mean = NaN; d.acc_rate_mean = NaN;  d.n_tau_clip_mean = NaN;
if isfield(mc, 'diag')
    D = mc.diag;
    if strcmp(which_est, 'block')
        if isfield(D, 'ess_beta_block'),   d.ess_beta_mean   = mean(D.ess_beta_block);   end
        if isfield(D, 'ess_logtau_block'), d.ess_logtau_mean = mean(D.ess_logtau_block); end
        if isfield(D, 'lag1_block'),       d.lag1_beta_mean  = mean(D.lag1_block);       end
        if isfield(D, 'n_tau_clip'),       d.n_tau_clip_mean = mean(D.n_tau_clip);       end
    else
        if isfield(D, 'pooled_ess_logtau'), d.ess_logtau_mean = mean(D.pooled_ess_logtau); end
        if isfield(D, 'pooled_acc_rate'),   d.acc_rate_mean   = mean(D.pooled_acc_rate);  end
    end
end
end
