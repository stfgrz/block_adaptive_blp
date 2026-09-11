function s = summarize_montecarlo(mc, early_H)
% PURPOSE
% -------
% Turn stored Monte Carlo output into the evaluation metrics of the
% study.  Every metric is defined explicitly below; horizons h = 1..H
% are evaluated (h = 0 is the shared identification step and is
% excluded from all metrics).
%
% MODEL / EQUATIONS (METRIC DEFINITIONS)
% --------------------------------------
% Let theta_hat(e,i,h,r) be estimator e's estimate of variable i's
% response at horizon h in replication r, theta0(i,h) the truth, and
% [lo, hi] the interval bounds.  With R replications:
%
%   bias(e,i,h)   = (1/R) sum_r [ theta_hat - theta0 ]
%   variance(e,i,h) = (1/R) sum_r [ theta_hat - mean_r theta_hat ]^2
%                   (Monte Carlo variance of the estimator; the 1/R
%                    convention makes the identity below exact)
%   mse(e,i,h)    = bias^2 + variance   (checked numerically below)
%   rmse(e,i,h)   = sqrt( mse ) = sqrt( (1/R) sum_r [theta_hat-theta0]^2 )
%   coverage(e,i,h) = (1/R) sum_r 1{ lo <= theta0 <= hi }
%   avg_len(e,i,h)  = (1/R) sum_r ( hi - lo )
%
% INTEGRATED RMSE -- TWO DEFINITIONS
%   irmse(e,i)    = (1/H) sum_{h=1}^{H} rmse(e,i,h)          [LEGACY]
%   irmse_h2(e,i) = (1/(H-1)) sum_{h=2}^{H} rmse(e,i,h)      [PREFERRED]
% The legacy metric is kept so every earlier number stays reproducible,
% but it is NOT a fair estimator comparison in FMAR mode with the
% default cfg.fmar.h1_mode = 'bvar': there the global baseline reports
% the Bayesian VAR at h = 1 while the block-adaptive estimators report
% an h = 1 local projection, so the h = 1 term mixes "VAR vs LP" into
% what should be "global vs adaptive".  irmse_h2 drops that horizon and
% is the metric the write-up should quote.  Setting
% cfg.fmar.h1_mode = 'lp' makes h = 1 like-for-like as well, in which
% case the two definitions answer the same question.
% Two further aggregates help separate "helps early" from "helps
% overall":
%   irmse_early(e,i) = mean of rmse over h = 2..early_H
%   irmse_late(e,i)  = mean of rmse over h = early_H+1..H
%
% Block-scale diagnostics (block-adaptive BLP only), scheme
% 'per_variable' with G = K blocks:
%
%   tau_bar(i,g,h) = (1/R) sum_r tau_mean(i,g,h,r)
%       average posterior-mean local scale, per equation i, block g,
%       horizon h.
%   tau_rep(g,r)   = mean over (i,h) of tau_mean(i,g,h,r)
%       one aggregated scale per block and replication.
%   detect_prob    = (1/R) sum_r 1{ argmax_g tau_rep(g,r) = g* }
%       probability that the TRULY misspecified block g* receives the
%       LARGEST aggregated local scale (only defined when the DGP has a
%       unique misspecified block, mc.misspec_block).
%   argmax_freq(g) = (1/R) sum_r 1{ argmax_g' tau_rep(g',r) = g }
%       full distribution of which block "wins" (reported always).
%
% PER-EQUATION detection (added after the first full experiment showed
% the aggregate metric dilutes a localised signal):
%   tau_eq(i,g,r)   = mean over h of tau_mean(i,g,h,r)      (all horizons)
%   tau_eqe(i,g,r)  = mean over h <= early_H of tau_mean    (early window)
%   detect_prob_by_eq(i)       = (1/R) sum_r 1{ argmax_g tau_eq(i,g,r)  = g* }
%   detect_prob_by_eq_early(i) = (1/R) sum_r 1{ argmax_g tau_eqe(i,g,r) = g* }
%   argmax_freq_by_eq(i,g)     = share of reps where block g wins in eq i
%       (computed on the all-horizon version; defined for every DGP, so
%        on the CORRECT DGP these same numbers are false-positive rates).
% The early window matters because transmission-channel misspecification
% concentrates at short horizons; early_H defaults to min(6, H) and can
% be passed as a second argument.  NOTE for reporting: quoting the
% detection rate of the equation known to be misspecified uses oracle
% knowledge of the location; present these as a diagnostic MAP (per
% equation x block x horizon) together with the correct-DGP rates.
%
% Optional flag diagnostic (crude detection rule):
%   flag(r) = 1{ max_g tau_rep(g,r) / median_g tau_rep(g,r) > 1.5 }
%   * On a correctly specified DGP, mean(flag) is a FALSE-POSITIVE rate.
%   * On the sparse DGP, mean(flag AND argmax = g*) is a joint
%     "flag fired AND pointed at the right block" detection rate.
% The 1.5 threshold is arbitrary and only illustrative.
%
% INPUTS
% ------
% mc : output struct of run_montecarlo.
%
% OUTPUTS
% -------
% s : struct with fields
%   .est_names, .dgp_name, .misspec_block
%   .bias, .variance, .mse, .rmse, .coverage, .avg_len : (nE x K x H)
%   .coverage_post, .avg_len_post : (nE x K x H) the same for the
%                          posterior-quantile bands of the sampled
%                          estimators (NaN where the Monte Carlo did not
%                          store them)
%   .irmse, .irmse_h2, .irmse_early, .irmse_late       : (nE x K)
%   .ibias_h2, .ivar_h2  : (nE x K) integrated |bias| and variance over
%                          h = 2..H, the bias-variance counterpart of
%                          irmse_h2
%   .tau_stats           : struct produced by tau_diagnostics for each
%                          available scale estimator ('block', 'pooled')
%   .meta                : the mc.meta reproducibility record, passed
%                          through so a summary is self-describing
%   .tau_bar                          : (K x G x H)
%   .argmax_freq                      : (G x 1)
%   .detect_prob                      : scalar or NaN
%   .false_pos_rate, .joint_detect    : scalar or NaN
%   .var_stable_share, .mean_lag1_*   : diagnostics
%
% DIMENSIONS
% ----------
% As stored by run_montecarlo (4 estimators, K variables, H+1 horizons,
% R replications, G blocks).
%
% NOTES
% -----
% Deterministic given mc.

[nE, K, Hp1, R] = size(mc.theta);
H = Hp1 - 1;
G = size(mc.tau_mean, 2);

theta0 = mc.theta_true;                 % (K x (H+1))

bias = zeros(nE, K, H);  rmse = zeros(nE, K, H);
vari = zeros(nE, K, H);  mse  = zeros(nE, K, H);
cover = zeros(nE, K, H); alen = zeros(nE, K, H);
for e = 1:nE
    for i = 1:K
        for h = 1:H
            est = squeeze(mc.theta(e, i, h + 1, :));
            err = est - theta0(i, h + 1);
            bias(e, i, h) = mean(err);
            vari(e, i, h) = mean((est - mean(est)).^2);   % 1/R convention
            mse(e, i, h)  = mean(err.^2);
            rmse(e, i, h) = sqrt(mse(e, i, h));
            l = squeeze(mc.lo(e, i, h + 1, :));
            u = squeeze(mc.hi(e, i, h + 1, :));
            cover(e, i, h) = mean(l <= theta0(i, h + 1) & ...
                                  theta0(i, h + 1) <= u);
            alen(e, i, h)  = mean(u - l);
        end
    end
end
% Secondary (posterior-quantile) bands, where the Monte Carlo stored
% them: the PRIMARY intervals are FMAR's quasi-Bayesian Newey-West
% sandwich, and the gap between the two coverages is the empirical price
% of that choice (docs/APPROXIMATIONS.md, item 3).
cover_post = nan(nE, K, H);  alen_post = nan(nE, K, H);
if isfield(mc, 'lo_post') && isfield(mc, 'hi_post')
    for e = 1:nE
        for i = 1:K
            for h = 1:H
                l = squeeze(mc.lo_post(e, i, h + 1, :));
                u = squeeze(mc.hi_post(e, i, h + 1, :));
                if all(isnan(l)), continue; end
                cover_post(e, i, h) = mean(l <= theta0(i, h + 1) & ...
                                           theta0(i, h + 1) <= u);
                alen_post(e, i, h)  = mean(u - l);
            end
        end
    end
end

% The decomposition MSE = bias^2 + variance is exact with the 1/R
% variance convention; assert it rather than trusting it.
decomp_err = max(abs(mse(:) - (bias(:).^2 + vari(:))));
assert(decomp_err < 1e-8 * max(1, max(mse(:))), ...
    'summarize_montecarlo: bias-variance decomposition off by %.3g.', decomp_err);

irmse = mean(rmse, 3);                  % (nE x K), LEGACY h = 1..H

% --- Block-scale summaries ---------------------------------------------
tau_bar = mean(mc.tau_mean, 4);         % (K x G x H)

tau_rep = zeros(G, R);                  % aggregated scale per block/rep
for r = 1:R
    tm = mc.tau_mean(:, :, :, r);       % (K x G x H)
    tau_rep(:, r) = squeeze(mean(mean(tm, 3), 1));
end
[~, winner] = max(tau_rep, [], 1);      % (1 x R)

argmax_freq = zeros(G, 1);
for g = 1:G
    argmax_freq(g) = mean(winner == g);
end

% Per-equation winners (all horizons, and early window).
if nargin < 2 || isempty(early_H)
    early_H = min(6, H);
end

% --- integrated metrics ------------------------------------------------
if H >= 2
    irmse_h2 = mean(rmse(:, :, 2:H), 3);
    ibias_h2 = mean(abs(bias(:, :, 2:H)), 3);
    ivar_h2  = mean(vari(:, :, 2:H), 3);
else
    irmse_h2 = nan(nE, K);  ibias_h2 = nan(nE, K);  ivar_h2 = nan(nE, K);
end
h_early = 2:min(early_H, H);
h_late  = min(early_H, H) + 1:H;
if isempty(h_early), irmse_early = nan(nE, K);
else,                irmse_early = mean(rmse(:, :, h_early), 3); end
if isempty(h_late),  irmse_late = nan(nE, K);
else,                irmse_late = mean(rmse(:, :, h_late), 3); end
win_eq  = zeros(K, R);   win_eqe = zeros(K, R);
for r = 1:R
    for i = 1:K
        t_all  = squeeze(mean(mc.tau_mean(i, :, :, r), 3));         % 1 x G
        t_earl = squeeze(mean(mc.tau_mean(i, :, 1:early_H, r), 3)); % 1 x G
        [~, win_eq(i, r)]  = max(t_all);
        [~, win_eqe(i, r)] = max(t_earl);
    end
end
argmax_freq_by_eq = zeros(K, G);
for i = 1:K
    for g = 1:G
        argmax_freq_by_eq(i, g) = mean(win_eq(i, :) == g);
    end
end

ratio = max(tau_rep, [], 1) ./ (median(tau_rep, 1) + eps);
flag  = ratio > 1.5;

if ~isempty(mc.misspec_block)
    gstar = mc.misspec_block;
    s.detect_prob  = mean(winner == gstar);
    s.detect_prob_by_eq       = mean(win_eq  == gstar, 2);   % (K x 1)
    s.detect_prob_by_eq_early = mean(win_eqe == gstar, 2);   % (K x 1)
    s.joint_detect = mean(flag & (winner == gstar));
    s.false_pos_rate = NaN;             % not defined under misspecification
else
    s.detect_prob  = NaN;
    s.detect_prob_by_eq       = nan(K, 1);
    s.detect_prob_by_eq_early = nan(K, 1);
    s.joint_detect = NaN;
    if strcmp(mc.dgp_name, 'correct')
        s.false_pos_rate = mean(flag);  % correct DGP: any flag is false
    else
        s.false_pos_rate = NaN;         % dense DGP: flags are not "false"
    end
end
s.flag_rate = mean(flag);               % raw flag frequency, any DGP

% --- Pack ---------------------------------------------------------------
s.est_names = mc.est_names;
s.dgp_name  = mc.dgp_name;
s.misspec_block = mc.misspec_block;
s.bias = bias;  s.rmse = rmse;  s.coverage = cover;  s.avg_len = alen;
s.variance = vari;  s.mse = mse;
s.coverage_post = cover_post;   % posterior-quantile bands (NaN where absent)
s.avg_len_post  = alen_post;
s.irmse = irmse;                 % LEGACY   (h = 1..H)
s.irmse_h2 = irmse_h2;           % PREFERRED (h = 2..H)
s.irmse_early = irmse_early;     % h = 2..early_H
s.irmse_late  = irmse_late;      % h = early_H+1..H
s.ibias_h2 = ibias_h2;  s.ivar_h2 = ivar_h2;
s.tau_bar = tau_bar;
s.argmax_freq = argmax_freq;
s.argmax_freq_by_eq = argmax_freq_by_eq;
s.early_H = early_H;
s.var_stable_share = mean(mc.diag.var_stable);
s.mean_lag1_glob   = mean(mc.diag.lag1_glob);
s.mean_lag1_block  = mean(mc.diag.lag1_block);
s.mean_n_tau_clip  = mean(mc.diag.n_tau_clip);
s.H = H;  s.K = K;  s.G = G;  s.R = R;

% --- per-estimator tau diagnostics -------------------------------------
% tau_diagnostics returns the full localisation / false-positive /
% chain-diagnostic block for one set of scales; run it for every scale
% estimator present in mc.
s.tau_stats = struct();
if isfield(mc, 'tau') && isstruct(mc.tau)
    f = fieldnames(mc.tau);
    for k = 1:numel(f)
        s.tau_stats.(f{k}) = tau_diagnostics(mc, f{k}, early_H);
    end
elseif isfield(mc, 'tau_mean')
    s.tau_stats.block = tau_diagnostics(mc, 'block', early_H);
end

if isfield(mc, 'meta'), s.meta = mc.meta; end
if isfield(mc, 'dgp_params'), s.dgp_params = mc.dgp_params; end
if isfield(mc, 'description'), s.description = mc.description; end
end