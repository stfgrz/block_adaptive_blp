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
%   rmse(e,i,h)   = sqrt( (1/R) sum_r [ theta_hat - theta0 ]^2 )
%   coverage(e,i,h) = (1/R) sum_r 1{ lo <= theta0 <= hi }
%   avg_len(e,i,h)  = (1/R) sum_r ( hi - lo )
%   irmse(e,i)    = (1/H) sum_{h=1}^{H} rmse(e,i,h)
%                   ("integrated RMSE": average RMSE across horizons)
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
%   .bias, .rmse, .coverage, .avg_len : (4 x K x H)
%   .irmse                            : (4 x K)
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
cover = zeros(nE, K, H); alen = zeros(nE, K, H);
for e = 1:nE
    for i = 1:K
        for h = 1:H
            err = squeeze(mc.theta(e, i, h + 1, :)) - theta0(i, h + 1);
            bias(e, i, h) = mean(err);
            rmse(e, i, h) = sqrt(mean(err.^2));
            l = squeeze(mc.lo(e, i, h + 1, :));
            u = squeeze(mc.hi(e, i, h + 1, :));
            cover(e, i, h) = mean(l <= theta0(i, h + 1) & ...
                                  theta0(i, h + 1) <= u);
            alen(e, i, h)  = mean(u - l);
        end
    end
end
irmse = mean(rmse, 3);                  % (nE x K)

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
s.irmse = irmse;
s.tau_bar = tau_bar;
s.argmax_freq = argmax_freq;
s.argmax_freq_by_eq = argmax_freq_by_eq;
s.early_H = early_H;
s.var_stable_share = mean(mc.diag.var_stable);
s.mean_lag1_glob   = mean(mc.diag.lag1_glob);
s.mean_lag1_block  = mean(mc.diag.lag1_block);
s.mean_n_tau_clip  = mean(mc.diag.n_tau_clip);
s.H = H;  s.K = K;  s.G = G;  s.R = R;
end