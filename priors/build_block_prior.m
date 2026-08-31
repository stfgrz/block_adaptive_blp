function bp = build_block_prior(Y, cfg)
% PURPOSE
% -------
% Build the pieces of the prior that do NOT depend on the horizon:
%   (i)  the block membership of every LP coefficient, and
%   (ii) the fixed diagonal scaling matrix D (one column per equation).
% The horizon-specific prior CENTRE comes from var_implied_lp_prior.m;
% the two are combined inside the BLP estimators.
%
% MODEL / EQUATIONS
% -----------------
% Prior on the deviation from the VAR centre, equation i, horizon h:
%
%   beta_{g,h} - mu_{g,h}^VAR ~ N(0, lambda_h^2 * tau_{g,h}^2 * D_g),
%
% where D_g is the diagonal sub-block of D for block g.  D fixes the
% UNITS of the problem so that lambda and tau are unit-free:
%
%   d_j (equation i, regressor j a lag of variable v):
%       d_j = ( std(y_i) / std(y_v) )^2,
%   d_1 (intercept): ( cfg.blp.intercept_scale * std(y_i) )^2.
%
% This is a Minnesota-style relative scaling: a prior standard
% deviation of lambda*tau in "coefficient units" corresponds to a
% deviation of roughly lambda*tau standard deviations of y_i per
% standard deviation of the regressor.  It is a pragmatic prototype
% choice, clearly separated here so it can be replaced by the exact
% scaling of the published BLP.
%
% Block scheme 'per_variable': within each equation, block g collects
% the p coefficients multiplying variable g at lags 0,...,p-1, i.e.
% z-columns 1+g, 1+K+g, ..., 1+(p-1)K+g.  The intercept has block id 0
% and is NEVER multiplied by lambda or tau (loose fixed prior).
%
% INPUTS
% ------
% Y   : (T x K) data (only sample standard deviations are used).
% cfg : configuration struct (cfg.p, cfg.blocks.scheme,
%       cfg.blp.intercept_scale).
%
% OUTPUTS
% -------
% bp : struct with fields
%   .block_id (m x 1) integer block labels; 0 = unshrunk intercept
%   .G        scalar  number of tau-carrying blocks (= K here)
%   .d_all    (m x K) column i = diagonal of D for equation i
%   .p_g      (G x 1) number of coefficients per block (= p each)
%
% DIMENSIONS
% ----------
% m = 1 + K*p.
%
% NOTES
% -----
% Only marginal standard deviations of Y are used, no estimated
% dynamics, so D is "fixed" relative to the Bayesian machinery.

[~, K] = size(Y);
p = cfg.p;
m = 1 + K * p;

assert(strcmp(cfg.blocks.scheme, 'per_variable'), ...
    'build_block_prior: only scheme ''per_variable'' is implemented.');

% --- Block membership -------------------------------------------------
block_id = zeros(m, 1);                 % intercept keeps id 0
for lag = 0:p-1
    for v = 1:K
        block_id(1 + lag*K + v) = v;    % lag of variable v -> block v
    end
end
G   = K;
p_g = zeros(G, 1);
for g = 1:G
    p_g(g) = sum(block_id == g);
end
assert(all(p_g == p), 'build_block_prior: unexpected block sizes.');

% --- Diagonal scaling D per equation ----------------------------------
sd_y = std(Y, 0, 1)';                   % (K x 1) sample std devs
assert(all(sd_y > 0), 'build_block_prior: a variable has zero variance.');

d_all = zeros(m, K);
for i = 1:K
    d = zeros(m, 1);
    d(1) = (cfg.blp.intercept_scale * sd_y(i))^2;   % loose intercept
    for j = 2:m
        v = block_id(j);
        d(j) = (sd_y(i) / sd_y(v))^2;
    end
    d_all(:, i) = d;
end

bp.block_id = block_id;
bp.G        = G;
bp.d_all    = d_all;
bp.p_g      = p_g;
end
