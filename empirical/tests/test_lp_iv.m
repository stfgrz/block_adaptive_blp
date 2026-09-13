function test_lp_iv()
% TEST_LP_IV  Checks for the frequentist LP-IV benchmark (estimate_lp_iv).
%
% PURPOSE
% -------
% (a) Consistency: on a stationary VAR(2), K = 3, T = 6000, with an
%     instrument correlated only with structural shock 1, LP-IV recovers
%     the true IRFs within 3 HAC se at every (i, h) and within 0.05 in
%     absolute value at h <= 3.
% (b) Estimand equivalence with the route used by the Bayesian estimators:
%     Gamma_h * b_z with b_z = cov(u, z) / cov(u_1, z) from OLS VAR(2)
%     innovations and Gamma_h the y_t block of the OLS LP on
%     [1, y_t, ..., y_{t-p+1}] agrees with theta_lpiv within 0.03.
% (c) The AR set contains the 2SLS point; with a WEAK instrument it is
%     wider than the Wald band or flagged unbounded in most cells.
% (d) NaN rows in z are dropped without error and N_h shrinks accordingly.
% (e) theta_s(0) == 1 exactly.
%
% MODEL / EQUATIONS
% -----------------
%   y_t = c + A1 y_{t-1} + A2 y_{t-2} + B0 eta_t,  eta_t ~ N(0, I),
%   z_t = eta_{1,t} + 0.7 xi_t  (strong) or 0.05 eta_{1,t} + xi_t (weak),
%   true theta(:, h+1) = Psi_h B0(:,1) / B0(1,1),  Psi_h = J F^h J'.
%
% NOTES
% -----
% Seeded; runtime well under a minute.

fprintf('--- test_lp_iv ------------------------------------------\n');
rng(20260913, 'twister');

% ---- true VAR(2) (the repo's DGP-1 parameters) ------------------------
A1 = [0.55 0.10 0.05; 0.10 0.50 0.10; 0.05 0.10 0.45];
A2 = [0.15 0.00 0.00; 0.00 0.10 0.05; 0.00 0.05 0.15];
c  = [0.20; 0.10; 0.10];
B0 = [1.00 0.00 0.00; 0.30 0.90 0.00; 0.15 0.25 0.80];
K = 3;  p = 2;  H = 6;
F = [A1 A2; eye(K) zeros(K)];
assert(max(abs(eig(F))) < 0.95, 'test_lp_iv: true VAR not stable.');
theta_true = zeros(K, H + 1);
Fh = eye(2 * K);
for h = 0:H
    Psi = Fh(1:K, 1:K);
    theta_true(:, h + 1) = Psi * B0(:, 1) / B0(1, 1);
    Fh = Fh * F;
end

cfg = struct('p', p, 'H', H, 'ci_level', 0.90, 'shock_var', 1);
zc  = normal_quantile(0.95);

% ---- (a) strong instrument, T = 6000 ----------------------------------
T = 6000;
[Y, eta] = sim_var2(A1, A2, c, B0, T, 500);
z = eta(1, :)' + 0.7 * randn(T, 1);
tic;
lp = estimate_lp_iv(Y, z, cfg);
fprintf('    estimate_lp_iv (T = %d, H = %d): %.2f s\n', T, H, toc);

err  = abs(lp.theta - theta_true);
tsd  = err ./ max(lp.se, eps);  tsd(1, 1) = 0;
assert(all(tsd(:) <= 3), 'test_lp_iv (a): |theta - true| > 3 se, max = %.2f', max(tsd(:)));
assert(all(all(err(:, 1:4) < 0.05)), ...
    'test_lp_iv (a): |theta - true| >= 0.05 at h <= 3, max = %.3f', max(max(err(:, 1:4))));
assert(all(isfinite(lp.se(:))) && all(lp.se(:) >= 0));
assert(all(lp.lo(:) <= lp.theta(:)) && all(lp.hi(:) >= lp.theta(:)));
fprintf('    (a) max |err| h<=3 = %.4f, max |err|/se = %.2f, first-stage F_eff = %.1f\n', ...
    max(max(err(:, 1:4))), max(tsd(:)), lp.first_stage.F_eff);
assert(lp.first_stage.F_eff > 37.42, 'test_lp_iv (a): strong instrument not strong?');
% lag-augmented variant agrees too
err_la = abs(lp.theta_la - theta_true);
tsd_la = err_la ./ max(lp.se_la, eps);  tsd_la(1, 1) = 0;
assert(all(tsd_la(:) <= 3.5), 'test_lp_iv (a): LA variant off by > 3.5 se.');
assert(all(lp.N_h == T - p - (0:H)), 'test_lp_iv (a): N_h wrong without NaNs.');

% ---- (e) unit effect ---------------------------------------------------
assert(lp.theta(1, 1) == 1 && lp.se(1, 1) == 0, 'test_lp_iv (e): theta_s(0) ~= 1.');
assert(lp.theta_la(1, 1) == 1, 'test_lp_iv (e): LA theta_s(0) ~= 1.');
assert(lp.ar_lo(1, 1) == 1 && lp.ar_hi(1, 1) == 1);

% ---- (b) Gamma_h * b_z route (shared estimand) --------------------------
Zall = build_lp_regressors(Y, p);                 % [1, y_t, ..., y_{t-p+1}]
Z1 = Zall(1:end - 1, :);  Y1 = Y(p + 1:T, :);
Bv = Z1 \ Y1;  U = Y1 - Z1 * Bv;                  % OLS VAR(2) innovations
zv = z(p + 1:T);                                  % z_t aligned with u_t (t = p+1..T)
cuz = (U - mean(U))' * (zv - mean(zv)) / (numel(zv) - 1);
b_z = cuz / cuz(1);
idx = 2:K + 1;
dev = zeros(1, H);
for h = 1:H
    Zh = Zall(1:end - h, :);  Yh = Y(p + h:T, :);
    Bh = Zh \ Yh;
    Gam = Bh(idx, :)';                            % K x K block on y_t
    dev(h) = max(abs(Gam * b_z - lp.theta(:, h + 1)));
end
fprintf('    (b) max |Gamma_h b_z - theta_lpiv| over h = 1..%d: %.4f\n', H, max(dev));
assert(max(dev) < 0.03, 'test_lp_iv (b): Gamma_h b_z and LP-IV disagree (%.3f).', max(dev));
assert(max(abs(b_z - lp.theta(:, 1))) < 0.03, 'test_lp_iv (b): impact vectors disagree.');

% ---- (c) AR set: contains the point; weak instrument -> wide/unbounded --
inside = lp.ar_lo <= lp.theta + 1e-12 & lp.ar_hi >= lp.theta - 1e-12;
assert(all(inside(~lp.ar_empty)), 'test_lp_iv (c): AR set excludes the 2SLS point.');
assert(~any(lp.ar_empty(:)), 'test_lp_iv (c): AR set empty with a strong instrument.');
% strong case: AR and Wald should be of similar width away from h = 0
wald_w = 2 * zc * lp.se;  ar_w = lp.ar_hi - lp.ar_lo;
ratio = ar_w(:, 2:end) ./ wald_w(:, 2:end);
assert(all(ratio(:) > 0.6 & ratio(:) < 1.8), ...
    'test_lp_iv (c): strong-IV AR/Wald width ratio outside [0.6, 1.8].');

Tw = 300;
[Yw, etaw] = sim_var2(A1, A2, c, B0, Tw, 500);
zw = 0.05 * etaw(1, :)' + randn(Tw, 1);
lpw = estimate_lp_iv(Yw, zw, cfg);
fprintf('    (c) weak instrument: first-stage F_eff = %.2f\n', lpw.first_stage.F_eff);
inside = lpw.ar_lo <= lpw.theta + 1e-12 & lpw.ar_hi >= lpw.theta - 1e-12;
assert(all(inside(~lpw.ar_empty)), 'test_lp_iv (c): weak-IV AR set excludes the point.');
wider = (lpw.ar_hi - lpw.ar_lo) > 2 * zc * lpw.se + 1e-12;
flag  = wider | lpw.ar_unbounded;
cells = true(size(flag));  cells(1, 1) = false;   % drop the identity cell
frac  = mean(flag(cells));
fprintf('    (c) weak instrument: %.0f%% of cells AR wider than Wald or unbounded\n', 100 * frac);
assert(frac >= 0.5, 'test_lp_iv (c): weak-IV AR sets not wider/unbounded enough.');

% ---- (d) NaN rows in z are dropped ------------------------------------
zn = z;  drop = 1000:1010;  zn(drop) = NaN;
lpn = estimate_lp_iv(Y, zn, cfg);
assert(all(lpn.N_h == lp.N_h - numel(drop)), 'test_lp_iv (d): N_h did not drop by %d.', numel(drop));
assert(all(isfinite(lpn.theta(:))) && all(isfinite(lpn.se(:))));
assert(max(abs(lpn.theta(:) - lp.theta(:))) < 0.02, 'test_lp_iv (d): NaN drop moved theta a lot.');
assert(lpn.first_stage.N == lp.first_stage.N - numel(drop));
% a NaN in the initial-lag rows (never a regression row) changes nothing
zn2 = z;  zn2(1:p) = NaN;
lpn2 = estimate_lp_iv(Y, zn2, cfg);
assert(all(lpn2.N_h == lp.N_h) && max(abs(lpn2.theta(:) - lp.theta(:))) == 0);

% ---- optional flag off -------------------------------------------------
lp0 = estimate_lp_iv(Y, z, cfg, struct('lag_augment', false, 'ar_n_grid', 100));
assert(isempty(lp0.theta_la) && lp0.opts.ar_n_grid == 101);
assert(max(abs(lp0.theta(:) - lp.theta(:))) == 0);

fprintf('PASS  test_lp_iv\n\n');
end

% =====================================================================
function [Y, eta] = sim_var2(A1, A2, c, B0, T, burn)
K = size(A1, 1);
n = T + burn;
eta_all = randn(K, n);
Yall = zeros(n, K);
for t = 3:n
    Yall(t, :) = (c + A1 * Yall(t - 1, :)' + A2 * Yall(t - 2, :)' + B0 * eta_all(:, t))';
end
Y   = Yall(burn + 1:end, :);
eta = eta_all(:, burn + 1:end);
end
