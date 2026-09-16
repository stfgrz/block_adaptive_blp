function rel = ea_relevance_iv(Y, bvar, z, cfg, opts)
% EA_RELEVANCE_IV  Instrument relevance, timing and influence diagnostics
% for the external-instrument design (docs/DESIGN.md Sec. 5).
%
% Everything is computed on the VAR(p) innovations u_t = y_t - B' z_{t-1}
% at the BVAR posterior mean (the object the proxy identification uses),
% over the months t = p+1..T with a finite instrument value.
%
% BLOCKS OF OUTPUT (rel.<block>)
% ------------------------------
%   fs          first stage u_{s,t} = a + b z_t + e_t: b, se/t (EHW and HAC
%               with opts.hac_lag lags), F_eff = t_hac^2 (single instrument:
%               the Montiel Olea-Pflueger effective F; thresholds 12.05 /
%               15.06 / 23.11 / 37.42 for 30/20/10/5 % Nagar bias, and the
%               Staiger-Stock 10), R2 (partial R2 given the lags, since u
%               is already orthogonal to them), N, n_nonzero, sd_z, sd_us,
%               corr, and the MOP verdict string
%   naive       Delta y_{s,t} on z_t (AGKL Fig. 9 analogue), same fields
%   placebo     u_{s,t} on z_{t+k}, k in opts.leads (default -3..3): b, t_hac
%               per k.  For a correctly dated, valid instrument only k = 0
%               is significant.  A significant k = -1 (u_s(t) on z_{t-1})
%               is the signature of a summed surprise against a monthly-
%               AVERAGE indicator (part of last month's surprise sits in
%               this month's innovation; Kilian 2024); a significant LEAD
%               (k > 0) means the instrument is recorded after the indicator
%               reacts, or that lead-lag exogeneity fails (Stock & Watson
%               2018) -- both are misdating/exogeneity problems to report
%   predict     z_t on pre-event information X_t (opts.predictors, T x q,
%               e.g. ea_pre_event_info) over event months (opts.event_months
%               logical, default z ~= 0): HAC Wald F, p-value, R2, N, the
%               coefficients and t's; and the AR(1) coefficient of z (all
%               months, and event months only)
%   influence   leave-one-month-out b and F over event months; DFBETA =
%               (b - b_(-t)) / se(b); the opts.n_top most influential months
%               (ym, z, u_s, dfbeta, F without it); the first stage with the
%               crisis window (opts.crisis, default 2008m9-2009m6) excluded;
%               with z winsorised at the 1st/99th percentiles of its
%               non-zero values
%   subsample   first stage per window in opts.subsamples (default 2001-08,
%               2009-11, 2012-19)
%   anatomy     WHERE the first-stage covariance comes from: the signed
%               share of sum_t (z_t - zbar)(u_st - ubar) contributed by each
%               year (share_by_year), by each opts.subsamples window
%               (share_window), by the five largest months (share_top5), the
%               number of months that carry half of the ABSOLUTE covariance
%               mass (n_half_abs), the ten largest months (top_ym, top_share)
%               and sd(u_s) per window (sd_us_window).  An instrument whose
%               covariance sits in a handful of crisis months identifies the
%               IRF from those months only, whatever the full-sample F says
%   predict.lags   z_t on the FULL lag set the innovations condition on
%               (constant + K*p lags, event months): R2, HAC Wald/q, p-value,
%               with q = K*p restrictions -- a rough guide only -- and
%               .pred_change, the one-regressor version: z_t on the BVAR
%               one-step forecast of the indicator CHANGE (slope, HAC t, R2).
%               Both test the lead-lag exogeneity that the proxy route needs
%               CONDITIONAL on the lags: a significant slope means the raw
%               surprise is partly predictable from the VAR's own past
%
% INPUTS
% ------
%   Y (T x K), bvar (estimate_bvar_niw output; .B used), z (T x 1, NaN =
%   missing), cfg (.p, .shock_var default 1), opts: .ym (T x 1, = 12*year +
%   month; needed for the calendar-based blocks), .z_name, .hac_lag (1),
%   .leads (-3:3), .predictors, .pred_names, .event_months, .crisis ([2008
%   9; 2009 6]), .subsamples (n x 4: y0 m0 y1 m1), .n_top (5), .verbose,
%   .out_csv (long-format csv: block,name,value).
%
% OUTPUT
% ------
%   rel : struct of the blocks above plus .z_name, .N, .p.

if nargin < 5, opts = struct(); end
opts = sd(opts, 'z_name', 'z');  opts = sd(opts, 'hac_lag', 1);  opts = sd(opts, 'leads', -3:3);
opts = sd(opts, 'crisis', [2008 9; 2009 6]);  opts = sd(opts, 'n_top', 5);  opts = sd(opts, 'verbose', true);
opts = sd(opts, 'subsamples', [2001 1 2008 12; 2009 1 2011 12; 2012 1 2019 12]);
opts = sd(opts, 'out_csv', '');
[T, K] = size(Y);  p = cfg.p;
s = 1;  if isfield(cfg, 'shock_var') && ~isempty(cfg.shock_var), s = cfg.shock_var; end
z = z(:);  L = opts.hac_lag;
assert(numel(z) == T, 'ea_relevance_iv: z must have %d rows.', T);
if ~isfield(opts, 'ym') || isempty(opts.ym), opts.ym = (1:T)'; end
ym = opts.ym(:);

Zall = build_lp_regressors(Y, p);
U = Y(p + 1:T, :) - Zall(1:end - 1, :) * bvar.B;      % innovations, rows t = p+1..T
tt = (p + 1:T)';                                     % row t of Y for each innovation
us = U(:, s);  zz = z(tt);  ymm = ym(tt);
ok = isfinite(zz) & isfinite(us);
rel = struct('z_name', opts.z_name, 'p', p, 'N', sum(ok), 's', s);

% --- first stage --------------------------------------------------------------
rel.fs = fs_block(zz(ok), us(ok), L);
% --- naive (AGKL Fig. 9 analogue) ------------------------------------------------
di = [NaN; diff(Y(:, s))];
okn = isfinite(z) & isfinite(di);
rel.naive = fs_block(z(okn), di(okn), L);
% --- lead / lag placebo -------------------------------------------------------------
lead = opts.leads(:)';  pb = struct('k', lead, 'b', nan(size(lead)), 't_hac', nan(size(lead)), 'N', nan(size(lead)));
for j = 1:numel(lead)
    k = lead(j);
    idx = tt + k;  good = idx >= 1 & idx <= T;
    zk = nan(size(tt));  zk(good) = z(idx(good));
    o = isfinite(zk) & isfinite(us);
    if sum(o) > 20
        st = slope_stats(zk(o), us(o), L);
        pb.b(j) = st.b;  pb.t_hac(j) = st.t_hac;  pb.N(j) = sum(o);
    end
end
rel.placebo = pb;
% --- predictability -------------------------------------------------------------------
pr = struct('available', false);
if isfield(opts, 'predictors') && ~isempty(opts.predictors)
    X = opts.predictors;  assert(size(X, 1) == T, 'ea_relevance_iv: predictors must have T rows.');
    if isfield(opts, 'event_months') && ~isempty(opts.event_months), ev = logical(opts.event_months(:));
    else, ev = z ~= 0; end
    o = ev & isfinite(z) & all(isfinite(X), 2);
    if sum(o) > size(X, 2) + 20
        Xo = [ones(sum(o), 1), X(o, :)];  yo = z(o);
        b = Xo \ yo;  u = yo - Xo * b;
        V = hac_cov(Xo, u, L);
        q = size(X, 2);  R = [zeros(q, 1), eye(q)];
        W = (R * b)' * ((R * V * R') \ (R * b));
        pr.available = true;  pr.N = sum(o);  pr.q = q;
        pr.F_hac = W / q;  pr.p_value = 1 - chi2cdf_(W, q);
        pr.R2 = 1 - sum(u.^2) / sum((yo - mean(yo)).^2);
        pr.b = b(2:end)';  pr.t_hac = (b(2:end) ./ sqrt(diag(V(2:end, 2:end))))';
        if isfield(opts, 'pred_names'), pr.names = opts.pred_names; end
        % orthogonalised instrument: residual on event months, zero elsewhere
        zb = zeros(T, 1);  zb(o) = u;  zb(~isfinite(z)) = NaN;
        pr.z_orth = zb;
    end
end
oz = isfinite(z);
zc = z(oz);
pr.ar1_all = ar1(zc);
ze = z(oz & z ~= 0);  pr.ar1_events = ar1(ze);
% z on the FULL lag set of the VAR (rows t = p+1..T; K*p lag columns)
Xl = Zall(1:end - 1, 2:end);
evt = zz ~= 0;
if isfield(opts, 'event_months') && ~isempty(opts.event_months), evt = logical(opts.event_months(tt)); end
ol = ok & evt;
lg = struct('available', false);
if sum(ol) > size(Xl, 2) + 20
    Xo = [ones(sum(ol), 1), Xl(ol, :)];  yo = zz(ol);
    bl = Xo \ yo;  ul = yo - Xo * bl;  Vl = hac_cov(Xo, ul, L);
    ql = size(Xl, 2);  Rl = [zeros(ql, 1), eye(ql)];
    Wl = (Rl * bl)' * ((Rl * Vl * Rl') \ (Rl * bl));
    lg.available = true;  lg.N = sum(ol);  lg.q = ql;
    lg.R2 = 1 - sum(ul.^2) / sum((yo - mean(yo)).^2);
    lg.F_hac = Wl / ql;  lg.p_value = 1 - chi2cdf_(Wl, ql);
    lg.note = 'q = K*p restrictions on ~200 event months: chi-square p-value is a rough guide; see pred_change';
end
pred_d = Zall(1:end - 1, :) * bvar.B(:, s) - Y(p:T - 1, s);       % E[i_t | lags] - i_{t-1}
o2 = ok & evt;
if sum(o2) > 20
    Xp = [ones(sum(o2), 1), pred_d(o2)];  yo = zz(o2);  bp = Xp \ yo;  up = yo - Xp * bp;
    Vp = hac_cov(Xp, up, L);
    lg.pred_change = struct('b', bp(2), 't_hac', bp(2) / sqrt(Vp(2, 2)), ...
                            'R2', 1 - sum(up.^2) / sum((yo - mean(yo)).^2), 'N', sum(o2));
else
    lg.pred_change = struct('b', NaN, 't_hac', NaN, 'R2', NaN, 'N', sum(o2));
end
pr.lags = lg;
rel.predict = pr;
% --- influence -------------------------------------------------------------------------
zo = zz(ok);  uo = us(ok);  ymo = ymm(ok);  n = numel(zo);
X1 = [ones(n, 1), zo];  b_full = X1 \ uo;  se_full = rel.fs.se_ehw;
ev_idx = find(zo ~= 0);
b_loo = nan(n, 1);  F_loo = nan(n, 1);
for j = ev_idx'
    keep = true(n, 1);  keep(j) = false;
    st = slope_stats(zo(keep), uo(keep), L);
    b_loo(j) = st.b;  F_loo(j) = st.t_hac^2;
end
dfb = (b_full(2) - b_loo) / se_full;
[~, ord] = sort(abs(dfb), 'descend');  ord = ord(isfinite(dfb(ord)));
top = ord(1:min(opts.n_top, numel(ord)));
inf_ = struct('b_loo', b_loo, 'F_loo', F_loo, 'dfbeta', dfb, 'ym', ymo);
inf_.top = struct('ym', ymo(top), 'year', floor((ymo(top) - 1) / 12), 'month', ymo(top) - 12 * floor((ymo(top) - 1) / 12), ...
                  'z', zo(top), 'u_s', uo(top), 'dfbeta', dfb(top), 'F_without', F_loo(top));
cw = 12 * opts.crisis(1, 1) + opts.crisis(1, 2) : 12 * opts.crisis(2, 1) + opts.crisis(2, 2);
kc = ~ismember(ymo, cw);
inf_.crisis_window = opts.crisis;
inf_.fs_excl_crisis = fs_block(zo(kc), uo(kc), L);
inf_.n_crisis_months = sum(~kc);
nz = zo(zo ~= 0);
if numel(nz) > 20
    lo = empirical_quantile(nz, 0.01);  hi = empirical_quantile(nz, 0.99);
    zw = min(max(zo, lo), hi);
    inf_.fs_winsorised = fs_block(zw, uo, L);
    inf_.winsor_bounds = [lo hi];
end
rel.influence = inf_;
% --- subsamples ----------------------------------------------------------------------------
ss = opts.subsamples;  sub = struct('window', {}, 'fs', {});
for j = 1:size(ss, 1)
    w0 = 12 * ss(j, 1) + ss(j, 2);  w1 = 12 * ss(j, 3) + ss(j, 4);
    o = ymo >= w0 & ymo <= w1;
    if sum(o) > 24
        sub(end + 1) = struct('window', ss(j, :), 'fs', fs_block(zo(o), uo(o), L));  %#ok<AGROW>
    end
end
rel.subsample = sub;
% --- anatomy: where the first-stage covariance comes from -----------------------------
zcn = zo - mean(zo);  ucn = uo - mean(uo);
contrib = zcn .* ucn;  tot = sum(contrib);
an = struct('total', tot, 'n', n);
yrs = floor((ymo - 1) / 12);  an.year = unique(yrs)';
an.share_by_year = zeros(size(an.year));
for j = 1:numel(an.year), an.share_by_year(j) = sum(contrib(yrs == an.year(j))) / tot; end
an.window = ss;  an.share_window = nan(size(ss, 1), 1);  an.sd_us_window = nan(size(ss, 1), 1);
for j = 1:size(ss, 1)
    w0 = 12 * ss(j, 1) + ss(j, 2);  w1 = 12 * ss(j, 3) + ss(j, 4);
    o = ymo >= w0 & ymo <= w1;
    if any(o), an.share_window(j) = sum(contrib(o)) / tot;  an.sd_us_window(j) = std(uo(o)); end
end
[~, ord] = sort(abs(contrib), 'descend');
nt = min(10, n);
an.top_ym = ymo(ord(1:nt));  an.top_share = contrib(ord(1:nt)) / tot;
an.share_top5 = sum(contrib(ord(1:min(5, n)))) / tot;
cabs = cumsum(abs(contrib(ord))) / sum(abs(contrib));
an.n_half_abs = find(cabs >= 0.5, 1);
rel.anatomy = an;

% --- print -------------------------------------------------------------------------------------
if opts.verbose
    f = rel.fs;
    fprintf('ea_relevance_iv [%s], N = %d (non-zero %d), sd(z) = %.2f, sd(u_s) = %.2f\n', ...
            opts.z_name, f.N, f.n_nonzero, f.sd_z, f.sd_us);
    fprintf('  first stage  b = %8.4f  t_ehw = %5.2f  t_hac = %5.2f  F_eff = %6.2f  R2 = %.3f  [%s]\n', ...
            f.b, f.t_ehw, f.t_hac, f.F_eff, f.R2, f.mop_verdict);
    f = rel.naive;
    fprintf('  naive d i_t  b = %8.4f  t_hac = %5.2f  F = %6.2f  R2 = %.3f  (N = %d)\n', f.b, f.t_hac, f.F_eff, f.R2, f.N);
    fprintf('  placebo t_hac by lead k:');  fprintf(' k=%+d:%5.2f', [pb.k; pb.t_hac]);  fprintf('\n');
    if pr.available
        fprintf('  predictability: F_hac = %.2f (p = %.3f), R2 = %.3f, N = %d; AR(1) of z: %.2f all, %.2f events\n', ...
                pr.F_hac, pr.p_value, pr.R2, pr.N, pr.ar1_all, pr.ar1_events);
    else
        fprintf('  predictability: (no predictors given); AR(1) of z: %.2f all, %.2f events\n', pr.ar1_all, pr.ar1_events);
    end
    fprintf('  influence: top months by |DFBETA|:');
    for j = 1:numel(inf_.top.ym)
        fprintf(' %d-%02d(z=%.1f,dfb=%.2f,F-=%.1f)', inf_.top.year(j), inf_.top.month(j), inf_.top.z(j), inf_.top.dfbeta(j), inf_.top.F_without(j));
    end
    fprintf('\n  excl. %d crisis months: b = %.4f, F_eff = %.2f; winsorised: F_eff = %.2f\n', ...
            inf_.n_crisis_months, inf_.fs_excl_crisis.b, inf_.fs_excl_crisis.F_eff, ...
            tern(isfield(inf_, 'fs_winsorised'), getf(inf_, 'fs_winsorised', 'F_eff'), NaN));
    for j = 1:numel(sub)
        fprintf('  %d-%02d..%d-%02d: b = %.4f, t_hac = %.2f, F_eff = %.2f, N = %d\n', sub(j).window, ...
                sub(j).fs.b, sub(j).fs.t_hac, sub(j).fs.F_eff, sub(j).fs.N);
    end
    fprintf('  anatomy: covariance share by window');
    for j = 1:size(an.window, 1)
        fprintf(' | %d-%d: %.0f%% (sd u_s %.3f)', an.window(j, 1), an.window(j, 3), 100 * an.share_window(j), an.sd_us_window(j));
    end
    fprintf('\n           top 5 months carry %.0f%%; %d months carry half the absolute mass; largest:', 100 * an.share_top5, an.n_half_abs);
    for j = 1:min(5, nt)
        fprintf(' %d-%02d(%.0f%%)', floor((an.top_ym(j) - 1) / 12), an.top_ym(j) - 12 * floor((an.top_ym(j) - 1) / 12), 100 * an.top_share(j));
    end
    fprintf('\n');
    if lg.available
        fprintf('  predictability from the VAR lags (q = %d, N = %d): R2 = %.3f, HAC Wald/q = %.2f, p = %.4f (rough); on the BVAR-predicted change: slope %.2f, t_hac = %.2f, R2 = %.3f\n', ...
                lg.q, lg.N, lg.R2, lg.F_hac, lg.p_value, lg.pred_change.b, lg.pred_change.t_hac, lg.pred_change.R2);
    end
end
if ~isempty(opts.out_csv), write_csv(rel, opts.out_csv); end
end

% =====================================================================================================
function f = fs_block(x, y, L)
f = slope_stats(x, y, L);
f.F_eff = f.t_hac^2;  f.N = numel(x);  f.n_nonzero = sum(x ~= 0);
f.sd_z = std(x);  f.sd_us = std(y);
c = corrcoef(x, y);  f.corr = c(1, 2);
if f.F_eff > 37.42, f.mop_verdict = 'strong (bias < 5%)';
elseif f.F_eff > 23.11, f.mop_verdict = 'bias < 10%';
elseif f.F_eff > 15.06, f.mop_verdict = 'bias < 20%';
elseif f.F_eff > 12.05, f.mop_verdict = 'bias < 30%';
elseif f.F_eff > 10, f.mop_verdict = 'above 10 only';
else, f.mop_verdict = 'WEAK'; end
end

function st = slope_stats(x, y, L)
X = [ones(numel(x), 1), x(:)];
b = X \ y(:);  u = y(:) - X * b;
V_ehw = hac_cov(X, u, 0);  V_hac = hac_cov(X, u, L);
st.b = b(2);  st.se_ehw = sqrt(V_ehw(2, 2));  st.t_ehw = b(2) / st.se_ehw;
st.se_hac = sqrt(V_hac(2, 2));  st.t_hac = b(2) / st.se_hac;
st.R2 = 1 - sum(u.^2) / sum((y - mean(y)).^2);
end

function V = hac_cov(X, u, L)
XtXi = inv(X' * X);                                                    %#ok<MINV>
G = X .* u;  S = G' * G;  n = numel(u);
for l = 1:min(L, n - 1)
    Gl = G(l + 1:end, :)' * G(1:end - l, :);
    S = S + (1 - l / (L + 1)) * (Gl + Gl');
end
V = XtXi * S * XtXi;
end

function r = ar1(v)
v = v(:);  r = NaN;
if numel(v) > 10
    c = corrcoef(v(2:end), v(1:end - 1));  r = c(1, 2);
end
end

function pv = chi2cdf_(x, k)
% regularised lower incomplete gamma, base MATLAB
pv = gammainc(x / 2, k / 2);
end

function v = getf(s, f, g)
v = s.(f).(g);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = sd(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function write_csv(rel, path)
fid = fopen(path, 'w');
assert(fid > 0, 'ea_relevance_iv: cannot write %s', path);
fprintf(fid, 'instrument,block,name,value\n');
put = @(blk, nm, v) fprintf(fid, '%s,%s,%s,%.6g\n', rel.z_name, blk, nm, v);
for blk = {'fs', 'naive'}
    f = rel.(blk{1});
    for nm = {'b', 'se_ehw', 't_ehw', 'se_hac', 't_hac', 'F_eff', 'R2', 'N', 'n_nonzero', 'sd_z', 'sd_us', 'corr'}
        put(blk{1}, nm{1}, f.(nm{1}));
    end
end
for j = 1:numel(rel.placebo.k)
    put('placebo', sprintf('t_hac_lead%+d', rel.placebo.k(j)), rel.placebo.t_hac(j));
    put('placebo', sprintf('b_lead%+d', rel.placebo.k(j)), rel.placebo.b(j));
end
if rel.predict.available
    put('predict', 'F_hac', rel.predict.F_hac);  put('predict', 'p_value', rel.predict.p_value);
    put('predict', 'R2', rel.predict.R2);  put('predict', 'N', rel.predict.N);
end
put('predict', 'ar1_all', rel.predict.ar1_all);  put('predict', 'ar1_events', rel.predict.ar1_events);
put('influence', 'F_eff_excl_crisis', rel.influence.fs_excl_crisis.F_eff);
put('influence', 'b_excl_crisis', rel.influence.fs_excl_crisis.b);
if isfield(rel.influence, 'fs_winsorised'), put('influence', 'F_eff_winsorised', rel.influence.fs_winsorised.F_eff); end
for j = 1:numel(rel.influence.top.ym)
    put('influence', sprintf('top%d_ym', j), rel.influence.top.ym(j));
    put('influence', sprintf('top%d_dfbeta', j), rel.influence.top.dfbeta(j));
    put('influence', sprintf('top%d_F_without', j), rel.influence.top.F_without(j));
end
for j = 1:numel(rel.subsample)
    w = rel.subsample(j).window;
    put('subsample', sprintf('F_eff_%d%02d_%d%02d', w), rel.subsample(j).fs.F_eff);
    put('subsample', sprintf('b_%d%02d_%d%02d', w), rel.subsample(j).fs.b);
end
an = rel.anatomy;
for j = 1:size(an.window, 1)
    put('anatomy', sprintf('share_%d_%d', an.window(j, 1), an.window(j, 3)), an.share_window(j));
    put('anatomy', sprintf('sd_us_%d_%d', an.window(j, 1), an.window(j, 3)), an.sd_us_window(j));
end
for j = 1:numel(an.year), put('anatomy', sprintf('share_%d', an.year(j)), an.share_by_year(j)); end
put('anatomy', 'share_top5', an.share_top5);  put('anatomy', 'n_half_abs', an.n_half_abs);
for j = 1:numel(an.top_ym), put('anatomy', sprintf('top%d_ym', j), an.top_ym(j));  put('anatomy', sprintf('top%d_share', j), an.top_share(j)); end
if rel.predict.lags.available
    put('predict', 'lags_R2', rel.predict.lags.R2);  put('predict', 'lags_F_hac', rel.predict.lags.F_hac);
    put('predict', 'lags_p_value', rel.predict.lags.p_value);  put('predict', 'lags_q', rel.predict.lags.q);
end
put('predict', 'predchg_b', rel.predict.lags.pred_change.b);  put('predict', 'predchg_t_hac', rel.predict.lags.pred_change.t_hac);
put('predict', 'predchg_R2', rel.predict.lags.pred_change.R2);
fclose(fid);
end
