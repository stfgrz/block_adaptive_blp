function chk = ea_check_series(name, ym, val, need_ym0, need_ym1)
% EA_CHECK_SERIES  Validate one monthly outcome series before it is used.
%
% PURPOSE
% -------
% The empirical chapter is only as good as its inputs, and the two ways
% an input can be silently wrong are (i) it does not actually cover the
% estimation window, and (ii) it is not real data at all -- a placeholder
% or synthetic series left over from wiring the pipeline up.  Both fail
% quietly: the estimators happily produce a full set of IRFs, tau maps
% and CSVs from fabricated numbers, and nothing in the output says so.
% This function makes both failures loud.
%
% CHECKS
% ------
% A. COVERAGE.  The series must span [need_ym0, need_ym1] with no missing
%    month inside it (a hole would surface much later as a NaN error in
%    assemble_dataset, without saying which series caused it).
% B. NON-DEGENERACY.  It must actually vary.
% C. SYNTHETIC-DATA DETECTION.  Three fingerprints, any of which is
%    essentially impossible for a genuine monthly macro/financial series
%    and typical of a generated one:
%
%      roughness  r = sd(diff2 y) / sd(diff y).
%          For any real monthly series the month-to-month change is close
%          to unforecastable, so diff2 has a LARGER standard deviation
%          than diff and r is around 1.4 (exactly sqrt(2) for a random
%          walk).  A smooth generated path (sine, polynomial, slow AR
%          with no innovations) has r far below 1.  Measured on the
%          placeholder files this package originally shipped: 0.00 (a
%          perfectly linear HICP ramp), 0.09, 0.17, 0.17.  The threshold
%          is 0.50, i.e. roughly three times the largest placeholder
%          value and less than half the value real data produces.
%
%      modal-step share = fraction of first differences exactly equal to
%          the single most common first difference.  Real data are
%          continuous-valued, so this is ~0.  A constant-increment ramp
%          gives 1.0; a series clipped at a floor gives a large share of
%          exact zeros.  Threshold 0.20.
%
%      longest flat run = longest run of exactly-repeated LEVELS.  A
%          monthly average of a market rate never repeats to four
%          decimals for years.  Threshold 12 months.
%
% INPUTS
% ------
%   name     : series label used in messages (e.g. 'rate1y_ea.csv').
%   ym       : (n x 1) month index, 12*year + month, ascending.
%   val      : (n x 1) observations.
%   need_ym0 : first month the pipeline requires (12*year + month).
%   need_ym1 : last month the pipeline requires.  Pass [] for both to
%              skip the coverage check (used when only the shape of the
%              series is being checked).
%
% OUTPUT
% ------
%   chk.ok        : true if every check passed.
%   chk.fatal     : true if a check failed that makes the series unusable
%                   (coverage, degeneracy, or a synthetic fingerprint).
%   chk.msgs      : cellstr of human-readable problem descriptions.
%   chk.roughness : the r statistic above (NaN if n < 3).
%   chk.modal_share, chk.max_flat_run, chk.n_obs, chk.ym0, chk.ym1.
%
% NOTES
% -----
% * Base MATLAB / Octave only.
% * The thresholds are deliberately far from both regimes; this is a
%   tripwire for fabricated inputs, not a statistical test.  If a genuine
%   series ever trips it (a policy rate literally fixed for over a year
%   is the plausible case), the message says which fingerprint fired so
%   the check can be overridden knowingly -- see fetch_outcome_data's
%   opts.skip_plausibility.

chk = struct('ok', true, 'fatal', false, 'msgs', {{}}, ...
             'roughness', NaN, 'modal_share', NaN, 'max_flat_run', NaN, ...
             'n_obs', numel(ym), 'ym0', NaN, 'ym1', NaN);
if isempty(ym)
    chk.ok = false;  chk.fatal = true;
    chk.msgs{end + 1} = sprintf('%s: no observations parsed.', name);
    return
end
chk.ym0 = ym(1);  chk.ym1 = ym(end);

% --- A. coverage ----------------------------------------------------------
if nargin >= 5 && ~isempty(need_ym0) && ~isempty(need_ym1)
    if ym(1) > need_ym0 || ym(end) < need_ym1
        chk.ok = false;  chk.fatal = true;
        chk.msgs{end + 1} = sprintf( ...
            '%s: covers %s..%s but the sample needs %s..%s.', name, ...
            ym_str(ym(1)), ym_str(ym(end)), ym_str(need_ym0), ym_str(need_ym1));
    end
    inwin = ym >= need_ym0 & ym <= need_ym1;
    n_need = need_ym1 - need_ym0 + 1;
    if sum(inwin) < n_need
        chk.ok = false;  chk.fatal = true;
        chk.msgs{end + 1} = sprintf( ...
            '%s: %d of %d months missing inside %s..%s.', name, ...
            n_need - sum(inwin), n_need, ym_str(need_ym0), ym_str(need_ym1));
    end
end
if any(diff(ym) == 0)
    chk.ok = false;  chk.fatal = true;
    chk.msgs{end + 1} = sprintf( ...
        ['%s: duplicated TIME_PERIOD values -- the file probably holds ' ...
         'more than one series (several geo/unit dimensions). Re-export ' ...
         'a single series.'], name);
end

% --- B/C. shape statistics, computed on the required window when given ---
v = val(:);
if nargin >= 5 && ~isempty(need_ym0) && ~isempty(need_ym1)
    sel = ym >= need_ym0 & ym <= need_ym1;
    if sum(sel) >= 3, v = val(sel); end
end
n = numel(v);

if n >= 2
    d1 = diff(v);
    if std_(d1) <= 0
        chk.ok = false;  chk.fatal = true;
        chk.msgs{end + 1} = sprintf('%s: series is constant or a perfect ramp.', name);
    end
    chk.modal_share = modal_share(d1);
    chk.max_flat_run = max_run(v);
end
if n >= 3
    d1 = diff(v);  d2 = diff(d1);
    if std_(d1) > 0
        chk.roughness = std_(d2) / std_(d1);
    end
end

% --- C. synthetic fingerprints -------------------------------------------
if ~isnan(chk.roughness) && chk.roughness < 0.50
    chk.ok = false;  chk.fatal = true;
    chk.msgs{end + 1} = sprintf( ...
        ['%s: LOOKS SYNTHETIC -- roughness sd(diff2)/sd(diff) = %.3f, far ' ...
         'below the ~1.4 of any real monthly series. This is the signature ' ...
         'of a generated smooth path, not observed data.'], name, chk.roughness);
end
if ~isnan(chk.modal_share) && chk.modal_share > 0.20
    chk.ok = false;  chk.fatal = true;
    chk.msgs{end + 1} = sprintf( ...
        ['%s: LOOKS SYNTHETIC -- %.0f%% of month-to-month changes are ' ...
         'exactly identical (a constant-increment ramp or a clipped ' ...
         'series).'], name, 100 * chk.modal_share);
end
if ~isnan(chk.max_flat_run) && chk.max_flat_run > 12
    chk.ok = false;  chk.fatal = true;
    chk.msgs{end + 1} = sprintf( ...
        ['%s: LOOKS SYNTHETIC -- the level is exactly unchanged for %d ' ...
         'consecutive months.'], name, chk.max_flat_run);
end
end

% -------------------------------------------------------------------------
function s = ym_str(v)
s = sprintf('%d-%02d', floor((v - 1) / 12), v - 12 * floor((v - 1) / 12));
end

function v = std_(x)
if numel(x) < 2, v = 0; return; end
v = sqrt(sum((x - mean(x)).^2) / (numel(x) - 1));
end

function s = modal_share(d)
% Fraction of entries equal to the most frequent value (exact equality,
% after rounding away pure floating-point noise).
d = round(d(:) * 1e10) / 1e10;
u = unique(d);
best = 0;
for k = 1:numel(u)
    best = max(best, sum(d == u(k)));
end
s = best / numel(d);
end

function r = max_run(v)
% Longest run of exactly repeated consecutive values.
r = 1;  cur = 1;
for i = 2:numel(v)
    if v(i) == v(i - 1), cur = cur + 1; else, cur = 1; end
    r = max(r, cur);
end
end
