function cfg = mc_preset(name, overrides)
% PURPOSE
% -------
% Named Monte Carlo configurations, so that every reported experiment is
% produced by a preset with a NAME rather than by ad-hoc edits to a
% script.  A stored result records the preset it came from (via
% mc.meta.preset), which together with mc.meta is enough to re-run it.
%
% PRESETS
% -------
% 'final'  : the headline design.  cfg.mode = 'fmar', K = 3, p = 2,
%            T = 200, H = 20, R = 500, all five estimators, chains of
%            300 + 700 draws, cfg.fmar.h1_mode = 'lp' so the h = 1
%            comparison is like-for-like.  This is the configuration
%            whose numbers the thesis quotes.
% 'grid'   : the EXPLORATORY design used to sweep misspecification
%            strength, sample size, lag order and sparsity.  Same
%            estimators and the same machinery, but R = 40, H = 12 and
%            chains of 200 + 400 draws.  Cheap enough to run a dozen
%            cells; NOT precise enough to quote a single RMSE from.  Its
%            job is to say WHERE to spend the R = 500 budget.
% 'legacy' : reproduces the pre-existing four-estimator results exactly:
%            cfg.fmar.h1_mode = 'bvar' (the FMAR convention), no pooled
%            estimator, and cfg.blp.point_estimate = 'draw_mean' (the
%            IRF averaged over draws rather than Rao-Blackwellised).
%            Kept so the earlier numbers in results/ stay reproducible.
% 'smoke'  : a few-minute version of 'final' for interface checks.
%
% MONTE CARLO PRECISION (read before quoting a difference)
% --------------------------------------------------------
% With R replications, the standard error of an estimated RMSE is
% roughly rmse / sqrt(2R): about 3.2% at R = 500 and about 11% at
% R = 40.  Differences between estimators are more precise than that
% because they are computed on the SAME simulated samples (common random
% numbers: every estimator sees the same data in every replication), but
% a grid cell at R = 40 still cannot resolve a 5% RMSE gap.
%
% INPUTS
% ------
% name      : 'final' | 'grid' | 'legacy' | 'smoke'.
% overrides : OPTIONAL struct of dotted field overrides applied last,
%             e.g. struct('T', 400, 'dgp', struct('sparse_a31', -0.15)).
%             Nested structs are merged field by field.
%
% OUTPUTS
% -------
% cfg : configuration struct, with cfg.mc.preset = name.

if nargin < 2, overrides = struct(); end

cfg = default_config();
cfg.mode = 'fmar';
cfg.fmar.h1_mode = 'lp';        % like-for-like at h = 1
cfg.mc.estimators = {'lp', 'var', 'global', 'block', 'pooled'};

switch name
    case 'final'
        cfg.K = 3;  cfg.p = 2;  cfg.T = 200;  cfg.H = 20;
        cfg.mc.n_rep = 500;
        cfg.mc.gibbs_n_burn = 300;
        cfg.mc.gibbs_n_keep = 700;
        cfg.fmar.n_niw_draws = 500;
    case 'grid'
        cfg.K = 3;  cfg.p = 2;  cfg.T = 200;  cfg.H = 12;
        cfg.mc.n_rep = 40;
        cfg.mc.gibbs_n_burn = 200;
        cfg.mc.gibbs_n_keep = 400;
        cfg.fmar.n_niw_draws = 200;
    case 'legacy'
        cfg.K = 3;  cfg.p = 2;  cfg.T = 200;  cfg.H = 20;
        cfg.mc.n_rep = 500;
        cfg.mc.gibbs_n_burn = 300;
        cfg.mc.gibbs_n_keep = 700;
        cfg.fmar.n_niw_draws = 500;
        cfg.fmar.h1_mode = 'bvar';
        cfg.mc.estimators = {'lp', 'var', 'global', 'block'};
        cfg.blp.point_estimate = 'draw_mean';
    case 'smoke'
        cfg.K = 3;  cfg.p = 2;  cfg.T = 140;  cfg.H = 5;
        cfg.mc.n_rep = 3;
        cfg.mc.gibbs_n_burn = 50;
        cfg.mc.gibbs_n_keep = 100;
        cfg.fmar.n_niw_draws = 50;
    otherwise
        error('mc_preset: unknown preset "%s".', name);
end

cfg = merge_struct(cfg, overrides);
cfg.mc.preset = name;
end

% =====================================================================
function a = merge_struct(a, b)
% Recursive field-by-field merge: b wins wherever it defines something.
if ~isstruct(b), return; end
f = fieldnames(b);
for k = 1:numel(f)
    if isstruct(b.(f{k})) && isfield(a, f{k}) && isstruct(a.(f{k}))
        a.(f{k}) = merge_struct(a.(f{k}), b.(f{k}));
    else
        a.(f{k}) = b.(f{k});
    end
end
end
