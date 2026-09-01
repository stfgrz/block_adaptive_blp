function test_isrw_vector()
% PURPOSE
% -------
% Guard the vectorised cfg.fmar.isrw prior mean in estimate_bvar_niw.
%
% The empirical application (Ch. 7) mixes a white-noise-ish policy
% surprise with persistent levels and therefore needs a PER-VARIABLE
% random-walk flag, cfg.fmar.isrw = [0 1 1 1 1].  The original code read
% the flag as `if cfg.fmar.isrw`, which in MATLAB/Octave is all() of the
% elements: a vector starting with 0 evaluates FALSE and would silently
% drop the random-walk prior centre for EVERY variable -- a wrong prior
% with no error message.  This test pins down the fixed behaviour.
%
% CHECKS
% ------
% 1. BACKWARD COMPATIBILITY.  Scalar isrw = false / true must give
%    bit-identical results to the pre-patch code, i.e. prior mean 0 and
%    eye(K) on the first own lag respectively.  Checked by comparing the
%    scalar runs against explicit all-zeros / all-ones vectors.
% 2. VECTOR SEMANTICS.  isrw = [0 1 1 ...] must differ from both scalar
%    cases (it is genuinely a mixed centre) and must equal a hand-built
%    reference fit obtained by the same estimator with the equivalent
%    vector transposed / logical / double forms.
% 3. THE SILENT-FAILURE CASE.  isrw = [0 1 1] must NOT reproduce the
%    isrw = false fit; if it does, the truthiness bug is back.
% 4. INPUT VALIDATION.  A wrong-length isrw must error, not silently
%    mis-broadcast.
%
% INPUTS / OUTPUTS
% ----------------
% None; asserts and prints a PASS line.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'), fullfile(root, 'dgp'), ...
        fullfile(root, 'estimators'), fullfile(root, 'priors'), ...
        fullfile(root, 'samplers'), fullfile(root, 'utils'));

cfg = default_config();
cfg.mode = 'fmar';
cfg.p = 2;  cfg.H = 4;  cfg.T = 160;
cfg.fmar.n_niw_draws = 20;          % bands are irrelevant here; keep it fast

rng(9091, 'twister');
d = simulate_var_dgp(cfg);
Y = d.Y;
K = size(Y, 2);

fit = @(flag) run_bvar(Y, cfg, flag);

% --- 1. scalar backward compatibility -----------------------------------
B_false_s = fit(false);
B_false_v = fit(zeros(1, K));
assert(max(abs(B_false_s(:) - B_false_v(:))) == 0, ...
       'test_isrw_vector: scalar false must equal an all-zeros vector.');

B_true_s = fit(true);
B_true_v = fit(ones(1, K));
assert(max(abs(B_true_s(:) - B_true_v(:))) == 0, ...
       'test_isrw_vector: scalar true must equal an all-ones vector.');
fprintf('  scalar isrw backward compatibility: OK\n');

% The two scalar centres must actually differ, otherwise checks 2-3 are
% vacuous (e.g. if the prior mean were ignored altogether).
assert(max(abs(B_true_s(:) - B_false_s(:))) > 1e-8, ...
       'test_isrw_vector: isrw true vs false must change the fit.');

% --- 2. vector semantics -------------------------------------------------
v = [0, ones(1, K - 1)];
B_mix = fit(v);
B_col = fit(v(:));                   % column orientation
B_log = fit(logical(v));             % logical type
assert(max(abs(B_mix(:) - B_col(:))) == 0, ...
       'test_isrw_vector: row and column isrw must agree.');
assert(max(abs(B_mix(:) - B_log(:))) == 0, ...
       'test_isrw_vector: logical and double isrw must agree.');
fprintf('  vector isrw orientation/type invariance: OK\n');

% --- 3. the silent-failure case -----------------------------------------
% Pre-patch, `if [0 1 1]` was FALSE, so this fit collapsed onto isrw=false.
gap_false = max(abs(B_mix(:) - B_false_s(:)));
gap_true  = max(abs(B_mix(:) - B_true_s(:)));
assert(gap_false > 1e-8, ...
       ['test_isrw_vector: isrw = [0 1 ... 1] reproduced the isrw = false ' ...
        'fit -- the `if isrw` truthiness bug is back.']);
assert(gap_true > 1e-8, ...
       'test_isrw_vector: isrw = [0 1 ... 1] must differ from isrw = true.');
fprintf('  mixed centre distinct from both scalar cases (gaps %.2e / %.2e): OK\n', ...
        gap_false, gap_true);

% The mixed centre must differ from the all-ones centre ONLY through
% variable 1's own-lag prior mean, so the difference is a rank-limited
% object: check the fit moves in the expected direction (variable 1's
% first own-lag coefficient shrinks towards 0 rather than towards 1).
i1 = 1 + 1;                          % row of y_1(t) in the regressor block
assert(B_mix(i1, 1) < B_true_s(i1, 1), ...
       'test_isrw_vector: white-noise centre must pull A1(1,1) below the RW centre.');
fprintf('  A1(1,1) pulled towards white noise (%.4f < %.4f): OK\n', ...
        B_mix(i1, 1), B_true_s(i1, 1));

% --- 4. input validation -------------------------------------------------
threw = false;
try
    fit(ones(1, K + 1));
catch
    threw = true;
end
assert(threw, 'test_isrw_vector: a wrong-length isrw must raise an error.');
fprintf('  wrong-length isrw rejected: OK\n');

fprintf('PASS: test_isrw_vector\n');
end

% -------------------------------------------------------------------------
function B = run_bvar(Y, cfg, flag)
% Posterior-mean coefficient matrix at a given isrw flag.  The seed is
% reset so the (irrelevant here) band draws cannot make runs differ.
cfg.fmar.isrw = flag;
rng(4321, 'twister');
bv = estimate_bvar_niw(Y, cfg);
B  = bv.B;
end
