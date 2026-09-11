function cases = experiment_grid_spec()
% PURPOSE
% -------
% The EXPLORATORY experiment grid: the cells over which the study asks
% "WHEN is block adaptation useful?".  Each cell names a DGP and a set of
% configuration overrides; montecarlo/run_experiment_grid.m runs them.
%
% DESIGN OF THE GRID (and why it is not a full factorial)
% -------------------------------------------------------
% Four things plausibly decide whether a local escape from the VAR prior
% pays for itself:
%   (a) MISSPECIFICATION STRENGTH -- how far the VAR-implied prior centre
%       is from the truth.  Too weak and there is nothing to find; too
%       strong and every estimator that shrinks to the VAR loses, so the
%       adaptive one wins for a reason that has nothing to do with
%       locality.
%   (b) SAMPLE SIZE -- adaptation buys bias reduction at the cost of
%       variance, and the variance cost is what shrinks with T.
%   (c) FITTED LAG ORDER -- p = 1 makes the VAR(2)-based prior wrong in
%       more places; p = 4 nests the truth of the sparse DGP entirely, so
%       there is nothing left to escape from.  This is the sharpest test
%       of "does the diagnostic point at a real defect".
%   (d) SPARSITY of the misspecification -- correct / sparse (one entry,
%       one equation) / intermediate (one column, all equations) / dense
%       (VARMA, every block of every equation).  This is the axis the
%       whole thesis is about.
% A full factorial over these is 4 x 3 x 3 x 4 = 144 cells, which at any
% honest chain length is out of reach.  The grid below is instead a
% ONE-FACTOR-AT-A-TIME sweep around a baseline (sparse DGP, T = 200,
% p = 2, a31 = -0.30), plus the full sparsity axis at the baseline: 13
% cells, each cheap enough to run.  A one-factor sweep cannot identify
% interactions, and this is stated wherever the grid is reported.
%
% EXPLORATORY MEANS EXPLORATORY
% -----------------------------
% The grid runs at the 'grid' preset (R = 40, H = 12, shorter chains).
% The standard error of an RMSE at R = 40 is roughly 11% of the RMSE
% itself, so a grid cell can locate a LARGE, SYSTEMATIC effect but must
% never be quoted as a precise number.  Its purpose is to decide where to
% spend the R = 500 budget.
%
% OUTPUTS
% -------
% cases : (1 x n) struct array with fields
%   .name      short identifier, used for filenames
%   .dgp       DGP name for run_montecarlo
%   .axis      which question this cell belongs to
%   .label     human-readable description of the cell
%   .overrides cfg overrides passed to mc_preset

c = {};

% --- (d) the sparsity axis, at the baseline ---------------------------
c{end+1} = mk('sparsity_correct', 'correct', 'sparsity', ...
    'correctly specified VAR(2): nothing to escape from', struct());
c{end+1} = mk('sparsity_sparse', 'sparse', 'sparsity', ...
    'sparse: one wrong entry A3(3,1) = -0.30, one equation', struct());
c{end+1} = mk('sparsity_interm', 'intermediate', 'sparsity', ...
    'intermediate: A3(:,1) wrong in every equation, one block', struct());
c{end+1} = mk('sparsity_dense', 'dense', 'sparsity', ...
    'dense: VARMA(2,1), every block of every equation wrong', struct());

% --- (a) misspecification strength (sparse DGP) -----------------------
for a = [-0.15, -0.45, -0.60]
    c{end+1} = mk(sprintf('strength_a%03d', round(-100 * a)), 'sparse', ...
        'strength', sprintf('sparse DGP with A3(3,1) = %.2f', a), ...
        struct('dgp', struct('sparse_a31', a)));
end
% and a weaker/stronger dense case, to see whether the dense conclusion
% is about density or simply about the size of the conflict
c{end+1} = mk('strength_dense_half', 'dense', 'strength', ...
    'dense DGP with the MA matrix halved', ...
    struct('dgp', struct('dense_scale', 0.5)));

% --- (b) sample size (sparse DGP) --------------------------------------
for T = [100, 400]
    c{end+1} = mk(sprintf('T%03d', T), 'sparse', 'sample_size', ...
        sprintf('sparse DGP at T = %d', T), struct('T', T));
end

% --- (c) fitted lag order (sparse DGP; truth is a VAR(3)) --------------
for p = [1, 3, 4]
    if p >= 3
        note = sprintf(['fitted p = %d NESTS the true VAR(3): the prior ' ...
                        'centre is correct, so any escape is a false positive'], p);
    else
        note = sprintf('fitted p = %d: the VAR prior is wrong in more places', p);
    end
    c{end+1} = mk(sprintf('p%d', p), 'sparse', 'lag_order', note, ...
        struct('p', p));
end

cases = [c{:}];
end

% =====================================================================
function s = mk(name, dgp, axis, label, overrides)
s = struct('name', name, 'dgp', dgp, 'axis', axis, 'label', label, ...
           'overrides', overrides);
end
