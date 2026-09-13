function [fixed_tau_blocks, held_mask, fixed_by_eq, eq_set] = resolve_blp_cell_options(cfg, K, G, caller)
% PURPOSE
% -------
% Resolve the three options that restrict WHICH (equation, block) cells of
% a block-adaptive BLP are sampled.  Shared by estimate_blp_blockadaptive
% and estimate_blp_blockpooled so that both read them identically:
%
%   cfg.blocks.fixed_tau       blocks held at tau = 1 in EVERY equation;
%   cfg.blocks.fixed_tau_mask  (K x G) logical: cells held at tau = 1;
%   cfg.blp.equations          equations to estimate ([] = all).
%
% MODEL
% -----
% Cell (i, g) is held at tau_{i,g,h} = 1 for every h iff
%     fixed_tau_mask(i, g)  OR  g in fixed_tau          (union).
% The samplers take one block list per call (opts.fixed_blocks), and both
% estimators call their sampler once per equation, so the union mask is
% resolved here into one list per equation.
%
% INPUTS
% ------
% cfg    : configuration struct; any of the fields may be missing/empty.
% K, G   : number of equations and of blocks (G = K under 'per_variable').
% caller : name used in error messages.
%
% OUTPUTS
% -------
% fixed_tau_blocks : (1 x n) validated cfg.blocks.fixed_tau ([] = none).
% held_mask        : (K x G) logical, the effective union mask.
% fixed_by_eq      : (K x 1) cell; fixed_by_eq{i} = row vector of the
%                    blocks held in equation i (the sampler's
%                    opts.fixed_blocks), [] when none.
% eq_set           : (1 x n_eq) sorted equation indices to estimate.
%
% NOTES
% -----
% With every option empty: fixed_tau_blocks = [], held_mask all false,
% every fixed_by_eq{i} = [] and eq_set = 1:K, so the samplers receive
% exactly the inputs they received before these options existed (an empty
% opts.fixed_blocks is treated as "absent" by both samplers).

if nargin < 4, caller = 'resolve_blp_cell_options'; end

% --- cfg.blocks.fixed_tau: whole blocks -----------------------------------
fixed_tau_blocks = [];
if isfield(cfg, 'blocks') && isfield(cfg.blocks, 'fixed_tau') && ~isempty(cfg.blocks.fixed_tau)
    fixed_tau_blocks = cfg.blocks.fixed_tau(:)';
    assert(all(fixed_tau_blocks >= 1 & fixed_tau_blocks <= G), ...
        '%s: cfg.blocks.fixed_tau must index blocks 1..%d.', caller, G);
end
held_mask = false(K, G);
held_mask(:, fixed_tau_blocks) = true;

% --- cfg.blocks.fixed_tau_mask: single cells (union with the above) ------
if isfield(cfg, 'blocks') && isfield(cfg.blocks, 'fixed_tau_mask') && ~isempty(cfg.blocks.fixed_tau_mask)
    M = cfg.blocks.fixed_tau_mask;
    assert(isequal(size(M), [K, G]), ...
        '%s: cfg.blocks.fixed_tau_mask must be %d x %d (equation x block), got %d x %d.', ...
        caller, K, G, size(M, 1), size(M, 2));
    assert(islogical(M) || all(M(:) == 0 | M(:) == 1), ...
        '%s: cfg.blocks.fixed_tau_mask must be logical (or 0/1).', caller);
    held_mask = held_mask | logical(M);
end

fixed_by_eq = cell(K, 1);
for i = 1:K
    fb = find(held_mask(i, :));
    if isempty(fb), fb = []; end            % normalise 1x0 to []
    fixed_by_eq{i} = fb;
end

% --- cfg.blp.equations: equation subset ------------------------------------
eq_set = 1:K;
if isfield(cfg, 'blp') && isfield(cfg.blp, 'equations') && ~isempty(cfg.blp.equations)
    eq_set = unique(cfg.blp.equations(:)');
    assert(all(eq_set >= 1 & eq_set <= K & eq_set == round(eq_set)), ...
        '%s: cfg.blp.equations must index equations 1..%d.', caller, K);
end
end
