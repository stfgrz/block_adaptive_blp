function [key, info] = ea_design_key(ds_or_varnames, cfg, extra)
% EA_DESIGN_KEY  Canonical identifier of a tau-diagnostic design.
%
% PURPOSE
% -------
% tau_bar is measured relative to lambda_h, which depends on the system,
% the sample, the lag order, the horizon window and the prior options, so a
% null calibration is only valid for the design it was simulated under
% (DESIGN.md Sec. 7).  This function renders every design-relevant
% choice into ONE canonical string and appends a short hash of it.  The
% null file stores the key; ea_apply_protocol refuses a null whose key
% differs from the run's, so no threshold can be read against the wrong
% design by accident.
%
% INPUTS
% ------
% ds_or_varnames : dataset struct (uses .varnames, size(.Y, 1) and .isrw)
%                  or a cellstr of variable names.
% cfg            : config struct; uses .p, .H, .fmar.h1_mode,
%                  .fmar.psi_floor, .blocks.fixed_tau, .blocks.fixed_tau_mask
%                  (optional), .fmar.isrw (if the dataset has no .isrw) and,
%                  only when extra.include_chains is true, .gibbs.n_burn and
%                  .gibbs.n_keep.
% extra          : optional struct
%     .h_early        [h1 h2] window of the tau_bar statistic; defaults to
%                     cfg.null.h_early when present, else an error
%     .T              sample length, overrides size(ds.Y, 1) / cfg.T
%     .isrw           prior-centre vector, overrides ds.isrw / cfg.fmar.isrw
%     .include_chains false (default): the sampler chain lengths are NOT
%                     part of the key -- the null deliberately uses shorter
%                     chains than the real-data run
%     .ident          identification label (proxy vs recursive, instrument
%                     name); recorded in info.ident but NOT in the key,
%                     because tau does not depend on it
%
% OUTPUTS
% -------
% key  : char, '<canonical string>;hash=<8 hex chars>'
% info : struct of the fields that entered the string (.varnames, .T, .p,
%        .H, .h_early, .h1_mode, .psi_floor, .fixed_tau, .fixed_tau_mask,
%        .isrw, .include_chains, .chains), plus .canonical, .hash and
%        .ident (empty unless given)
%
% NOTES
% -----
% * The hash is 32-bit FNV-1a computed in plain double arithmetic (exact:
%   the 16-bit split keeps every intermediate below 2^53), so it is
%   identical in MATLAB and Octave and needs no Java or toolbox.
% * The canonical string itself is kept in the key, so a mismatch can be
%   read off by eye; the hash only makes the key short to compare.

if nargin < 3 || isempty(extra), extra = struct(); end

% --- variable names, T, isrw ------------------------------------------------
if isstruct(ds_or_varnames)
    ds = ds_or_varnames;
    assert(isfield(ds, 'varnames'), 'ea_design_key: dataset has no .varnames');
    varnames = ds.varnames;
    T = [];  if isfield(ds, 'Y'), T = size(ds.Y, 1); end
    isrw = [];  if isfield(ds, 'isrw'), isrw = ds.isrw; end
else
    varnames = ds_or_varnames;
    T = [];  isrw = [];
end
if ~iscell(varnames), varnames = cellstr(varnames); end
varnames = reshape(varnames, 1, []);
if isfield(extra, 'T') && ~isempty(extra.T), T = extra.T; end
if isempty(T) && isfield(cfg, 'T'), T = cfg.T; end
assert(~isempty(T), 'ea_design_key: sample length T unknown (pass a dataset, extra.T or cfg.T)');
if isfield(extra, 'isrw') && ~isempty(extra.isrw), isrw = extra.isrw; end
if isempty(isrw) && isfield(cfg, 'fmar') && isfield(cfg.fmar, 'isrw'), isrw = cfg.fmar.isrw; end
assert(~isempty(isrw), 'ea_design_key: isrw unknown (pass a dataset, extra.isrw or cfg.fmar.isrw)');

% --- lag order, horizons, window -------------------------------------------
assert(isfield(cfg, 'p'), 'ea_design_key: cfg needs .p');
if ~isfield(cfg, 'H'), cfg.H = NaN; end
% NOTE: H is recorded in info but is NOT part of the key.  The tau_bar
% statistic lives on h <= h_early(2); the independent estimator's tau at
% those horizons does not depend on H, and the null is deliberately run at
% H = h_early(2) while the real-data run uses a longer H for the IRFs.  For
% the horizon-pooled estimator the log-tau random walk spans all H, so a
% null at a shorter H is an approximation, stated in the protocol output.
if isfield(extra, 'h_early') && ~isempty(extra.h_early)
    h_early = extra.h_early;
elseif isfield(cfg, 'null') && isfield(cfg.null, 'h_early')
    h_early = cfg.null.h_early;
else
    error('ea_design_key: h_early unknown (pass extra.h_early or cfg.null.h_early)');
end
assert(numel(h_early) == 2, 'ea_design_key: h_early must be a 2-vector');

% --- prior options ------------------------------------------------------------
assert(isfield(cfg, 'fmar') && isfield(cfg.fmar, 'h1_mode'), 'ea_design_key: cfg.fmar.h1_mode missing');
h1_mode = cfg.fmar.h1_mode;
psi_floor = false;
if isfield(cfg.fmar, 'psi_floor'), psi_floor = logical(cfg.fmar.psi_floor); end
fixed_tau = [];
if isfield(cfg, 'blocks') && isfield(cfg.blocks, 'fixed_tau'), fixed_tau = cfg.blocks.fixed_tau; end
fixed_tau = sort(reshape(fixed_tau, 1, []));
mask_str = '-';
if isfield(cfg, 'blocks') && isfield(cfg.blocks, 'fixed_tau_mask') && ...
        ~isempty(cfg.blocks.fixed_tau_mask) && any(cfg.blocks.fixed_tau_mask(:))
    m = logical(cfg.blocks.fixed_tau_mask);
    mask_str = sprintf('%dx%d:%s', size(m, 1), size(m, 2), char('0' + reshape(m', 1, [])));
end

% --- chains (optional) --------------------------------------------------------
include_chains = isfield(extra, 'include_chains') && logical(extra.include_chains);
chains = [];
if include_chains
    assert(isfield(cfg, 'gibbs'), 'ea_design_key: cfg.gibbs missing but include_chains = true');
    chains = [cfg.gibbs.n_burn, cfg.gibbs.n_keep];
end

% --- canonical string ---------------------------------------------------------
canon = sprintf('vars=%s;T=%d;p=%d;h_early=%d-%d;h1_mode=%s;psi_floor=%d;fixed_tau=[%s];fixed_tau_mask=%s;isrw=%s', ...
                strjoin(varnames, '|'), T, cfg.p, h_early(1), h_early(2), ...
                h1_mode, psi_floor, strjoin(arrayfun(@(v) sprintf('%d', v), fixed_tau, ...
                'UniformOutput', false), ','), mask_str, ...
                char('0' + (reshape(isrw, 1, []) ~= 0)));
if include_chains
    canon = sprintf('%s;chains=%d+%d', canon, chains(1), chains(2));
end
hash = fnv1a_hex(canon);
key = sprintf('%s;hash=%s', canon, hash);

info = struct();
info.varnames = varnames;  info.T = T;  info.p = cfg.p;  info.H = cfg.H;
info.h_early = reshape(h_early, 1, 2);  info.h1_mode = h1_mode;
info.psi_floor = psi_floor;  info.fixed_tau = fixed_tau;
info.fixed_tau_mask = mask_str;  info.isrw = reshape(isrw, 1, []);
info.include_chains = include_chains;  info.chains = chains;
info.ident = '';
if isfield(extra, 'ident'), info.ident = extra.ident; end
info.canonical = canon;  info.hash = hash;
end

% -------------------------------------------------------------------------
function hx = fnv1a_hex(s)
% 32-bit FNV-1a over the bytes of s, in exact double arithmetic.
b = double(uint8(s));
h = 2166136261;                 % offset basis
prime = 16777619;
for k = 1:numel(b)
    h = bitxor(h, b(k));
    lo = mod(h, 65536);  hi = (h - lo) / 65536;
    h = mod(mod(hi * prime, 65536) * 65536 + lo * prime, 4294967296);
end
hx = lower(dec2hex(h, 8));
end
