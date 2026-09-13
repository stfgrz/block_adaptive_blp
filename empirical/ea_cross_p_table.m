function tab = ea_cross_p_table(spec, out_csv)
% EA_CROSS_P_TABLE  Escape counts of the tau protocol across designs.
%
% PURPOSE
% -------
% One row per (design, scale estimator) summarising the reading protocol
% (ea_apply_protocol) of several real-data runs, each against ITS OWN null
% calibration.  Raw tau_bar is never compared across designs: tau is
% measured relative to lambda_h, which differs by p and by system, so only
% null-relative quantities are tabulated -- ratio to q95, null percentiles
% and escape counts (CH7_REDESIGN.md Sec. 7).
%
% INPUTS
% ------
% spec    : struct array, one entry per design
%     .label        e.g. 'lev4 p=2'
%     .result_file  .mat with blpb (and optionally blpp), h_early, varnames
%                   or ds_info.varnames, cfg  (a RUN_EMPIRICAL result)
%     .null_file    .mat written by run_null_calibration for that design
%     .held_blocks  indices of blocks held at tau = 1 ([] = none)
% out_csv : csv path; '' = do not write
%
% OUTPUTS
% -------
% tab : struct array, one row per (design, estimator), fields
%     label, estimator, p, n_free_cells, n_escape_q95, n_escape_holm,
%     n_eq_flag, max_ratio_q95, median_ratio_q95, top1, top2, top3
%     ('eq<-block(ratio)' strings, by ratio_q95 over the free cells),
%     key_status (the design-key check of ea_apply_protocol)
%
% NOTES
% -----
% * The run's design key is rebuilt from the result file (cfg, varnames,
%   h_early, T) and compared with the null's; a mismatch is reported in
%   key_status rather than aborting the table (ea_apply_protocol is called
%   with opts.force = true), so the reader sees WHICH rows rest on a
%   mismatched null.  T is taken from the result file when it carries one
%   (T, ds_info.T, Y, cfg.T) and from the null otherwise.
% * The per-cell csvs are not written here; call ea_apply_protocol for
%   those.

if nargin < 2, out_csv = ''; end
tab = struct('label', {}, 'estimator', {}, 'p', {}, 'n_free_cells', {}, ...
             'n_escape_q95', {}, 'n_escape_holm', {}, 'n_eq_flag', {}, ...
             'max_ratio_q95', {}, 'median_ratio_q95', {}, ...
             'top1', {}, 'top2', {}, 'top3', {}, 'key_status', {});

fprintf('\nea_cross_p_table: %d designs.  Raw tau_bar is NOT comparable across designs;\n', numel(spec));
fprintf('  only ratio_q95, null percentiles and counts are (each design against its own null).\n');

for d = 1:numel(spec)
    sp = spec(d);
    assert(exist(sp.result_file, 'file') == 2, 'ea_cross_p_table: %s not found', sp.result_file);
    assert(exist(sp.null_file, 'file') == 2, 'ea_cross_p_table: %s not found', sp.null_file);
    R = load(sp.result_file);
    N = load(sp.null_file);
    held = [];  if isfield(sp, 'held_blocks'), held = sp.held_blocks; end

    vn = first_of(R, {'varnames', 'ds_info.varnames'}, first_of(N, {'varnames'}, {}));
    assert(~isempty(vn), 'ea_cross_p_table: %s carries no variable names', sp.result_file);
    h_early = first_of(R, {'h_early'}, first_of(N, {'h_early'}, []));
    assert(numel(h_early) == 2, 'ea_cross_p_table: h_early unknown for %s', sp.label);
    p = NaN;
    if isfield(R, 'cfg') && isfield(R.cfg, 'p'), p = R.cfg.p; elseif isfield(N, 'p'), p = N.p; end

    % the run's design key, when the result file lets us build it
    run_key = '';
    if isfield(R, 'cfg')
        T = first_of(R, {'T', 'ds_info.T'}, []);
        if isempty(T) && isfield(R, 'Y'), T = size(R.Y, 1); end
        if isempty(T) && isfield(R.cfg, 'T'), T = R.cfg.T; end
        if isempty(T), T = first_of(N, {'T'}, []); end
        try
            run_key = ea_design_key(vn, R.cfg, struct('h_early', h_early, 'T', T));
        catch err
            fprintf('  [%s] design key not built (%s); null left unchecked\n', sp.label, err.message);
        end
    end

    assert(isfield(R, 'blpb'), 'ea_cross_p_table: %s has no blpb', sp.result_file);
    sets = {'block', R.blpb};
    if isfield(R, 'blpp') && isfield(N, 'pooled'), sets(end + 1, :) = {'pooled', R.blpp}; end

    fprintf('\n[%s] (p = %g, null %s)\n', sp.label, p, sp.null_file);
    prot = ea_apply_protocol(sets, N, vn, h_early, held, '', ...
                             struct('design_key', run_key, 'force', true));

    for e = 1:numel(prot.names)
        nm = prot.names{e};  E = prot.(nm);
        K = numel(vn);
        rf = E.ratio_q95;  rf(:, ~prot.free) = NaN;      % free cells only
        idx = find(~isnan(rf(:)));                       % (sort puts NaN FIRST in 'descend')
        r = rf(idx);
        [rs, o] = sort(r, 'descend');
        ord = idx(o);
        top = {'', '', ''};
        for t = 1:min(3, numel(r))
            [ii, gg] = ind2sub([K, K], ord(t));
            top{t} = sprintf('%s<-%s(%.2f)', vn{ii}, vn{gg}, rs(t));
        end
        row = struct('label', sp.label, 'estimator', nm, 'p', p, ...
                     'n_free_cells', E.n_free_cells, 'n_escape_q95', E.n_escape_q95, ...
                     'n_escape_holm', E.n_escape_holm, 'n_eq_flag', E.n_eq_flag, ...
                     'max_ratio_q95', max_or_nan(r), 'median_ratio_q95', median_or_nan(r), ...
                     'top1', top{1}, 'top2', top{2}, 'top3', top{3}, ...
                     'key_status', prot.key_status);
        tab(end + 1) = row;                                             %#ok<AGROW>
    end
end

% --- print ------------------------------------------------------------------
fprintf('\n%-16s %-7s %3s %5s %6s %6s %6s %8s %8s  %-22s %-22s %-22s %s\n', ...
        'design', 'est', 'p', 'free', 'q95', 'holm', 'eqflag', 'maxrat', 'medrat', ...
        'top1', 'top2', 'top3', 'key');
for k = 1:numel(tab)
    t = tab(k);
    fprintf('%-16s %-7s %3g %5d %6d %6s %6d %8.3f %8.3f  %-22s %-22s %-22s %s\n', ...
            t.label, t.estimator, t.p, t.n_free_cells, t.n_escape_q95, ...
            num2str(t.n_escape_holm), t.n_eq_flag, t.max_ratio_q95, t.median_ratio_q95, ...
            t.top1, t.top2, t.top3, t.key_status);
end

% --- csv ----------------------------------------------------------------------
if ~isempty(out_csv)
    fid = fopen(out_csv, 'w');
    assert(fid > 0, 'ea_cross_p_table: cannot write %s', out_csv);
    fprintf(fid, ['# raw tau_bar is not comparable across designs: only ratio_q95, ' ...
                  'null percentiles and counts are (each design against its own null)\n']);
    fprintf(fid, ['design,scale_estimator,p,n_free_cells,n_escape_q95,n_escape_holm,n_eq_flag,' ...
                  'max_ratio_q95,median_ratio_q95,top1,top2,top3,key_status\n']);
    for k = 1:numel(tab)
        t = tab(k);
        fprintf(fid, '"%s",%s,%g,%d,%d,%s,%d,%.5f,%.5f,%s,%s,%s,%s\n', ...
                strrep(t.label, '"', ''''), t.estimator, t.p, t.n_free_cells, t.n_escape_q95, ...
                num2str(t.n_escape_holm), t.n_eq_flag, t.max_ratio_q95, t.median_ratio_q95, ...
                t.top1, t.top2, t.top3, t.key_status);
    end
    fclose(fid);
    fprintf('  wrote %s\n', out_csv);
end
end

% -------------------------------------------------------------------------
function v = first_of(s, paths, d)
% First present field among dotted paths, else the default d.
v = d;
for k = 1:numel(paths)
    parts = strsplit(paths{k}, '.');
    x = s;  ok = true;
    for j = 1:numel(parts)
        if isstruct(x) && isfield(x, parts{j}), x = x.(parts{j}); else, ok = false; break; end
    end
    if ok && ~isempty(x), v = x; return; end
end
end

function m = max_or_nan(r)
if isempty(r), m = NaN; else, m = max(r); end
end

function m = median_or_nan(r)
if isempty(r), m = NaN; else, m = median(r); end
end
