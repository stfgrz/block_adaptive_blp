function mc = merge_montecarlo(chunks)
% PURPOSE
% -------
% Glue together Monte Carlo CHUNKS produced by run_montecarlo with
% cfg.mc.rep_range set, so one experiment can be spread over several
% processes (this project's machine has 4 cores and the samplers are the
% bottleneck) and still yield a single result that is BIT-IDENTICAL to
% the serial run.
%
% WHY THIS IS EXACT
% -----------------
% Replication r's random stream is seeded by cfg.seed + 100000*dgp_id + r
% and nothing else: no state carries over between replications.  A chunk
% covering replications r0..r1 therefore contains exactly the columns the
% full run would have produced for those replications, and merging is
% concatenation in replication order -- no re-randomisation, no
% approximation.  The function ASSERTS that the chunks agree on the
% design (DGP, mode, estimators, dimensions, master seed) and that the
% union of their replication indices is exactly 1..R with no gaps and no
% duplicates, so a silently incomplete merge is impossible.
%
% INPUTS
% ------
% chunks : cell array of mc structs (or a struct array), in any order.
%
% OUTPUTS
% -------
% mc : merged struct with the same fields as a serial run,
%      .rep_index = 1:R sorted, .meta.is_complete set accordingly, and
%      .meta.merged_from recording the chunk ranges.
%
% NOTES
% -----
% Deterministic.  Timing fields are summed; .meta.created is taken from
% the first chunk.

if isstruct(chunks), chunks = num2cell(chunks); end
assert(iscell(chunks) && ~isempty(chunks), ...
    'merge_montecarlo: expects a non-empty cell array of mc structs.');

ref = chunks{1};
key = {'dgp_name', 'est_names', 'theta_true', 'misspec_block'};
for c = 2:numel(chunks)
    for k = 1:numel(key)
        assert(isequaln(ref.(key{k}), chunks{c}.(key{k})), ...
            'merge_montecarlo: chunks disagree on %s.', key{k});
    end
    assert(ref.meta.master_seed == chunks{c}.meta.master_seed, ...
        'merge_montecarlo: chunks used different master seeds.');
    assert(ref.meta.R == chunks{c}.meta.R, ...
        'merge_montecarlo: chunks belong to designs with different R.');
    assert(strcmp(ref.meta.mode, chunks{c}.meta.mode), ...
        'merge_montecarlo: chunks used different modes.');
end

% --- replication bookkeeping -------------------------------------------
idx = [];  ranges = {};
for c = 1:numel(chunks)
    idx = [idx, chunks{c}.rep_index(:)'];  %#ok<AGROW>
    ranges{end + 1} = sprintf('%d-%d', chunks{c}.rep_index(1), ...
                              chunks{c}.rep_index(end));  %#ok<AGROW>
end
[sorted, order] = sort(idx);
assert(numel(unique(sorted)) == numel(sorted), ...
    'merge_montecarlo: overlapping replication indices.');
R_total = ref.meta.R;
complete = isequal(sorted(:)', 1:R_total);
if ~complete
    missing = setdiff(1:R_total, sorted);
    fprintf(2, ['merge_montecarlo: WARNING -- merged result is INCOMPLETE ' ...
                '(%d of %d replications; first missing: %d).\n'], ...
            numel(sorted), R_total, missing(1));
end

mc = ref;
mc.theta = cat_reps(chunks, 'theta', order);
mc.lo    = cat_reps(chunks, 'lo', order);
mc.hi    = cat_reps(chunks, 'hi', order);
if isfield(ref, 'lambda'),   mc.lambda   = cat_reps(chunks, 'lambda', order);   end
if isfield(ref, 'lo_post'),  mc.lo_post  = cat_reps(chunks, 'lo_post', order);  end
if isfield(ref, 'hi_post'),  mc.hi_post  = cat_reps(chunks, 'hi_post', order);  end
if isfield(ref, 'tau_mean'), mc.tau_mean = cat_reps(chunks, 'tau_mean', order); end

if isfield(ref, 'tau') && isstruct(ref.tau)
    f = fieldnames(ref.tau);
    for k = 1:numel(f)
        sub = fieldnames(ref.tau.(f{k}));
        for j = 1:numel(sub)
            A = [];
            for c = 1:numel(chunks)
                A = cat(ndims(ref.tau.(f{k}).(sub{j})), A, ...
                        chunks{c}.tau.(f{k}).(sub{j}));
            end
            mc.tau.(f{k}).(sub{j}) = reorder_last(A, order);
        end
    end
end

% Per-replication diagnostics are column vectors.
df = fieldnames(ref.diag);
for k = 1:numel(df)
    v = [];
    for c = 1:numel(chunks), v = [v; chunks{c}.diag.(df{k})(:)]; end  %#ok<AGROW>
    if numel(v) == numel(idx)
        mc.diag.(df{k}) = v(order);
    end
end

sd = [];
for c = 1:numel(chunks), sd = [sd; chunks{c}.seeds(:)]; end  %#ok<AGROW>
mc.seeds = sd(order);
mc.rep_index = sorted(:)';
mc.meta.rep_index   = mc.rep_index;
mc.meta.R_done      = numel(sorted);
mc.meta.is_complete = complete;
mc.meta.merged_from = ranges;
el = 0;
for c = 1:numel(chunks), el = el + chunks{c}.meta.elapsed_seconds; end
mc.meta.elapsed_seconds = el;
mc.meta.seeds = mc.seeds(:)';
end

% =====================================================================
function A = cat_reps(chunks, field, order)
A = [];
nd = ndims(chunks{1}.(field));
for c = 1:numel(chunks)
    A = cat(nd, A, chunks{c}.(field));
end
A = reorder_last(A, order);
end

function A = reorder_last(A, order)
% Reorder along the LAST dimension (the replication dimension).
nd = ndims(A);
switch nd
    case 2, A = A(:, order);
    case 3, A = A(:, :, order);
    case 4, A = A(:, :, :, order);
    otherwise
        error('reorder_last: unsupported %d-d array.', nd);
end
end
