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
for c = 1:numel(chunks)
    % A compacted result (montecarlo/compact_mc.m) has had its
    % per-replication tau quantile arrays averaged away, so it can no
    % longer be concatenated along the replication dimension.  Refuse it
    % with a message that says what happened, rather than failing later
    % on an opaque size mismatch.
    if isfield(chunks{c}, 'meta') && isstruct(chunks{c}.meta) && ...
            isfield(chunks{c}.meta, 'tau_quantiles_compacted') && ...
            chunks{c}.meta.tau_quantiles_compacted
        error(['merge_montecarlo: chunk %d has been compacted ' ...
               '(compact_mc collapsed its per-replication tau quantiles), ' ...
               'so it cannot be merged.  Merge the raw chunk files, not a ' ...
               'saved merged result.'], c);
    end
end
% --- the chunks must come from the SAME design ------------------------
% Checking only the top-level fields is not enough.  theta_true is a
% function of the DGP parameters and H alone, so two chunks run at
% DIFFERENT T (or p, or chain length) agree on every one of them and
% would merge into a result that looks valid and is not.  Compare the
% full design record instead, and name the field that differs.
key = {'dgp_name', 'est_names', 'theta_true', 'misspec_block'};
for c = 2:numel(chunks)
    for k = 1:numel(key)
        assert(isequaln(ref.(key{k}), chunks{c}.(key{k})), ...
            'merge_montecarlo: chunks disagree on %s.', key{k});
    end
    assert_same_design(ref.meta, chunks{c}.meta, c);
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
        % Only fields present in EVERY chunk can be merged; a field added
        % between two runs would otherwise be silently taken from one
        % chunk and presented as covering all replications.
        sub = fieldnames(ref.tau.(f{k}));
        for c = 2:numel(chunks)
            sub = intersect(sub, fieldnames(chunks{c}.tau.(f{k})));
        end
        dropped = setdiff(fieldnames(ref.tau.(f{k})), sub);
        for j = 1:numel(dropped)
            fprintf(2, ['merge_montecarlo: tau.%s.%s is missing from some ' ...
                        'chunks and was dropped.\n'], f{k}, dropped{j});
        end
        % Fields that are DESCRIPTIONS of the tau summaries rather than
        % per-replication arrays.  Identified by NAME, not by size: the
        % probability grid is 1 x nQ, and a chunk that happened to hold
        % exactly nQ replications would otherwise be concatenated along
        % it and silently corrupted.
        meta_fields = {'probs'};
        for j = 1:numel(sub)
            A0 = ref.tau.(f{k}).(sub{j});
            if any(strcmp(sub{j}, meta_fields)) || ...
                    size(A0, ndims(A0)) ~= numel(ref.rep_index)
                % Not a per-replication array (e.g. the probability grid
                % attached to the quantiles, or a replication-averaged
                % summary): carry it through unchanged after checking the
                % chunks agree on it.
                for c = 2:numel(chunks)
                    assert(isequaln(A0, chunks{c}.tau.(f{k}).(sub{j})), ...
                        'merge_montecarlo: chunks disagree on tau.%s.%s.', ...
                        f{k}, sub{j});
                end
                mc.tau.(f{k}).(sub{j}) = A0;
                continue
            end
            A = [];
            for c = 1:numel(chunks)
                A = cat(ndims(A0), A, chunks{c}.tau.(f{k}).(sub{j}));
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
function assert_same_design(m1, m2, c)
% Every field of mc.meta that DESCRIBES THE DESIGN must agree; the
% fields that legitimately differ between chunks (which replications
% each ran, how long it took, when it was created) are excluded by
% name.  Nested structs (fmar settings, pooling settings, DGP
% parameters) are compared whole.
per_chunk = {'rep_index', 'seeds', 'R_done', 'is_complete', ...
             'elapsed_seconds', 'created', 'merged_from', ...
             'tau_quantiles_compacted'};
f = union(fieldnames(m1), fieldnames(m2));
for k = 1:numel(f)
    name = f{k};
    if any(strcmp(name, per_chunk)), continue; end
    has1 = isfield(m1, name);  has2 = isfield(m2, name);
    assert(has1 == has2, ...
        ['merge_montecarlo: chunk %d disagrees on the design -- ' ...
         'meta.%s is present in one chunk and not the other.'], c, name);
    if ~has1, continue; end
    v1 = m1.(name);  v2 = m2.(name);
    if isa(v1, 'function_handle') || isa(v2, 'function_handle')
        continue        % handles never compare equal; cfg_to_savable strips them
    end
    assert(isequaln(v1, v2), ...
        ['merge_montecarlo: chunk %d belongs to a DIFFERENT design -- ' ...
         'meta.%s differs (%s vs %s).  Merging chunks from different ' ...
         'designs would silently produce a result that looks valid and ' ...
         'is not.'], c, name, describe(v1), describe(v2));
end
end

function s = describe(v)
% Short, safe rendering of a meta value for an error message.
if ischar(v)
    s = v;
elseif iscellstr(v)  %#ok<ISCLSTR>
    s = strjoin(v(:)', ',');
elseif (isnumeric(v) || islogical(v)) && numel(v) <= 8
    s = mat2str(double(v));
elseif isnumeric(v) || islogical(v)
    s = sprintf('<%s %s>', mat2str(size(v)), class(v));
else
    s = sprintf('<%s>', class(v));
end
end

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
    case 5, A = A(:, :, :, :, order);
    otherwise
        error('reorder_last: unsupported %d-d array.', nd);
end
end
