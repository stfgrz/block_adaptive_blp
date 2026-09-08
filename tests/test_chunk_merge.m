function test_chunk_merge()
% PURPOSE
% -------
% Check the two pieces of machinery that make a long Monte Carlo
% practical AND reproducible on a 4-core machine:
%   * replication CHUNKING (cfg.mc.rep_range) plus
%     montecarlo/merge_montecarlo.m, and
%   * the configuration round trip through utils/cfg_to_savable.m /
%     cfg_from_saved.m, without which a saved cfg loses (Octave: refuses
%     to write) the anonymous function it carries.
%
% WHAT IS CHECKED
% ---------------
% 1. EXACTNESS. A design split into three chunks and merged must be
%    IDENTICAL to the serial run -- every estimate, every interval
%    bound, every tau array and every seed.  This is what licenses
%    running the R = 500 experiments in parallel and reporting them as
%    one experiment.  The comparison uses isequaln because the interval
%    arrays legitimately contain NaN for estimators that do not report a
%    band at a given horizon.
% 2. BOOKKEEPING. The merged struct must know it is complete, must carry
%    the replication indices in sorted order, and must record which
%    chunks it came from.
% 3. REFUSALS. Overlapping chunks must be rejected outright (a silent
%    double-count would bias every metric), and chunks from different
%    designs must be rejected.
% 4. INCOMPLETE MERGES must still work but must be MARKED incomplete,
%    so a partial result cannot be mistaken for a finished one.
% 5. CONFIGURATION ROUND TRIP. cfg_to_savable must strip the function
%    handle, the stripped struct must actually save and load, and
%    cfg_from_saved must restore a handle that computes the same values
%    as the original (checked against the horizon-dependent hyperprior
%    sd rule at every horizon).

fprintf('test_chunk_merge:\n');

cfg = mc_preset('smoke');
cfg.mc.n_rep = 6;
cfg.H = 4;

serial = run_montecarlo(cfg, 'sparse');

ch = cell(1, 3);
ranges = {[1 2], [3 4], [5 6]};
for k = 1:3
    c = cfg;  c.mc.rep_range = ranges{k};
    ch{k} = run_montecarlo(c, 'sparse');
    assert(~ch{k}.meta.is_complete, 'a 2-of-6 chunk claimed to be complete');
end

% --- 1. exactness (chunks deliberately handed over out of order) -------
merged = merge_montecarlo(ch([3 1 2]));
fields = {'theta', 'lo', 'hi', 'lo_post', 'hi_post', 'lambda'};
for k = 1:numel(fields)
    assert(isequaln(serial.(fields{k}), merged.(fields{k})), ...
        'merged %s differs from the serial run', fields{k});
end
assert(isequaln(serial.tau.block.mean,  merged.tau.block.mean), ...
    'merged block tau differs from the serial run');
assert(isequaln(serial.tau.pooled.mean, merged.tau.pooled.mean), ...
    'merged pooled tau differs from the serial run');
assert(isequal(serial.seeds(:), merged.seeds(:)), 'merged seeds differ');
fprintf('  3 chunks merged out of order == serial run, exactly: OK\n');

% --- 2. bookkeeping -----------------------------------------------------
assert(merged.meta.is_complete, 'a complete merge was not marked complete');
assert(isequal(merged.rep_index, 1:6), 'merged rep_index is %s', ...
       mat2str(merged.rep_index));
assert(numel(merged.meta.merged_from) == 3, 'merged_from does not list 3 chunks');
fprintf('  completeness, replication index and provenance recorded: OK\n');

% --- 3. refusals --------------------------------------------------------
threw = false;
try
    merge_montecarlo({ch{1}, ch{1}});
catch err
    threw = true;
    assert(~isempty(strfind(err.message, 'overlapping')), ...
        'wrong error for overlapping chunks: %s', err.message);
end
assert(threw, 'overlapping chunks were silently double-counted');

c_other = cfg;  c_other.seed = cfg.seed + 1;  c_other.mc.rep_range = [3 4];
other = run_montecarlo(c_other, 'sparse');
threw = false;
try
    merge_montecarlo({ch{1}, other});
catch
    threw = true;
end
assert(threw, 'chunks from different designs were merged');
fprintf('  overlapping chunks and mismatched designs both refused: OK\n');

% --- 4. incomplete merges are marked ------------------------------------
partial = merge_montecarlo(ch(1:2));
assert(~partial.meta.is_complete, 'an incomplete merge was marked complete');
assert(partial.meta.R_done == 4 && partial.meta.R == 6, ...
    'incomplete merge reports %d of %d', partial.meta.R_done, partial.meta.R);
fprintf('  incomplete merge kept, and marked incomplete (%d of %d): OK\n', ...
        partial.meta.R_done, partial.meta.R);

% --- 5. configuration round trip ----------------------------------------
raw = default_config();
assert(isa(raw.fmar.hyper_sd_rule, 'function_handle'), ...
    'the configuration no longer carries the hyperprior rule as a handle');
safe = cfg_to_savable(raw);
assert(ischar(safe.fmar.hyper_sd_rule), 'the handle was not stripped');
f = fullfile(tempdir, sprintf('cfg_roundtrip_%d.mat', round(rand * 1e6)));
save(f, '-struct', 'safe');                 % this is what used to fail
L = load(f);
delete(f);
back = cfg_from_saved(L);
assert(isa(back.fmar.hyper_sd_rule, 'function_handle'), ...
    'the handle was not restored');
h = 1:24;
d = max(abs(back.fmar.hyper_sd_rule(h) - raw.fmar.hyper_sd_rule(h)));
assert(d < 1e-15, 'the restored hyperprior rule differs by %.3g', d);
assert(strcmp(back.mode, raw.mode) && strcmp(back.blocks.scheme, raw.blocks.scheme), ...
    'ordinary string settings were mangled by the round trip');
fprintf('  cfg saves, loads and restores its function handle exactly: OK\n');

fprintf('PASS: test_chunk_merge\n\n');
end
