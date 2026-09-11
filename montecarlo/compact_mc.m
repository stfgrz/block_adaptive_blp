function mc = compact_mc(mc)
% PURPOSE
% -------
% Shrink a MERGED Monte Carlo result before it is written to disk, by
% collapsing the per-replication posterior-quantile arrays of tau to
% their replication averages.
%
% WHY THIS IS LOSSLESS FOR EVERY CONSUMER
% ---------------------------------------
% mc.tau.<est>.q is stored per replication (K x G x H x nQ x R) for ONE
% reason: so that chunks of a parallel run merge by plain concatenation
% along the replication dimension.  Once merged, a result is never
% re-merged, and every consumer of that field averages it over
% replications -- tau_diagnostics computes mean(q, 5), which is the
% identity on an already-averaged 4-d array, so nothing downstream can
% tell the difference.  At R = 500 the array is a few megabytes per
% estimator; collapsing it keeps a stored headline result small enough
% to live in the repository next to the earlier baselines.
%
% Everything else -- every estimate, interval bound, tau mean/median,
% P(tau > 1), seed and diagnostic -- is left exactly as it is.
%
% INPUTS / OUTPUTS
% ----------------
% mc : a run_montecarlo / merge_montecarlo output, returned compacted.

if ~isfield(mc, 'tau') || ~isstruct(mc.tau), return; end
f = fieldnames(mc.tau);
for k = 1:numel(f)
    if isfield(mc.tau.(f{k}), 'q') && ~isempty(mc.tau.(f{k}).q)
        mc.tau.(f{k}).q = mean(mc.tau.(f{k}).q, 5);
    end
end
mc.meta.tau_quantiles_compacted = true;
end
