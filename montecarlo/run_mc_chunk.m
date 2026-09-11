function out_file = run_mc_chunk(preset, dgp_name, r0, r1, out_file, overrides)
% PURPOSE
% -------
% Run ONE chunk of one Monte Carlo experiment and save it, so a single
% design can be spread over several processes on a multi-core machine.
% Replication seeds are pure functions of the replication index, so the
% chunks reassemble (montecarlo/merge_montecarlo.m) into exactly the
% result a serial run would have produced.
%
% USAGE (one process per chunk; the shell driver launches four at a time)
% -----------------------------------------------------------------------
%   octave --eval "addpath(genpath('<repo>')); \
%       run_mc_chunk('final','sparse',1,125,'<repo>/results/chunks/final_sparse_1.mat')"
%
% INPUTS
% ------
% preset    : name passed to montecarlo/mc_preset.m.
% dgp_name  : 'correct' | 'sparse' | 'intermediate' | 'dense'.
% r0, r1    : first and last replication index of this chunk.
% out_file  : where to save (directories are created).
% overrides : OPTIONAL struct of cfg overrides (see mc_preset).
%
% OUTPUTS
% -------
% out_file : the path written.  The file contains `mc` and `cfg`.
%
% NOTES
% -----
% Prints a one-line header naming the preset, DGP and chunk, so a log of
% several parallel processes stays readable.

if nargin < 6, overrides = struct(); end
cfg = mc_preset(preset, overrides);
cfg.mc.rep_range = [r0, r1];

fprintf('[run_mc_chunk] preset=%s dgp=%s reps=%d-%d T=%d p=%d H=%d R=%d\n', ...
        preset, dgp_name, r0, r1, cfg.T, cfg.p, cfg.H, cfg.mc.n_rep);

mc = run_montecarlo(cfg, dgp_name);
mc.meta.preset = preset;
mc.meta.overrides = overrides;

d = fileparts(out_file);
if ~isempty(d) && exist(d, 'dir') ~= 7
    [ok, msg] = mkdir(d);
    assert(ok == 1, 'run_mc_chunk: cannot create %s (%s)', d, msg);
end
cfg = cfg_to_savable(cfg);      % handles cannot be written to a .mat
save(out_file, 'mc', 'cfg', '-v7');
fprintf('[run_mc_chunk] wrote %s\n', out_file);
end
