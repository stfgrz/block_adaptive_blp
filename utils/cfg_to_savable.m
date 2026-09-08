function c = cfg_to_savable(c)
% PURPOSE
% -------
% Make a configuration struct SAFE TO SAVE by replacing every function
% handle it contains with the equivalent source string (func2str).
%
% WHY THIS EXISTS
% ---------------
% cfg.fmar.hyper_sd_rule is an anonymous function.  Octave's `save`
% refuses function handles outright ("wrong type argument 'function
% handle'") and MATLAB stores them in a form that silently loses the
% closure's workspace.  Either way, a configuration snapshot containing
% a handle is a reproducibility hazard: the very object that documents
% how a result was produced is the one part of it that does not survive
% being written to disk.  Everything that saves a cfg (run_montecarlo's
% mc.cfg, run_mc_chunk, run_experiment_grid) passes it through here
% first, and utils/cfg_from_saved.m converts the strings back.
%
% The conversion is recursive through nested structs and cell arrays, so
% it does not need to know where the handles are.
%
% INPUTS
% ------
% c : any struct / cell / value.
%
% OUTPUTS
% -------
% c : the same value with every function handle replaced by func2str(.)
%     of it (a char row vector beginning with '@').
%
% NOTES
% -----
% Round trip: cfg_from_saved(cfg_to_savable(cfg)) restores the handles.
% Anonymous functions that capture variables are reconstructed from
% their source text, so a handle that closed over a workspace variable
% would come back bound to whatever that name means at restore time;
% none of this project's handles capture anything (they are pure
% functions of their argument), and that is asserted in
% tests/test_mc_metrics.m by checking the restored rule reproduces the
% stored hyper_sd values.

if isa(c, 'function_handle')
    c = func2str(c);
    return
end
if isstruct(c)
    for e = 1:numel(c)
        f = fieldnames(c);
        for k = 1:numel(f)
            c(e).(f{k}) = cfg_to_savable(c(e).(f{k}));
        end
    end
    return
end
if iscell(c)
    for k = 1:numel(c)
        c{k} = cfg_to_savable(c{k});
    end
end
end
