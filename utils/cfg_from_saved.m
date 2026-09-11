function c = cfg_from_saved(c)
% PURPOSE
% -------
% Inverse of utils/cfg_to_savable.m: turn the stored source strings of
% anonymous functions back into callable handles, so a configuration
% loaded from a .mat file can be handed straight back to the estimators.
%
% A char value is converted only when it looks like an anonymous
% function, i.e. its first non-blank character is '@'.  No other string
% in a configuration struct starts that way (modes, filenames, estimator
% keys and block schemes are all plain words), so ordinary settings are
% left untouched.
%
% INPUTS / OUTPUTS
% ----------------
% c : any struct / cell / value; returned with '@...' strings converted
%     to function handles, recursively.
%
% NOTES
% -----
% If str2func fails on a string, the string is left as it is rather than
% raising: a configuration that cannot be fully restored should still be
% readable.

if ischar(c) && ~isempty(c)
    t = strtrim(c);
    if ~isempty(t) && t(1) == '@'
        try
            c = str2func(t);
        catch
            % leave the string in place; see NOTES
        end
    end
    return
end
if isstruct(c)
    for e = 1:numel(c)
        f = fieldnames(c);
        for k = 1:numel(f)
            c(e).(f{k}) = cfg_from_saved(c(e).(f{k}));
        end
    end
    return
end
if iscell(c)
    for k = 1:numel(c)
        c{k} = cfg_from_saved(c{k});
    end
end
end
