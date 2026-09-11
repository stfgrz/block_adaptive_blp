function [ok, why] = can_plot()
% PURPOSE
% -------
% Report whether figures can be created in this session.
%
% WHY THIS EXISTS
% ---------------
% The demonstration scripts do two things: they compute numbers and they
% draw figures.  On a headless machine (a container, a compute node, a
% CI runner) Octave may have NO graphics toolkit installed at all, and
% `figure` then raises "no graphics toolkits are available!".  Without a
% guard that error kills the script AFTER the expensive estimation and
% BEFORE anything is printed or saved -- the worst possible failure
% mode.  The scripts therefore ask this first and skip the figures with
% a message, keeping every number they were run for.
%
% MATLAB always has graphics, so this returns true there.
%
% OUTPUTS
% -------
% ok  : logical, true when a figure can be created.
% why : '' when ok, otherwise a short explanation to print.

ok = true;  why = '';
if exist('OCTAVE_VERSION', 'builtin') ~= 5
    return                      % MATLAB: always has graphics
end
tk = available_graphics_toolkits();
if isempty(tk)
    ok = false;
    why = ['no graphics toolkit is available in this Octave build ' ...
           '(install gnuplot or the fltk/qt widgets to get figures); ' ...
           'the numerical output is unaffected'];
end
end
