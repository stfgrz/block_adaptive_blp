function [X, names] = ea_pre_event_info(Y, varnames)
% EA_PRE_EVENT_INFO  Public information available BEFORE the events of
% month t, from the dataset itself (Bauer & Swanson 2023 logic).
%
% Columns (all known at the end of month t-1; NaN where lags are missing):
%   stoxx_ret3    100 log change of the stock index, t-4 -> t-1
%   ip_gr12       100 log change of industrial production, t-13 -> t-1
%   infl12        100 log change of the price level, t-13 -> t-1 (or the
%                 level of a year-on-year inflation variable at t-1)
%   ind_chg3      change of the policy indicator (column 1), t-4 -> t-1
% Used (i) to orthogonalise an instrument on pre-event information
% ('bs' variants) and (ii) in the predictability test of ea_relevance_iv.
% Variables are located by name: 'ip', 'hicp*', 'stoxx'; the indicator is
% column 1.  Missing variables are skipped.

[T, K] = size(Y);  %#ok<ASGLU>
X = zeros(T, 0);  names = {};
lagdiff = @(v, a, b) [nan(b, 1); v(b + 1 - a:end - a) - v(1:end - b)];  % v(t-a) - v(t-b), b > a
j = find(strcmp(varnames, 'stoxx'), 1);
if ~isempty(j), X(:, end + 1) = lagdiff(Y(:, j), 1, 4);  names{end + 1} = 'stoxx_ret3'; end
j = find(strcmp(varnames, 'ip'), 1);
if ~isempty(j), X(:, end + 1) = lagdiff(Y(:, j), 1, 13);  names{end + 1} = 'ip_gr12'; end
j = find(strncmp(varnames, 'hicp', 4), 1);
if ~isempty(j)
    if strcmp(varnames{j}, 'hicp_yoy')
        X(:, end + 1) = [NaN; Y(1:end - 1, j)];
    else
        X(:, end + 1) = lagdiff(Y(:, j), 1, 13);
    end
    names{end + 1} = 'infl12';
end
X(:, end + 1) = lagdiff(Y(:, 1), 1, 4);  names{end + 1} = 'ind_chg3';
end
