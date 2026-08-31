function fig = plot_rmse(s, cfg)
% PURPOSE
% -------
% Plot Monte Carlo RMSE by horizon, one subplot per response variable,
% one line per estimator.  Lower is better; the interesting question is
% WHERE (which horizons / variables) each estimator wins.
%
% MODEL / EQUATIONS
% -----------------
% rmse(e,i,h) as defined in summarize_montecarlo.m.
%
% INPUTS
% ------
% s   : summary struct from summarize_montecarlo.
% cfg : configuration struct (cfg.H used for the axis).
%
% OUTPUTS
% -------
% fig : figure handle.
%
% DIMENSIONS
% ----------
% s.rmse is (4 x K x H); horizons 1..H are plotted (h = 0 excluded by
% construction of the metrics).
%
% NOTES
% -----
% Colours in one editable list; solid lines throughout.

[nE, K, H] = size(s.rmse);
colors = [0.85 0.33 0.10;
          0.00 0.45 0.74;
          0.47 0.67 0.19;
          0.49 0.18 0.56];

fig = figure('Name', sprintf('RMSE (%s DGP)', s.dgp_name), ...
             'Color', 'w', 'Position', [100 100 380*K 320]);
for i = 1:K
    subplot(1, K, i); hold on;
    pos = get(gca, 'Position'); pos(4) = pos(4) * 0.88;
    set(gca, 'Position', pos);
    hnd = [];
    for e = 1:nE
        hnd(end + 1) = plot(1:H, squeeze(s.rmse(e, i, :)), '-', ...
             'Color', colors(e, :), 'LineWidth', 1.5);
    end
    xlabel('horizon h'); title(sprintf('RMSE, response of y_%d', i));
    if i == 1
        ylabel('RMSE across replications');
        legend(hnd, s.est_names, 'Location', 'northwest', 'Box', 'off');
    end
    xlim([1 H]); grid on; box on;
end
axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.97, sprintf('Monte Carlo RMSE by horizon -- DGP: %s (R = %d)', ...
     s.dgp_name, s.R), 'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold', 'Interpreter', 'none');
end
