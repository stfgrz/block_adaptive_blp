function fig = plot_block_scales(tau, cfg, ttl, misspec_block)
% PURPOSE
% -------
% Visualise the posterior local scales tau_{g,h} of the block-adaptive
% BLP: one subplot per LP EQUATION (response variable), one line per
% BLOCK (regressor variable), horizons on the x-axis.  A horizontal
% reference line at tau = 1 marks "no escape from the VAR prior"; the
% truly misspecified block (if any) is highlighted in the line width.
%
% MODEL / EQUATIONS
% -----------------
% tau_{g,h} is the block scale in
%   beta_{g,h} - mu_{g,h}^VAR ~ N(0, lambda_h^2 tau_{g,h}^2 D_g);
% plotted values are posterior means (single dataset) or Monte Carlo
% averages of posterior means (from summarize_montecarlo.tau_bar).
%
% INPUTS
% ------
% tau           : (K x G x H) array of scales (equation x block x horizon).
% cfg           : configuration struct (cfg.H).
% ttl           : figure title (char).
% misspec_block : [] or integer index of the truly misspecified block
%                 (drawn thicker; pass [] if unknown / none).
%
% OUTPUTS
% -------
% fig : figure handle.
%
% DIMENSIONS
% ----------
% K subplots, G lines each.
%
% NOTES
% -----
% Log y-scale: the horseshoe scale is multiplicative by nature.

if nargin < 4
    misspec_block = [];
end
[K, G, H] = size(tau);
colors = [0.00 0.45 0.74;
          0.47 0.67 0.19;
          0.85 0.33 0.10;
          0.49 0.18 0.56];

fig = figure('Name', ttl, 'Color', 'w', 'Position', [120 120 380*K 320]);
for i = 1:K
    subplot(1, K, i); hold on;
    pos = get(gca, 'Position'); pos(4) = pos(4) * 0.88;
    set(gca, 'Position', pos);
    hnd = [];
    for g = 1:G
        lw = 1.2;
        if ~isempty(misspec_block) && g == misspec_block
            lw = 2.4;                   % highlight the true culprit
        end
        hnd(end + 1) = plot(1:H, squeeze(tau(i, g, :)), '-', ...
             'Color', colors(1 + mod(g - 1, size(colors, 1)), :), ...
             'LineWidth', lw);
    end
    plot([1 H], [1 1], 'k:', 'LineWidth', 0.8);   % tau = 1 reference
    set(gca, 'YScale', 'log');
    xlabel('horizon h');
    title(sprintf('equation for y_%d', i));
    if i == 1
        ylabel('posterior block scale \tau_{g,h}');
        labels = cell(1, G);
        for g = 1:G
            labels{g} = sprintf('block %d (lags of y_%d)', g, g);
        end
        legend(hnd, labels, 'Location', 'northwest', 'Box', 'off');
    end
    xlim([1 H]); grid on; box on;
end
axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.97, ttl, 'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold', 'Interpreter', 'none');
end
