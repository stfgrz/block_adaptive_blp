function fig = plot_irfs(theta_true, est_list, cfg, ttl)
% PURPOSE
% -------
% Compare estimated impulse responses against the truth: one subplot
% per response variable, all estimators overlaid as lines, with an
% uncertainty band drawn for any estimator that supplies lo/hi bounds
% and has .show_band = true (by default only the block-adaptive BLP,
% to keep the figure readable).
%
% MODEL / EQUATIONS
% -----------------
% Everything plotted is the response of variable i at horizons
% h = 0..H to the recursive unit shock to cfg.shock_var.
%
% INPUTS
% ------
% theta_true : (K x (H+1)) true IRF, or [] to omit the truth line.
% est_list   : cell array of structs, each with fields
%                .name  : legend label (char)
%                .theta : (K x (H+1)) point estimates
%                .lo/.hi: (K x (H+1)) bounds or [] (optional)
%                .show_band : logical (optional, default false)
% cfg        : configuration struct (cfg.H, cfg.shock_var).
% ttl        : overall figure title (char).
%
% OUTPUTS
% -------
% fig : figure handle.
%
% DIMENSIONS
% ----------
% K subplots in one row (K is small in the prototype).
%
% NOTES
% -----
% Colours come from a small editable list below; nothing else about
% appearance is hard-coded.  Bands are drawn as shaded patches without
% transparency so the figure renders identically across MATLAB
% versions and graphics back-ends.

K = size(est_list{1}.theta, 1);
H = cfg.H;
hgrid = 0:H;

% Editable style choices.
colors = [0.85 0.33 0.10;   % 1 orange-red
          0.00 0.45 0.74;   % 2 blue
          0.47 0.67 0.19;   % 3 green
          0.49 0.18 0.56;   % 4 purple
          0.30 0.75 0.93];  % 5 cyan
band_fade = 0.80;           % 0 = full colour, 1 = white

fig = figure('Name', ttl, 'Color', 'w', 'Position', [80 80 380*K 340]);
for i = 1:K
    subplot(1, K, i); hold on;
    pos = get(gca, 'Position'); pos(4) = pos(4) * 0.88;
    set(gca, 'Position', pos);          % headroom for the overall title

    % Bands first so lines stay on top.
    for e = 1:numel(est_list)
        est = est_list{e};
        if isfield(est, 'show_band') && est.show_band && ...
           isfield(est, 'lo') && ~isempty(est.lo)
            ok = isfinite(est.lo(i, :)) & isfinite(est.hi(i, :));
            hh = hgrid(ok);
            cb = colors(e, :) * (1 - band_fade) + band_fade;
            fill([hh, fliplr(hh)], ...
                 [est.lo(i, ok), fliplr(est.hi(i, ok))], cb, ...
                 'EdgeColor', 'none');
        end
    end
    hnd = [];  labels = {};
    if ~isempty(theta_true)
        hnd(end + 1) = plot(hgrid, theta_true(i, :), 'k-', 'LineWidth', 2.2);
        labels{end + 1} = 'true';
    end
    for e = 1:numel(est_list)
        hnd(end + 1) = plot(hgrid, est_list{e}.theta(i, :), '-', ...
             'Color', colors(e, :), 'LineWidth', 1.4);
        labels{end + 1} = est_list{e}.name;
    end
    plot(hgrid, zeros(size(hgrid)), 'k:', 'LineWidth', 0.5);

    xlabel('horizon h');
    title(sprintf('response of y_%d', i));
    if i == 1
        ylabel(sprintf('IRF to unit shock in y_%d', cfg.shock_var));
        legend(hnd, labels, 'Location', 'northeast', 'Box', 'off');
    end
    xlim([0 H]); grid on; box on;
end
% Overall title: text on an invisible full-figure axes (renders centred
% in both MATLAB and Octave, unlike annotation under some toolkits).
axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.97, ttl, 'HorizontalAlignment', 'center', ...
     'FontWeight', 'bold', 'Interpreter', 'none');
end
