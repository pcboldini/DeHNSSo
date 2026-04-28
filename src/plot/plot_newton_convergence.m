function plot_newton_convergence(Opt, Stab, Mode)
%PLOT_NEWTON_CONVERGENCE  Cumulative per-mode ||F|| with AF annotations.

if isempty(findobj('Type','figure','Number',4))
    figure(4); set(gcf, 'Visible', 'on');
else
    set(0, 'CurrentFigure', 4); clf;
end

F_hist  = Opt.newton_F_history;   % nphys x total_iters
AF_hist = Opt.newton_AF_history;  % 1 x total_iters
[nphys, n_pts] = size(F_hist);
iters = 1:n_pts;

half = round(Mode.nf/2);
phys_modes = Mode.RunJ(Mode.RunJ >= half);

set(gca, 'YScale', 'log');
hold on;

% Per-mode ||F|| lines
all_phys = half:Mode.nf;
colors = lines(numel(all_phys));
modeslegend = {};
for kk = 1:min(nphys, size(F_hist,1))
    j = phys_modes(kk);
    cidx = j - half + 1;
    semilogy(iters, F_hist(kk,:), '-o', 'Color', colors(cidx,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4);
    modeslegend{end+1} = sprintf('(m,n) = (%d,%d)', Stab.Mvec(j), Stab.Nvec(j));
end

% Tolerance line
tol = 1e-10;
if isfield(Opt, 'newton_tol'); tol = Opt.newton_tol; end
semilogy(iters, tol .* ones(1, n_pts), 'k-', 'LineWidth', 1.0);
modeslegend{end+1} = 'Newton tolerance';

% AF change markers — vertical lines where AF changes
AF_changes = find(diff(AF_hist) ~= 0);
yl = ylim;
for ii = 1:numel(AF_changes)
    xline(AF_changes(ii) + 0.5, '--', sprintf('AF=%.3f', AF_hist(AF_changes(ii)+1)), ...
        'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'center', 'FontSize', 8);
end
% Label the first AF
if n_pts > 0
    text(1, yl(2)*0.5, sprintf('AF=%.3f', AF_hist(1)), ...
        'FontSize', 8, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'left');
end

xlim([0.5 n_pts+0.5])
xlabel('Newton sub-iteration', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\|F_k\|$', 'Interpreter', 'latex', 'FontSize', 12)
title('Newton convergence history', 'Interpreter', 'latex', 'FontSize', 13)
legend(modeslegend, 'Location', 'best', 'FontSize', 9)
grid on; box on;
set(gcf, 'Position', [810 550 700 400]);
drawnow;

end
