function plot_convergence(StabGrid, Mode, Stab, Opt)
%PLOT_CONVERGENCE  Plot convergence history.

if isempty(findobj('Type','figure','Number',3))
    figure(3); set(gcf, 'Visible', 'on');
else
    set(0, 'CurrentFigure', 3); clf;
end
set(gca, 'YScale', 'log');
hold on;
colors = lines(numel(Mode.RunJ));
k = 0;
modeslegend = {};
iters = 1:Opt.k1;

for j = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2))
    k = k + 1;
    semilogy(iters, Opt.dalsave(j, iters), '-o', 'Color', colors(k,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4);
    modeslegend{end+1} = sprintf('(m,n) = (%d,%d)', Stab.Mvec(j), Stab.Nvec(j));
end

% Convergence thresholds
semilogy(iters, Opt.Conv*Opt.ConvF .* ones(1,Opt.k1), 'k--', 'LineWidth', 0.8);
semilogy(iters, Opt.Conv .* ones(1,Opt.k1), 'k-', 'LineWidth', 1.0);
modeslegend{end+1} = 'Intermediate threshold';
modeslegend{end+1} = 'Final threshold';

xlim([0 Opt.k1+1])
xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Relative error $\delta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Convergence history', 'Interpreter', 'latex', 'FontSize', 13)
legend(modeslegend, 'Location', 'best', 'FontSize', 9)
grid on; box on;
set(gcf, 'Position', [100 550 700 400]);
drawnow;

end
