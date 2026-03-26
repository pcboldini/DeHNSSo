function plot_amplitudes(StabGrid, StabRes, ref_data)
%PLOT_AMPLITUDES  Final amplitude development plot.
%
%   plot_amplitudes(StabGrid, StabRes)
%   plot_amplitudes(StabGrid, StabRes, ref_data)
%
%   ref_data (optional) struct with:
%     .file   — path to .mat file with reference StabRes/StabGrid
%     .label  — legend prefix (e.g. 'NPSE', 'LPSE', 'DNS')

figure(4); clf; set(gcf, 'Position', [100 100 900 500]); hold on;

marker = {'o','+','*','s','^','v','d','x','<','>','p','h'};
xmark  = linspace(StabGrid.xun(1), StabGrid.xun(end), 15);
h_list = [];
legs   = {};

% HNS modes (black)
for j = 1:size(StabRes.A, 1)
    if ~any(StabRes.A(j,:)); continue; end
    m_idx = j - 1;
    if m_idx == 0
        h = semilogy(StabGrid.xun, StabRes.A(j,:)/sqrt(2), 'k--', 'LineWidth', 1.5);
    else
        semilogy(StabGrid.xun, StabRes.A(j,:)/sqrt(2), 'k-', 'LineWidth', 1.5);
        ym = interp1(StabGrid.xun, StabRes.A(j,:)/sqrt(2), xmark);
        mk = marker{min(m_idx, numel(marker))};
        h = semilogy(xmark, ym, mk, 'Color', 'k', 'MarkerSize', 5);
    end
    h_list(end+1) = h;
    legs{end+1} = sprintf('HNS (%d,0)', m_idx);
end

% Reference data (red, optional)
if nargin >= 3 && ~isempty(ref_data) && isfile(ref_data.file)
    S = load(ref_data.file);
    lbl = 'REF';
    if isfield(ref_data, 'label'); lbl = ref_data.label; end

    % Find StabRes and StabGrid in loaded file
    [ref_res, ref_grid] = find_stab_structs(S);

    if ~isempty(ref_res) && ~isempty(ref_grid)
        for j = 1:size(ref_res.A, 1)
            if ~any(ref_res.A(j,:)); continue; end
            m_idx = j - 1;
            if m_idx == 0
                h = semilogy(ref_grid.xun, ref_res.A(j,:)/sqrt(2), 'r--', 'LineWidth', 1.5);
            else
                semilogy(ref_grid.xun, ref_res.A(j,:)/sqrt(2), 'r-', 'LineWidth', 1.5);
                ym = interp1(ref_grid.xun, ref_res.A(j,:)/sqrt(2), xmark);
                mk = marker{min(m_idx, numel(marker))};
                h = semilogy(xmark, ym, mk, 'Color', 'r', 'MarkerSize', 5);
            end
            h_list(end+1) = h;
            legs{end+1} = sprintf('%s (%d,0)', lbl, m_idx);
        end
    else
        warning('Could not find StabRes/StabGrid in %s', ref_data.file);
    end
end

set(gca, 'YScale', 'log');
xlim(StabGrid.xun([1 end]))
ylim([1e-10 1e-1])
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u''_{\mathrm{rms}}$', 'Interpreter', 'latex', 'FontSize', 12)
title('Amplitude development', 'FontSize', 13)
legend(h_list, legs, 'Location', 'EastOutside', 'FontSize', 8)
grid on; box on;

end

function [res, grid] = find_stab_structs(S)
%FIND_STAB_STRUCTS  Find StabRes and StabGrid in a loaded .mat struct.
    res = []; grid = [];
    fnames = fieldnames(S);
    for k = 1:numel(fnames)
        v = S.(fnames{k});
        if isstruct(v) && isfield(v, 'A') && isfield(v, 'alpha')
            res = v;
        end
        if isstruct(v) && isfield(v, 'xun')
            grid = v;
        end
    end
end
