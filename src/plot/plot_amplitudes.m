function plot_amplitudes(StabGrid, StabRes, ref_data)
%PLOT_AMPLITUDES  Final amplitude development plot.
%
%   plot_amplitudes(StabGrid, StabRes)
%   plot_amplitudes(StabGrid, StabRes, ref_data)
%
%   ref_data (optional) struct with:
%     .file   — path to .mat file with reference StabRes/StabGrid
%     .label  — legend prefix (e.g. 'NPSE', 'LPSE', 'DNS')
%
%   Mode labels use StabRes.omegavec/betavec to determine (m,n) indices,
%   so 3D cases (N>0) are labelled correctly.

figure(4); clf; set(gcf, 'Position', [100 100 900 500]); hold on;

marker = {'o','+','*','s','^','v','d','x','<','>','p','h'};
xmark  = linspace(StabGrid.xun(1), StabGrid.xun(end), 15);
h_list = [];
legs   = {};

% Derive fundamental frequency/wavenumber for integer (m,n) labels
omega_vec = StabRes.omegavec;
beta_vec  = StabRes.betavec;
omega_f = min(omega_vec(omega_vec > 0));  % fundamental frequency
beta_f  = min(beta_vec(beta_vec > 0));    % fundamental wavenumber

A_all = [];  % collect all plotted amplitudes for auto ylim

% HNS modes (black)
mk_count = 0;
for j = 1:size(StabRes.A, 1)
    if ~any(StabRes.A(j,:)); continue; end
    Aj = StabRes.A(j,:) / sqrt(2);
    A_all = [A_all, Aj(Aj > 0)]; %#ok<AGROW>

    if ~isempty(omega_f); m = round(omega_vec(j)/omega_f); else; m = 0; end
    if ~isempty(beta_f);  n = round(beta_vec(j)/beta_f);   else; n = 0; end
    lbl = mode_label('HNS', m, n);

    if m == 0 && n == 0
        h = semilogy(StabGrid.xun, Aj, 'k--', 'LineWidth', 1.5);
    else
        mk_count = mk_count + 1;
        semilogy(StabGrid.xun, Aj, 'k-', 'LineWidth', 1.5);
        ym = interp1(StabGrid.xun, Aj, xmark);
        mk = marker{mod(mk_count-1, numel(marker))+1};
        h = semilogy(xmark, ym, mk, 'Color', 'k', 'MarkerSize', 5);
    end
    h_list(end+1) = h; %#ok<AGROW>
    legs{end+1}   = lbl; %#ok<AGROW>
end

% Reference data (red, optional)
if nargin >= 3 && ~isempty(ref_data) && isfile(ref_data.file)
    S = load(ref_data.file);
    lbl_prefix = 'REF';
    if isfield(ref_data, 'label'); lbl_prefix = ref_data.label; end

    [ref_res, ref_grid] = find_stab_structs(S);

    if ~isempty(ref_res) && ~isempty(ref_grid)
        % Reference (m,n) labels from its own omegavec/betavec if available
        if isfield(ref_res,'omegavec') && isfield(ref_res,'betavec')
            ref_ov = ref_res.omegavec;
            ref_bv = ref_res.betavec;
            ref_of = min(ref_ov(ref_ov > 0));
            ref_bf = min(ref_bv(ref_bv > 0));
        else
            ref_ov = []; ref_bv = []; ref_of = []; ref_bf = [];
        end

        mk_count_r = 0;
        for j = 1:size(ref_res.A, 1)
            if ~any(ref_res.A(j,:)); continue; end
            Aj = ref_res.A(j,:) / sqrt(2);
            A_all = [A_all, Aj(Aj > 0)]; %#ok<AGROW>

            if ~isempty(ref_of) && numel(ref_ov) >= j
                m = round(ref_ov(j)/ref_of);
            else
                m = j - 1;
            end
            if ~isempty(ref_bf) && numel(ref_bv) >= j
                n = round(ref_bv(j)/ref_bf);
            else
                n = 0;
            end
            lbl = mode_label(lbl_prefix, m, n);

            if m == 0 && n == 0
                h = semilogy(ref_grid.xun, Aj, 'r--', 'LineWidth', 1.5);
            else
                mk_count_r = mk_count_r + 1;
                semilogy(ref_grid.xun, Aj, 'r-', 'LineWidth', 1.5);
                ym = interp1(ref_grid.xun, Aj, xmark);
                mk = marker{mod(mk_count_r-1, numel(marker))+1};
                h = semilogy(xmark, ym, mk, 'Color', 'r', 'MarkerSize', 5);
            end
            h_list(end+1) = h; %#ok<AGROW>
            legs{end+1}   = lbl; %#ok<AGROW>
        end
    else
        warning('Could not find StabRes/StabGrid in %s', ref_data.file);
    end
end

set(gca, 'YScale', 'log');
xlim(StabGrid.xun([1 end]))
if ~isempty(A_all)
    ymin = 10^floor(log10(min(A_all)));
    ymax = 10^ceil(log10(max(A_all)));
    ylim([max(ymin, 1e-12), min(ymax, 1)]);
end
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u''_{\mathrm{rms}}$', 'Interpreter', 'latex', 'FontSize', 12)
title('Amplitude development', 'FontSize', 13)
if ~isempty(h_list)
    legend(h_list, legs, 'Location', 'EastOutside', 'FontSize', 8)
end
grid on; box on;

end

%% Local helpers

function lbl = mode_label(prefix, m, n)
if n == 0
    lbl = sprintf('%s (%d,0)', prefix, m);
else
    lbl = sprintf('%s (%d,%d)', prefix, m, n);
end
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
