function BF_new = resample_cartesian(BF, params)
%RESAMPLE_CARTESIAN  Interpolate a structured BF onto a Malik-Chebyshev η × equidistant ξ grid.
%
%   BF_new = resample_cartesian(BF, params)
%
%   Designed for BF data that already lives on a structured (nx_raw × ny_raw)
%   grid (e.g. heritage BF_SweptWing_flat.mat) but at a different resolution
%   than what DeHNSSo expects. resample_curvilinear / resample_unstructured
%   handle the body-fitted equivalents.
%
%   Inputs:
%     BF     — X, Y (1-D axes), U, V, W on a regular grid; or a flat scatter
%              (paired vectors), canonicalised on entry.
%     params — n_eta_new / n_xi_new   target eta / xi counts (default: source)
%              H, y_i                 wall-normal height / Malik height (H=0.99 span, y_i=H/10)
%              xi_range / xi_trim_*   restrict or trim the x-span
%              mug, sig, ag           streamwise clustering (see cluster_xi)
%              w_act                  actuator half-width, sizes sig (see prepare_forcing)
%              interp_method          interp2 method (default 'makima')
%              rescale, Uref, lref    non-dimensionalise U/X/Y before resampling
%   Output:
%     BF_new — DeHNSSo-ready (X, Y, U, V, W, H, y_i, cluster) + any
%              Re/lref/Uref/nu present in BF.

%% Canonicalise a flat-scatter input to gridded (ny x nx); no-op if already 2-D
BF = reshape_scatter(BF);

x_raw  = BF.X(:).';          % [1 × nx_raw]
y_raw  = BF.Y(:);            % [ny_raw × 1]
nx_raw = numel(x_raw);
ny_raw = numel(y_raw);
U_raw  = BF.U;   V_raw = BF.V;   W_raw = BF.W;
% Compressible fields (auto-detected by presence of ρ and T).
is_compr = isfield(BF,'rho') && isfield(BF,'T');
if is_compr
    rho_raw = BF.rho;
    T_raw   = BF.T;
end

%% Optional non-dim rescale: U/Uref, X/lref, Y/lref
do_rescale = isfield(params, 'rescale') && params.rescale;
if do_rescale
    assert(isfield(params,'Uref') && isfield(params,'lref'), ...
        'resample_cartesian: rescale=true requires params.Uref and params.lref.');
    U_raw = U_raw / params.Uref;
    V_raw = V_raw / params.Uref;
    W_raw = W_raw / params.Uref;
    x_raw = x_raw / params.lref;
    y_raw = y_raw / params.lref;
    fprintf('Rescaled: U/=%.3e (Uref),  X,Y/=%.3e (lref)\n', params.Uref, params.lref);
end

%% Target grid dimensions — default to source resolution
if isfield(params, 'n_eta_new') && ~isempty(params.n_eta_new)
    ny_new = params.n_eta_new;
else
    ny_new = ny_raw;
end

if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    nx_new = params.n_xi_new;
else
    nx_new = nx_raw;
end

%% Wall-normal domain (H, y_i)
if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
elseif isfield(BF, 'H') && ~isempty(BF.H)
    H = BF.H;
else
    H = abs(y_raw(end) - y_raw(1)) * 0.99;   % BF.Y may be descending (η: top→wall)
end
if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
elseif isfield(BF, 'y_i') && ~isempty(BF.y_i)
    y_i = BF.y_i;
else
    y_i = H / 10;
end
if y_i >= H/2
    y_i = H/3;
    warning('y_i clamped to H/3 = %.5f (must be < H/2 for Malik mapping).', y_i);
end

%% Malik-Chebyshev η — descending: η(1)=H (freestream), η(end)=0 (wall)
a = H * y_i / (H - 2*y_i);
b = 1 + 2*a/H;
k_vec   = (0 : ny_new - 1)';
z_cheb  = cos(pi * k_vec / (ny_new - 1));
eta_new = a * (1 + z_cheb) ./ (b - z_cheb);             % [ny_new × 1], η in [0, H]

%% Streamwise ξ — optionally trimmed, optionally clustered
% Priority: xi_range overrides xi_trim_*. xi_trim_* drops a fraction from
% each end of the x-span (mirrors curvilinear convention).
if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    x_lo = x_raw(1) + params.xi_range(1) * (x_raw(end) - x_raw(1));
    x_hi = x_raw(1) + params.xi_range(2) * (x_raw(end) - x_raw(1));
else
    trim_lo = 0;  trim_hi = 0;
    if isfield(params, 'xi_trim_inflow')  && ~isempty(params.xi_trim_inflow)
        trim_lo = params.xi_trim_inflow;
    end
    if isfield(params, 'xi_trim_outflow') && ~isempty(params.xi_trim_outflow)
        trim_hi = params.xi_trim_outflow;
    end
    x_lo = x_raw(1) + trim_lo       * (x_raw(end) - x_raw(1));
    x_hi = x_raw(end) - trim_hi     * (x_raw(end) - x_raw(1));
end

% xi is x here, so the range ends map directly onto the clustering axis.
% w_act (actuator half-width) sizes sig to cover ~3x the orifice.
mug   = get_opt(params, 'mug', []);
ag    = get_opt(params, 'ag',  1);
w_act = get_opt(params, 'w_act', []);
if ~isempty(w_act)
    sig = 3 * w_act ./ (0.5 * (x_hi - x_lo));
else
    sig = get_opt(params, 'sig', 0.2);
end
xi_new = cluster_xi(x_lo, x_hi, nx_new, mug, sig, ag, x_lo, x_hi);   % [1 × nx_new]
report_clustering(mug, sig, ag);

%% Interpolate onto the target grid ('makima': no overshoot on non-uniform data)
interp_method = 'makima';
if isfield(params, 'interp_method') && ~isempty(params.interp_method)
    interp_method = params.interp_method;
end
[Xq, Yq] = meshgrid(xi_new, eta_new);
U_new = interp2(x_raw, y_raw, U_raw, Xq, Yq, interp_method);
V_new = interp2(x_raw, y_raw, V_raw, Xq, Yq, interp_method);
W_new = interp2(x_raw, y_raw, W_raw, Xq, Yq, interp_method);
% ρ, T (already non-dim by the BL caller; not affected by params.rescale)
if is_compr
    rho_new = interp2(x_raw, y_raw, rho_raw, Xq, Yq, interp_method);
    T_new   = interp2(x_raw, y_raw, T_raw,   Xq, Yq, interp_method);
end

%% Package into DeHNSSo-ready BF (no mesh_type -> build_stab_grid Cartesian branch)
BF_new       = struct();
BF_new.X     = xi_new;
BF_new.Y     = eta_new;
BF_new.U     = U_new;
BF_new.V     = V_new;
BF_new.W     = W_new;
BF_new.H     = H;
BF_new.y_i   = y_i;
BF_new.xi1D  = xi_new;
BF_new.eta1D = eta_new;

% Clustering metadata for the caller (xi = x, so xc_xi = mug); w is the
% provisioned actuator half-width. See build_stab_grid / add_wall_forcing.
BF_new.cluster = struct('mug', mug, 'sig', sig, 'ag', ag, 'xc_xi', mug, 'w', w_act);

% Compressible fields + thermo metadata (build_stab_grid runs Sutherland).
if is_compr
    BF_new.rho = rho_new;
    BF_new.T   = T_new;
    for f_name = {'Ma','Pr','gamma','Sstar'}
        if isfield(BF, f_name{1}) && ~isempty(BF.(f_name{1}))
            BF_new.(f_name{1}) = BF.(f_name{1});
        end
    end
end

for dn = bf_deriv_fields()
    if isfield(BF, dn{1}) && ~isempty(BF.(dn{1}))
        BF_new.(dn{1}) = interp2(x_raw, y_raw, BF.(dn{1}), Xq, Yq, interp_method);
    end
end

fprintf('Cartesian resample done: %d × %d → %d × %d (Malik-Cheb η, %s interp)\n', ...
        ny_raw, nx_raw, ny_new, nx_new, interp_method);
end


function BF = reshape_scatter(BF)
%RESHAPE_SCATTER  Flat-lattice scatter → canonical (ny × nx) gridded form.
%   X, Y are paired column vectors lying on a rectangular lattice; this
%   helper turns them into 1-D axes (X = 1×nx, Y = ny×1) and reshapes any
%   per-point numeric field of length numel(X) into a (ny × nx) array.
%   No-op when BF is already 2-D or when X/Y lengths differ.
if ~isvector(BF.X) || ~isvector(BF.Y) || numel(BF.X) ~= numel(BF.Y)
    return;
end

tolX = max(1e-6 * range(BF.X), 1e-12);
tolY = max(1e-6 * range(BF.Y), 1e-12);
[Xu, ~, ix] = uniquetol(BF.X(:), tolX, 'DataScale', 1);
[Yu, ~, iy] = uniquetol(BF.Y(:), tolY, 'DataScale', 1);
nx = numel(Xu);
ny = numel(Yu);

fprintf('Reshaping flat scatter → (ny=%d × nx=%d)\n', ny, nx);

fnames = setdiff(fieldnames(BF), {'X','Y'});
for k = 1:numel(fnames)
    f   = fnames{k};
    val = BF.(f);
    if ~isnumeric(val) || numel(val) ~= numel(ix)
        continue
    end
    M   = nan(ny, nx);
    lin = sub2ind([ny, nx], iy, ix);
    M(lin) = val(:);
    if any(isnan(M(:)))
        warning('resample_cartesian:gaps', ...
            'Field %s has %d missing lattice points after reshape.', f, nnz(isnan(M)));
    end
    BF.(f) = M;
end

BF.X = Xu(:).';   % 1 × nx
BF.Y = Yu(:);     % ny × 1
end
