function names = bf_deriv_fields()
%BF_DERIV_FIELDS  Physical base-flow derivative fields that may be supplied on
%   the raw BF (analytic base flow) instead of being computed in compute_bf_derivatives.
names = {'dxU','dyU','dxxU','ddyU','dxyU', 'dxV','dyV','dxxV','ddyV','dxyV', ...
         'dxW','dyW','dxxW','ddyW', 'dxrho','dyrho', 'dxT','dyT','dxxT','ddyT'};
end
