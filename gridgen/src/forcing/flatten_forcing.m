function specs = flatten_forcing(forcing)
%FLATTEN_FORCING  Normalise params.forcing to a cell of single-source specs.
%
%   specs = flatten_forcing(forcing)
%
%   Accepts a struct, a cell of structs, or an actuator struct with a vector xc
%   (and scalar or matching-length w) describing several actuators at once.
%   Returns a flat cell where each entry is one actuator or one distribution.

if ~iscell(forcing); forcing = {forcing}; end

specs = {};
for i = 1:numel(forcing)
    s = forcing{i};
    if isfield(s, 'xc') && numel(s.xc) > 1          % vector actuator -> split
        xc = s.xc(:).';
        w  = s.w;  if isscalar(w); w = w * ones(size(xc)); end
        assert(numel(w) == numel(xc), 'flatten_forcing: w must be scalar or match xc.');
        for k = 1:numel(xc)
            sk = s;  sk.xc = xc(k);  sk.w = w(k);
            specs{end+1} = sk;
        end
    else
        specs{end+1} = s;
    end
end
end
