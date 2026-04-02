function [BF] = load_input(input)
%LOAD_INPUT  Load base-flow data from file.
%
%   BF = load_input(input)
%
%   Reads base-flow data from .mat, .csv, or external CFD output.
%   The input struct must contain:
%     input.folder   — directory containing the data file
%     input.format   — file format: 'mat', 'csv'
%     input.filename — (optional) specific file; auto-detects if empty
%
%   Output: BF struct with fields X, Y, U, V, W (and derivatives if available).
%   Field names from OpenFOAM-style CSV ('x','y','U:0') are auto-renamed.
if isfield(input, 'filename') && ~isempty(input.filename)
    filepath = fullfile(input.folder, input.filename);
else
    ext = ['*.' input.format];
    files = dir(fullfile(input.folder, ext));
    if isempty(files)
        error('No %s files found in %s', input.format, input.folder);
    end
    filepath = fullfile(input.folder, files(1).name);
end

switch lower(input.format)
    case 'mat'
        load(filepath, 'BF');   % assumes BF is stored in the .mat file
    case 'bin'
        error('load_input:unsupported', 'Binary format not yet implemented.');
    case 'csv'
        warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
        BF = table2struct(readtable(filepath), 'ToScalar', true);
        warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
    otherwise
        error("Unknown input.format = '%s'. Use 'mat' or 'bin'.", input.format);
end

%% For 3D OpenFOAM exports (thin slab), keep only one z-slice
z_field = '';
if isfield(BF, 'Points_2'),  z_field = 'Points_2';
elseif isfield(BF, 'z'),     z_field = 'z';
end
if ~isempty(z_field)
    z_vals  = BF.(z_field);
    if isfield(input, 'z_slice')
        z_target = input.z_slice;
    else
        z_target = min(z_vals);
    end
    z_tol = 1e-6 * max(abs(z_vals(:)));  % relative tolerance for z-slice
    if z_tol == 0; z_tol = 1e-10; end
    keep = abs(z_vals - z_target) < z_tol;
    fnames = fieldnames(BF);
    for i = 1:numel(fnames)
        BF.(fnames{i}) = BF.(fnames{i})(keep);
    end
    BF = rmfield(BF, z_field);
    fprintf('Filtered z-slice: kept %d of %d points\n', sum(keep), numel(z_vals));
end

%% Rename fields to standard names (handles OpenFOAM and generic CSV formats)
rename_map = {'x','X'; 'y','Y'; 'u','U'; 'v','V'; 'w','W'; ...
              'Points_0','X'; 'Points_1','Y'; ...
              'U_0','U'; 'U_1','V'; 'U_2','W'};
for i = 1:size(rename_map, 1)
    lo = rename_map{i,1};  up = rename_map{i,2};
    if isfield(BF, lo) && ~isfield(BF, up)
        BF.(up) = BF.(lo);
        BF = rmfield(BF, lo);
    end
end

%% Set W to zero if not present (2D problems)
if ~isfield(BF, 'W')
    BF.W = zeros(size(BF.X));
end

%% Validate BF structure
required_fields = {'X', 'Y', 'U', 'V'};
missing_fields = {};
for i = 1:length(required_fields)
    if ~isfield(BF, required_fields{i})
        missing_fields{end+1} = required_fields{i};
    end
end

if ~isempty(missing_fields)
    error('BF:missingFields', 'BF is missing required fields: %s', strjoin(missing_fields, ', '));
end

fprintf('BF loaded\n');

end

