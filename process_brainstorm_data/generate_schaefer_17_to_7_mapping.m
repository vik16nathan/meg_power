%Builds s600_17_to_7.mat (a region-name -> region-name dictionary, mapping
%each Schaefer-600 17-network scout to the SAME anatomical region's name
%under the Schaefer-600 7-network scheme) and s600_region_names.mat (the
%list of unique 7-network-scheme region names). This only needs to run once
%-- the mapping is a fixed correspondence between two atlas naming schemes,
%not something that varies per subject -- so it uses the FIRST subject that
%has both atlases registered (see cat_defaults.m: 'Schaefer2018_600P_7N' was
%added there so future batch_cat12.m runs compute it alongside the existing
%17-network atlas).
%
%Matches the logic of brainstorm3/brainpower/rest/convert_17_to_7_net.m
%(matching the two atlases' scouts by their first vertex, since both cover
%the exact same 600-region parcellation, just labeled under each scheme).
%
%NOT YET RUN/VERIFIED against a live Brainstorm session.

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';

s600_17_to_7 = containers.Map();
found = false;

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    surf_file = sprintf('/export02/data/vikramn/OmegaSubset/anat/%s/tess_cortex_central_low.mat', subject);
    if ~exist(surf_file, 'file')
        continue;
    end
    surface_file = load(surf_file);
    atlas_names = {surface_file.Atlas.Name};

    i17 = find(contains(lower(atlas_names), 'schaefer') & contains(atlas_names, '600') & ...
        contains(lower(atlas_names), '17'), 1);
    i7 = find(contains(lower(atlas_names), 'schaefer') & contains(atlas_names, '600') & ...
        contains(lower(atlas_names), '7') & ~contains(lower(atlas_names), '17'), 1);

    if isempty(i17) || isempty(i7)
        fprintf('%s: missing one or both Schaefer-600 atlases (has: %s) -- skipping.\n', ...
            subject, strjoin(atlas_names, ', '));
        continue;
    end

    fprintf('Using %s: "%s" (17-net) + "%s" (7-net)\n', subject, atlas_names{i17}, atlas_names{i7});

    s600_17 = surface_file.Atlas(i17);
    s600_7  = surface_file.Atlas(i7);
    num_scouts = numel(s600_17.Scouts);
    if numel(s600_7.Scouts) ~= num_scouts
        error('%s: 17-net (%d scouts) and 7-net (%d scouts) atlases have different scout counts.', ...
            subject, num_scouts, numel(s600_7.Scouts));
    end

    % Build region -> first-vertex dictionaries for each atlas, then sort
    % both by vertex so matching rows correspond to the same anatomical region.
    label_v_dict_7 = cell(num_scouts, 2);
    label_v_dict_17 = cell(num_scouts, 2);
    for i = 1:num_scouts
        label_v_dict_7(i,:)  = {s600_7.Scouts(i).Label,  s600_7.Scouts(i).Vertices(1)};
        label_v_dict_17(i,:) = {s600_17.Scouts(i).Label, s600_17.Scouts(i).Vertices(1)};
    end
    [~, sort_ind_7]  = sort(cell2mat(label_v_dict_7(:,2)));
    [~, sort_ind_17] = sort(cell2mat(label_v_dict_17(:,2)));
    label_v_dict_7_sort  = label_v_dict_7(sort_ind_7,:);
    label_v_dict_17_sort = label_v_dict_17(sort_ind_17,:);

    for i = 1:num_scouts
        % containers.Map keys/values must be char vectors, not MATLAB
        % string objects -- values(s600_17_to_7) later feeds unique(),
        % which rejects a cell array of string objects with "Cell array
        % input must be a cell array of character vectors."
        name_17 = char(label_v_dict_17_sort{i,1});
        name_7  = char(label_v_dict_7_sort{i,1});
        s600_17_to_7(name_17) = name_7; %#ok<SAGROW>
    end

    found = true;
    break;
end

if ~found
    error(['No subject with both Schaefer-600 17-network and 7-network atlases found yet. ' ...
        'Schaefer2018_600P_7N was added to cat_defaults.m, so subjects processed by ' ...
        'batch_cat12.m from now on should have both -- wait for at least one more subject ' ...
        'to finish, then re-run this script.']);
end

save(fullfile(OUTPUTS_DIR, 's600_17_to_7.mat'), 's600_17_to_7');

% region_names: the unique set of 7-network-scheme region names that
% parcellate_cortical_column_res_vol.m averages across subjects (excludes the medial
% wall, which is filtered out separately by s17_names_to_ignore there, but
% harmless to include here since it's just the list of dictionary keys used).
region_names = unique(values(s600_17_to_7));
save(fullfile(OUTPUTS_DIR, 's600_region_names.mat'), 'region_names');

fprintf('Saved s600_17_to_7.mat (%d entries) and s600_region_names.mat (%d unique regions).\n', ...
    s600_17_to_7.Count, numel(region_names));
