%Map of L/A (resistance) and L*A (volume) across the brains of the resting
%individuals

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
output_dir = '/export02/data/vikramn/hbm_manuscript_code/outputs/resistances/';
load(fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat'));
avg_gm_res = 3.5; %ohm*meters
% sub-0011 excluded: outlier.
SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

% freq_band_names = {'delta', 'theta', 'alpha', ...
%     'beta', 'gamma1', 'gamma2'};
%atlas='s600';
atlas='s600';
load(fullfile(OUTPUTS_DIR, sprintf('%s_17_to_7.mat', atlas)));
load(fullfile(OUTPUTS_DIR, sprintf('%s_region_names.mat', atlas)));
band=''; %default - when working with non-band-pass filtered data
%also remove the for loop below
dk_index = 2;
% schaeffer_index was hardcoded to 6, but that's actually
% "Schaefer_400_17net" (402 scouts), not "Schaefer_600_17net" (602
% scouts, matching s600_region_names.mat's 602 entries) -- confirmed by
% listing each subject's registered atlases directly. Atlas registration
% order isn't even stable across subjects (sub-0035 has an extra atlas
% inserted, shifting later indices), so a fixed position is fragile
% regardless of which number is "correct" today. Found by NAME per subject
% instead, matching generate_schaefer_17_to_7_mapping.m's own approach
% (which is how s600_17_to_7.mat/s600_region_names.mat were built, so the
% region names have to come from the SAME atlas selection method to match).

s17_names_to_ignore = {'Background+FreeSurfer_Defined_Medial_Wall L', ...
 'Background+FreeSurfer_Defined_Medial_Wall R'   };

%These are AVERAGES for each vertex in the parcel.
sub_rgn_vol_dicts = containers.Map(); %thickness x support area
sub_rgn_res_dicts = containers.Map(); %L/A, assuming uniform resistivity across all regions
sub_rgn_thick_dicts = containers.Map(); %mean cortical thickness
sub_rgn_sa_dicts = containers.Map();    %mean dipole support area (note: this is not the same as the parcel area)
for s=1:numel(SubjectNames)
    subject=SubjectNames{s}; 
    sub_anat_filename = sprintf(['/export02/data/vikramn/OmegaSubset/anat/' ...
        '%s/tess_cortex_central_low.mat'], subject);
    support_areas = subject_support_areas(subject);
    thickness = subject_thicknesses(subject);

    % Per-vertex L/A resistance assuming uniform gray-matter resistivity --
    % saved per subject so cell_res_avg.m can load it instead of
    % recomputing the same formula for its own uniform-resistivity fallback.
    vertex_res_map_gm = avg_gm_res * thickness ./ support_areas;
    save(fullfile(output_dir, sprintf('%s_vertex_res_gm.mat', subject)), 'vertex_res_map_gm');

    surface_file = load(sub_anat_filename);
    atlas_names = {surface_file.Atlas.Name};
    schaeffer_index = find(contains(lower(atlas_names), 'schaefer') & contains(atlas_names, '600') & ...
        contains(lower(atlas_names), '17'), 1);
    if isempty(schaeffer_index)
        error('%s: could not find a Schaefer-600 17-network atlas (has: %s).', subject, strjoin(atlas_names, ', '));
    end
    region_vol_dict = containers.Map();
    region_res_dict = containers.Map();
    region_thick_dict = containers.Map();
    region_sa_dict = containers.Map();
    for region=surface_file.Atlas(schaeffer_index).Scouts %change atlas in brainstorm
        %IGNORE MEDIAL WALL

        region_name_s17=region.Label;
        if ismember(region_name_s17, s17_names_to_ignore)
            continue;
        end

        % The reference subject used to build s600_17_to_7.mat
        % (generate_schaefer_17_to_7_mapping.m: first subject found with
        % both Schaefer-600 atlases) is missing one region
        % ('SalVentAttnB_OFC_1 L') that every other subject's own atlas
        % does have -- a known edge case where atlas-based parcellation
        % assigns zero vertices to a small template region for a
        % particular subject's anatomy, not a data error. Skip it the same
        % way the medial wall is skipped above, rather than crash on a
        % missing map key.
        if ~isKey(s600_17_to_7, region_name_s17)
            continue;
        end
        region_name = s600_17_to_7(region_name_s17);
        vertices_in_rgn = region.Vertices;
        num_v = length(vertices_in_rgn);
        region_res_dict(region_name) = avg_gm_res*sum(thickness(vertices_in_rgn) ./ ...
            support_areas(vertices_in_rgn))/num_v;
        %normalize by parcel area before correlating with power/current?
        region_vol_dict(region_name) = sum(support_areas(vertices_in_rgn) .* ...
            thickness(vertices_in_rgn))/num_v;
        region_thick_dict(region_name) = sum(thickness(vertices_in_rgn))/num_v;
        region_sa_dict(region_name) = sum(support_areas(vertices_in_rgn))/num_v;
        
    end
    sub_rgn_vol_dicts(subject) = region_vol_dict;
    sub_rgn_res_dicts(subject) = region_res_dict;
    sub_rgn_thick_dicts(subject) = region_thick_dict;
    sub_rgn_sa_dicts(subject) = region_sa_dict;
end

%Average across subjects

% region_names (built by generate_schaefer_17_to_7_mapping.m as
% unique(values(s600_17_to_7))) includes the medial wall's 7-net name,
% since that script doesn't apply s17_names_to_ignore when building it --
% but every subject's per-subject loop above deliberately skips the medial
% wall, so no subject's dict ever has this key. Exclude it here the same
% way, translated to the 7-net name space.
ignored_region_names = unique(cellfun(@(k) s600_17_to_7(k), s17_names_to_ignore, 'UniformOutput', false));

avg_rgn_vol_dict = containers.Map();
avg_rgn_res_dict = containers.Map();
avg_rgn_thick_dict = containers.Map();
avg_rgn_sa_dict = containers.Map();

for region_list=region_names
    region_name=region_list{1};
    if ismember(region_name, ignored_region_names)
        continue;
    end
    rgn_total_vol = 0;
    rgn_total_res = 0;
    rgn_total_thick=0;
    rgn_total_sa=0;
    n_subs_with_rgn = 0;

    for s=1:numel(SubjectNames)
        subject = SubjectNames{s};
        sub_vol_dict = sub_rgn_vol_dicts(subject);
        % A small number of subject/region combinations genuinely have no
        % vertices assigned (same anatomical-variation edge case as
        % 'SalVentAttnB_OFC_1 L' noted above, just discovered per-region
        % here instead of universally) -- average over whichever subjects
        % actually have the region rather than crash on a missing key.
        if ~isKey(sub_vol_dict, region_name)
            continue;
        end
        n_subs_with_rgn = n_subs_with_rgn + 1;
        rgn_vol = sub_vol_dict(region_name);
        rgn_total_vol = rgn_total_vol + rgn_vol;

        sub_res_dict = sub_rgn_res_dicts(subject);
        rgn_res = sub_res_dict(region_name);
        rgn_total_res = rgn_total_res + rgn_res;

        sub_thick_dict = sub_rgn_thick_dicts(subject);
        rgn_thick = sub_thick_dict(region_name);
        rgn_total_thick = rgn_total_thick + rgn_thick;

        sub_sa_dict = sub_rgn_sa_dicts(subject);
        rgn_sa = sub_sa_dict(region_name);
        rgn_total_sa = rgn_total_sa + rgn_sa;
    end

    sub_avg_vol = rgn_total_vol/n_subs_with_rgn;
    sub_avg_res = rgn_total_res/n_subs_with_rgn;
    sub_avg_thick = rgn_total_thick/n_subs_with_rgn;
    sub_avg_sa = rgn_total_sa/n_subs_with_rgn;

    avg_rgn_vol_dict(region_name) = sub_avg_vol;
    avg_rgn_res_dict(region_name) = sub_avg_res;
    avg_rgn_thick_dict(region_name) = sub_avg_thick;
    avg_rgn_sa_dict(region_name) = sub_avg_sa;
end

save(fullfile(OUTPUTS_DIR, sprintf('%s_avg_rgn_vol_res_dicts.mat',atlas)),'avg_rgn_vol_dict','avg_rgn_res_dict', ...
    'avg_rgn_thick_dict', 'avg_rgn_sa_dict');
produceOutputCSV(avg_rgn_vol_dict, fullfile(OUTPUTS_DIR, sprintf('%s_subavg_column_vol.csv', atlas)));
produceOutputCSV(avg_rgn_res_dict, fullfile(OUTPUTS_DIR, sprintf('%s_subavg_column_res.csv', atlas)));
produceOutputCSV(avg_rgn_thick_dict, fullfile(OUTPUTS_DIR, sprintf('%s_subavg_column_thick.csv', atlas)));
produceOutputCSV(avg_rgn_sa_dict, fullfile(OUTPUTS_DIR, sprintf('%s_subavg_column_sa.csv', atlas)));