%Parcellates primary dipole P and its underlying RMS current dipole moment
%(calculate_fem_column_power.m) into Schaeffer-600 (7-network) regions, at
%every SnrFixed value swept upstream. Also builds each subject's per-region
%support-area map (used to sanity-check parcel sizes).
%
%Uses calculate_fem_column_power.m's output instead of a secondary-current
%FEM power file, unlike the earlier version of this script.
%
%Idempotent: skips a subject's area-map export if its CSV already exists,
%and skips a subject/SnrFixed combination's power parcellation if its .mat
%already exists, so a restart resumes instead of redoing work.
%
%Prerequisites:
% - s600_17_to_7.mat / s600_region_names.mat (process_brainstorm_data/generate_schaefer_17_to_7_mapping.m)
% - thickness_areas_omega_fem.mat            (process_brainstorm_data/calculate_dipole_p_inputs.m)
% - <subject>_snr<N>_dip_p_i_rms.mat          (calculate_fem_column_power.m, this pipeline)

atlas = 's600';
% schaeffer_index was hardcoded to 6, but that's actually
% "Schaefer_400_17net" (402 scouts), not "Schaefer_600_17net" (602 scouts,
% matching s600_region_names.mat's 602 entries) -- confirmed by listing
% each subject's registered atlases directly (process_brainstorm_data hit
% the identical bug in parcellate_cortical_column_res_vol.m). Atlas
% registration order isn't even stable across subjects, so a fixed
% position is fragile regardless of which number is "correct" today.
% Found by NAME per subject below instead, matching
% generate_schaefer_17_to_7_mapping.m's own approach (which is how
% s600_17_to_7.mat/s600_region_names.mat were built, so the region names
% have to come from the SAME atlas selection method to match).

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

% SnrFixed=1/5 are now populated too -- see calculate_fem_column_power.m's
% matching comment.
SnrValues = [1, 3, 5];

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
POWER_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'fem_column_power'); % calculate_fem_column_power.m output
output_dir = fullfile(OUTPUTS_DIR, 'primary_p', 'cent_surf_fem_fwd/');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

s600MappingFile = fullfile(OUTPUTS_DIR, sprintf('%s_17_to_7.mat', atlas));
regionNamesFile = fullfile(OUTPUTS_DIR, sprintf('%s_region_names.mat', atlas));
if ~exist(s600MappingFile, 'file') || ~exist(regionNamesFile, 'file')
    error('%s_17_to_7.mat / %s_region_names.mat not found in %s -- run process_brainstorm_data/generate_schaefer_17_to_7_mapping.m first.', ...
        atlas, atlas, OUTPUTS_DIR);
end
load(s600MappingFile);  % s600_17_to_7
load(regionNamesFile);  % region_names

thicknessAreasFile = fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat');
if ~exist(thicknessAreasFile, 'file')
    error('thickness_areas_omega_fem.mat not found in %s -- run process_brainstorm_data/calculate_dipole_p_inputs.m first.', OUTPUTS_DIR);
end
thicknessAreas = load(thicknessAreasFile);

s17_names_to_ignore = {'Background+FreeSurfer_Defined_Medial_Wall L', ...
    'Background+FreeSurfer_Defined_Medial_Wall R'};

% Build each subject's region->vertices map and support-area CSV once
% (SnrFixed-independent), instead of recomputing/reloading the surface file
% for every swept SnrFixed value below.
region_v_dicts_by_subject = containers.Map();
for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};
    sub_anat_filename = sprintf('/export02/data/vikramn/OmegaSubset/anat/%s/tess_cortex_central_low.mat', subject); %NOT pial
    if ~exist(sub_anat_filename, 'file')
        fprintf('Skipping %s: no anatomy file yet.\n', subject);
        continue;
    end
    support_areas = thicknessAreas.subject_support_areas(subject);
    surface_file = load(sub_anat_filename);
    atlas_names = {surface_file.Atlas.Name};
    schaeffer_index = find(contains(lower(atlas_names), 'schaefer') & contains(atlas_names, '600') & ...
        contains(lower(atlas_names), '17'), 1);
    if isempty(schaeffer_index)
        error('%s: could not find a Schaefer-600 17-network atlas (has: %s).', subject, strjoin(atlas_names, ', '));
    end

    region_v_dict = containers.Map();
    region_area_dict = containers.Map();
    for region = surface_file.Atlas(schaeffer_index).Scouts %change atlas in brainstorm
        region_name_s17 = region.Label;
        if ismember(region_name_s17, s17_names_to_ignore) %IGNORE MEDIAL WALL
            continue;
        end
        % The reference subject used to build s600_17_to_7.mat
        % (generate_schaefer_17_to_7_mapping.m: first subject found with
        % both Schaefer-600 atlases) is missing one region
        % ('SalVentAttnB_OFC_1 L') that every other subject's own atlas
        % does have -- a known edge case where atlas-based parcellation
        % assigns zero vertices to a small template region for a
        % particular subject's anatomy, not a data error (matches
        % process_brainstorm_data/parcellate_cortical_column_res_vol.m's
        % identical fix). Skip it the same way the medial wall is skipped
        % above, rather than crash on a missing map key.
        if ~isKey(s600_17_to_7, region_name_s17)
            continue;
        end
        region_name = s600_17_to_7(region_name_s17); %7 network
        region_v_dict(region_name) = region.Vertices;
        region_area_dict(region_name) = sum(support_areas(region.Vertices));
    end
    region_v_dicts_by_subject(subject) = region_v_dict;

    areaCsv = fullfile(output_dir, sprintf('%s_rgn_area_dict.csv', subject));
    if ~exist(areaCsv, 'file')
        produceOutputCSV(region_area_dict, areaCsv);
    end
end

for snrIdx = 1:numel(SnrValues)
    snr = SnrValues(snrIdx);
    fprintf('=== SnrFixed=%d ===\n', snr);

    sub_rgn_p_dicts = containers.Map();
    sub_rgn_i_dicts = containers.Map();

    for s = 1:numel(SubjectNames)
        subject = SubjectNames{s};
        if ~isKey(region_v_dicts_by_subject, subject)
            continue; % no anatomy file (already reported above)
        end

        powerFile = fullfile(POWER_DIR, sprintf('%s_snr%d_dip_p_i_rms.mat', subject, snr));
        if ~exist(powerFile, 'file')
            fprintf('  Skipping %s: no primary P yet (run calculate_fem_column_power.m first).\n', subject);
            continue;
        end
        powerData = load(powerFile, 'dip_p', 'dip_i');
        dipole_p_v_t = powerData.dip_p;
        dipole_i_v_t = powerData.dip_i;

        region_v_dict = region_v_dicts_by_subject(subject);
        region_pvt_dict = containers.Map();
        region_ivt_dict = containers.Map();
        fprintf('  Subject: %s\n', subject);
        for region_list = region_names
            region_name = region_list{1};
            % A handful of 7-network region names have no vertices in this
            % subject's own atlas (same anatomical-variation edge case as
            % 'SalVentAttnB_OFC_1 L' above, i.e. region_v_dict just never
            % got an entry for them) -- skip rather than crash on a missing
            % key (matches process_brainstorm_data/parcellate_cortical_column_res_vol.m's
            % identical guard).
            if ~isKey(region_v_dict, region_name)
                continue;
            end
            vertices_in_rgn = region_v_dict(region_name);
            region_pvt_dict(region_name) = sum(dipole_p_v_t(vertices_in_rgn, :), 1); %NOTE - average power per trial
            region_ivt_dict(region_name) = sum(dipole_i_v_t(vertices_in_rgn, :), 1);
        end

        save(fullfile(output_dir, sprintf('%s_%s_snr%d_region_pri_p_i.mat', atlas, subject, snr)), ...
            'region_pvt_dict', 'region_ivt_dict');

        produceOutputCSV(region_ivt_dict, fullfile(output_dir, sprintf('%s_%s_snr%d_i_pri_2min_rest.csv', atlas, subject, snr)));
        produceOutputCSV(region_pvt_dict, fullfile(output_dir, sprintf('%s_%s_snr%d_p_pri_2min_rest.csv', atlas, subject, snr)));

        sub_rgn_p_dicts(subject) = region_pvt_dict;
        sub_rgn_i_dicts(subject) = region_ivt_dict;
    end

    if sub_rgn_p_dicts.Count == 0
        fprintf('  No subjects with primary P at SnrFixed=%d yet -- skipping subject average.\n', snr);
        continue;
    end

    %Subject average - to compare with PET averaged across individuals
    %(note that sample sizes are different)
    avg_rgn_p_dict = containers.Map();
    avg_rgn_i_dict = containers.Map();
    doneSubjects = keys(sub_rgn_p_dicts);
    for region_list = region_names
        region_name = region_list{1};
        rgn_total_power_vec = 0;
        rgn_total_current_vec = 0;
        n_subs_with_rgn = 0;

        for s = 1:numel(doneSubjects)
            sub_p_dict = sub_rgn_p_dicts(doneSubjects{s});
            % Same missing-region edge case as above -- average over
            % whichever subjects actually have the region rather than
            % crash on a missing key.
            if ~isKey(sub_p_dict, region_name)
                continue;
            end
            n_subs_with_rgn = n_subs_with_rgn + 1;
            rgn_total_power_vec = rgn_total_power_vec + sub_p_dict(region_name);

            sub_i_dict = sub_rgn_i_dicts(doneSubjects{s});
            rgn_total_current_vec = rgn_total_current_vec + sub_i_dict(region_name);
        end

        if n_subs_with_rgn == 0
            continue;
        end
        avg_rgn_p_dict(region_name) = rgn_total_power_vec / n_subs_with_rgn;
        avg_rgn_i_dict(region_name) = rgn_total_current_vec / n_subs_with_rgn;
    end

    produceOutputCSV(avg_rgn_p_dict, fullfile(output_dir, sprintf('%s_p_pri_2min_rest_subavg_snr%d.csv', atlas, snr)));
    produceOutputCSV(avg_rgn_i_dict, fullfile(output_dir, sprintf('%s_i_pri_2min_rest_subavg_snr%d.csv', atlas, snr)));
end
