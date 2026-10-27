%Calculates "primary dipole P": the power dissipated by each vertex's
%resting-state RMS current dipole moment (rms_idip_rest.m) flowing through
%its cortical column resistance R = rho*L/A (rho=3.5 Ohm.m, Fig. 1A),
%R computed by process_brainstorm_data/parcellate_cortical_column_res_vol.m
%from CAT12 cortical thickness (L) and surface support area (A). Computed
%at every SnrFixed value swept by rms_idip_rest.m.
%
%Renamed from compare_fem_column_power.m: this used to also correlate
%primary P against a *secondary*-current FEM power estimate (a separate
%DUNEuro/Python volume-conductor pipeline, only available for the original
%n=4 subject set -- see meg_power/secondary_p). That comparison now lives in
%compare_pri_sec_dip_p.m and is NOT part of this pipeline.
%
%Idempotent: skips any subject/SnrFixed combination whose output .mat
%already exists, so a restart resumes instead of redoing work.
%
%Prerequisites (produced by process_brainstorm_data):
% - thickness_areas_omega_fem.mat (calculate_dipole_p_inputs.m)
% - <subject>_vertex_res_gm.mat   (parcellate_cortical_column_res_vol.m)
%Prerequisite (produced by rms_idip_rest.m in this pipeline):
% - <subject>_snr<N>_idip_rms.csv

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

% rms_idip_rest.m has now finished SnrFixed=1/3/5 for all 31 subjects, so
% the full sensitivity-analysis sweep is back in scope (was narrowed to [3]
% only while SnrFixed=3 -- the value the main Fig. 1 pipeline and
% metabolism_correlations consume -- was still catching up for everyone).
SnrValues = [1, 3, 5];

% nrow is NOT fixed at 15002 for every subject: each subject's own low-res
% cortex surface has a slightly different actual vertex count (a normal
% side effect of Brainstorm's surface downsampling) -- derived per subject
% below from the loaded RMS CSV itself, guaranteed to match.
ncol = 3;

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
RMS_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'idip_rms');       % rms_idip_rest.m output
RES_DIR = fullfile(OUTPUTS_DIR, 'resistances');                  % process_brainstorm_data/parcellate_cortical_column_res_vol.m output
POWER_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'fem_column_power');
if ~exist(POWER_DIR, 'dir')
    mkdir(POWER_DIR);
end

thicknessAreasFile = fullfile(OUTPUTS_DIR, 'thickness_areas_omega_fem.mat');
if ~exist(thicknessAreasFile, 'file')
    error('thickness_areas_omega_fem.mat not found in %s -- run process_brainstorm_data/calculate_dipole_p_inputs.m first.', OUTPUTS_DIR);
end
thicknessAreas = load(thicknessAreasFile);

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    resFile = fullfile(RES_DIR, sprintf('%s_vertex_res_gm.mat', subject));
    if ~exist(resFile, 'file')
        fprintf('Skipping %s: no cortical column resistance file yet (run process_brainstorm_data/parcellate_cortical_column_res_vol.m first).\n', subject);
        continue;
    end
    resData = load(resFile, 'vertex_res_map_gm');
    R = resData.vertex_res_map_gm; % per-vertex R = rho*L/A (Ohms)

    thickness_in_m = thicknessAreas.subject_thicknesses(subject);

    for j = 1:numel(SnrValues)
        snr = SnrValues(j);
        outFile = fullfile(POWER_DIR, sprintf('%s_snr%d_dip_p_i_rms.mat', subject, snr));
        if exist(outFile, 'file')
            fprintf('Skipping %s SnrFixed=%d: already computed (%s).\n', subject, snr, outFile);
            continue;
        end

        rmsFile = fullfile(RMS_DIR, sprintf('%s_snr%d_idip_rms.csv', subject, snr));
        if ~exist(rmsFile, 'file')
            fprintf('Skipping %s SnrFixed=%d: no RMS I_dip yet (run rms_idip_rest.m first).\n', subject, snr);
            continue;
        end
        rms_idip = load(rmsFile);
        nrow = size(rms_idip, 1); % this subject's actual vertex count
        curr_xyz = reshape(rms_idip, [nrow, 1, ncol]);

        power_vs_t_all = calculatePower(curr_xyz, thickness_in_m, R);
        dip_p = power_vs_t_all{2};       % per-vertex power, in W
        dip_i = vecnorm(curr_xyz, 2, 3); % per-vertex RMS current dipole moment, in A.m (Fig. 1B(i))

        save(outFile, 'dip_p', 'dip_i');
        fprintf('  Saved %s (total P = %s W)\n', outFile, num2str(sum(dip_p)));
    end
end
