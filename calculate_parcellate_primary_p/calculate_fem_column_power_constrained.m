%Companion to calculate_fem_column_power.m: computes "primary dipole P" for
%the CONSTRAINED (fixed, normal-to-cortex) orientation kernel, at
%SnrFixed=3 only (see process_meg_compute_mn_source_fem_constrained.m /
%rms_idip_rest_constrained.m). Same formula as calculate_fem_column_power.m
%(P = (I_dip/thickness)^2 * R, calculatePower.m), but simplified for 1
%orientation component per vertex instead of 3 -- calculatePower.m
%hardcodes 3 components (vecnorm across x/y/z), so this recomputes the
%(much simpler) single-component version inline rather than reusing that
%function. Pure vectorized per-vertex arithmetic on already-computed
%inputs -- no Brainstorm session, no per-block loop -- so this runs in
%seconds per subject, unlike the kernel/RMS computation steps upstream.
%
%Prerequisites (produced by process_brainstorm_data):
% - thickness_areas_omega_fem.mat (calculate_dipole_p_inputs.m)
% - <subject>_vertex_res_gm.mat   (parcellate_cortical_column_res_vol.m)
%Prerequisite (produced by rms_idip_rest_constrained.m in this pipeline):
% - <subject>_snr3_idip_rms_constrained.csv

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

SNR = 3; % constrained orientation was only ever computed at SnrFixed=3

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
RMS_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'idip_rms');
RES_DIR = fullfile(OUTPUTS_DIR, 'resistances');
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

    outFile = fullfile(POWER_DIR, sprintf('%s_snr%d_dip_p_i_rms_constrained.mat', subject, SNR));
    if exist(outFile, 'file')
        fprintf('Skipping %s: already computed (%s).\n', subject, outFile);
        continue;
    end

    rmsFile = fullfile(RMS_DIR, sprintf('%s_snr%d_idip_rms_constrained.csv', subject, SNR));
    if ~exist(rmsFile, 'file')
        fprintf('Skipping %s: no constrained RMS I_dip yet (run rms_idip_rest_constrained.m first).\n', subject);
        continue;
    end
    rms_idip = load(rmsFile); % [nrow x 1], single (normal-direction) component

    I_equiv = rms_idip ./ thickness_in_m;
    dip_p = (I_equiv.^2) .* R; % per-vertex power, W
    dip_i = abs(rms_idip);     % per-vertex RMS current dipole moment magnitude, A.m

    save(outFile, 'dip_p', 'dip_i');
    fprintf('  Saved %s (total P = %s W)\n', outFile, num2str(sum(dip_p)));
end
