%Companion to rms_idip_rest.m: computes the RMS of each vertex's primary
%current dipole moment across the first 120 s of resting-state recording,
%using the CONSTRAINED (fixed, normal-to-cortex) MN kernel from
%process_brainstorm_data/process_meg_compute_mn_source_fem_constrained.m,
%instead of the main pipeline's unconstrained (free, 3D) kernel.
%
%Targeted sensitivity analysis: does the orientation constraint change the
%estimated I_dip RMS at SnrFixed=3 (see snr_sensitivity_idip_rms.R for the
%comparison against rms_idip_rest.m's unconstrained SnrFixed=3 output)?
%SnrFixed is fixed at 3 -- this is NOT a SnrFixed sweep.
%
%Constrained orientation means 1 component per vertex (a signed scalar:
%positive/negative = current flowing out of/into the cortical surface along
%its normal), not 3 (x,y,z) -- so unlike rms_idip_rest.m, no
%convertTo3DArray de-interleaving or per-axis sum-of-squares is needed here;
%RMS is just sqrt(mean(amp.^2, 2)) per vertex directly.
%
%Prerequisites (produced by process_brainstorm_data):
% - preproc_condition_manifest.mat            (preprocess_meg_notch_highpass.m)
% - mn_kernel_manifest_constrained.mat         (process_meg_compute_mn_source_fem_constrained.m)

addpath('/export02/data/vikramn/hbm_manuscript_code/process_brainstorm_data'); % for bst_headless_init()
bst_headless_init();
ProtocolInfo = bst_get('ProtocolInfo');
STUDIES = ProtocolInfo.STUDIES;

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

% Mirrors process_meg_compute_mn_source_fem_constrained.m's spot-check
% switch. Single-subject test (sub-0002) verified 2026-08-06. Widened to
% the full cohort.
TEST_SINGLE_SUBJECT = false;
if TEST_SINGLE_SUBJECT
    SubjectNames = {'sub-0002'};
end

SnrFixed = 3; % KEEP fixed at 3 -- see module docstring

num_seconds = 120; %first 120 s of resting-state recording, matches rms_idip_rest.m
ncol = 2400;         %samples per 1 s data block (2400 Hz sampling)

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';
RMS_DIR = fullfile(OUTPUTS_DIR, 'primary_p', 'idip_rms');
if ~exist(RMS_DIR, 'dir')
    mkdir(RMS_DIR);
end

manifestFile = fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat');
if ~exist(manifestFile, 'file')
    error('preproc_condition_manifest.mat not found in %s -- run process_brainstorm_data/preprocess_meg_notch_highpass.m first.', OUTPUTS_DIR);
end
manifest = load(manifestFile);
RestConditionManifest = manifest.RestConditionManifest;

kernelManifestFile = fullfile(OUTPUTS_DIR, 'mn_kernel_manifest_constrained.mat');
if ~exist(kernelManifestFile, 'file')
    error('mn_kernel_manifest_constrained.mat not found in %s -- run process_brainstorm_data/process_meg_compute_mn_source_fem_constrained.m first.', OUTPUTS_DIR);
end
kernelData = load(kernelManifestFile);
KernelManifestConstrained = kernelData.KernelManifestConstrained;

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    if ~isKey(RestConditionManifest, subject)
        fprintf('Skipping %s: not preprocessed yet (run preprocess_meg_notch_highpass.m first).\n', subject);
        continue;
    end
    condition = RestConditionManifest(subject);

    outFile = fullfile(RMS_DIR, sprintf('%s_snr%d_idip_rms_constrained.csv', subject, SnrFixed));
    if exist(outFile, 'file')
        fprintf('Skipping %s SnrFixed=%d (Constrained): already computed (%s).\n', subject, SnrFixed, outFile);
        continue;
    end

    kernelKey = sprintf('%s_snr%d', subject, SnrFixed);
    if ~isKey(KernelManifestConstrained, kernelKey)
        fprintf('Skipping %s SnrFixed=%d (Constrained): no MN kernel yet (run process_meg_compute_mn_source_fem_constrained.m first).\n', subject, SnrFixed);
        continue;
    end
    kernel_file = KernelManifestConstrained(kernelKey);

    % Same block-materialization convention as rms_idip_rest.m -- reuses
    % the SAME data_blockNNN.mat files (raw sensor data doesn't depend on
    % source orientation, only the kernel does).
    blockCondition = regexprep(condition, '^@raw', '');

    kernelData = load(fullfile(STUDIES, kernel_file), 'ImagingKernel', 'GoodChannel');

    fprintf('=== %s SnrFixed=%d (Constrained): RMS I_dip across %d s ===\n', subject, SnrFixed, num_seconds);
    sum_sq_idip = [];
    for t = 1:num_seconds
        trial_num = sprintf('%03d', t);
        data_file = sprintf('%s/%s/data_block%s.mat', subject, blockCondition, trial_num);
        blockData = load(fullfile(STUDIES, data_file), 'F');
        amp = kernelData.ImagingKernel * blockData.F(kernelData.GoodChannel, :);

        nrow = size(amp, 1); % = this subject's actual vertex count (constrained orientation: 1 component/vertex)
        if isempty(sum_sq_idip)
            sum_sq_idip = zeros(nrow, 1);
        end
        sum_sq_idip = sum_sq_idip + sum(amp.^2, 2);

        if mod(t, 10) == 0 || t == num_seconds
            fprintf('  ...%d/%d s\n', t, num_seconds);
        end
    end
    ms_idip = sum_sq_idip ./ (ncol * num_seconds);
    rms_idip = sqrt(ms_idip);
    csvwrite(outFile, rms_idip);
    fprintf('  Saved %s\n', outFile);
end
