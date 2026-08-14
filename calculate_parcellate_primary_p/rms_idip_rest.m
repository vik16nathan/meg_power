%GOAL: Calculate the RMS of each vertex's 3D primary current dipole moment
%across the first 120 s of resting-state recording (Fig. 1A), at every
%SnrFixed value swept by process_brainstorm_data/process_meg_compute_mn_source_fem.m
%(sensitivity analysis on the effect of the MNE regularization parameter on
%the reconstructed primary current dipoles).
%
%Mirrors process_brainstorm_data's headless conventions: same 32-subject
%list, OUTPUTS_DIR-relative I/O so this can run unattended via
%run_primary_p_pipeline.sh. Idempotent: skips any subject/SnrFixed
%combination whose output CSV already exists, so a restart resumes instead
%of redoing work. SnrValues is ordered [3, 1, 5] (not sweep order) so every
%subject's SnrFixed=3 -- the value the main Fig. 1 pipeline and
%metabolism_correlations actually consume -- finishes before the
%sensitivity-analysis-only SnrFixed=1/5 values are attempted.
%
%Computes each block's source estimate directly (ImagingKernel * F, loading
%the kernel ONCE per subject/SnrFixed instead of re-resolving it through
%Brainstorm on every one of the 120 per-second blocks) instead of calling
%in_bst_results('link|kernel|data', 1) per block. Numerically verified
%IDENTICAL (isequal, max abs diff = 0) against in_bst_results for the same
%block on 2026-08-05. bst_headless_init() is still called once, only to
%resolve ProtocolInfo.STUDIES (the on-disk root the manifests' relative
%paths are under) -- no other Brainstorm calls happen anywhere in this
%script, so it no longer needs to avoid running concurrently with a
%Brainstorm session touching the same protocol.
%
%Prerequisites (produced by process_brainstorm_data):
% - preproc_condition_manifest.mat (preprocess_meg_notch_highpass.m)
% - mn_kernel_manifest.mat         (process_meg_compute_mn_source_fem.m)

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

SnrValues = [3, 1, 5]; %SnrFixed=3 (main pipeline) prioritized ahead of the 1/5 sensitivity-analysis values

num_seconds = 120; %first 120 s of resting-state recording (Fig. 1A)
ncol = 2400;         %samples per 1 s data block (2400 Hz sampling)
% nrow/num_v are NOT fixed at 45006/15002 for every subject: each subject's
% own low-res cortex surface (tess_cortex_central_low.mat) has a slightly
% different actual vertex count (15002-15011 observed across the cohort --
% a normal side effect of Brainstorm's surface downsampling), so nrow (=3x
% vertex count, unconstrained orientation) is derived per subject below
% from the loaded dipole data itself, guaranteed to match.

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

kernelManifestFile = fullfile(OUTPUTS_DIR, 'mn_kernel_manifest.mat');
if ~exist(kernelManifestFile, 'file')
    error('mn_kernel_manifest.mat not found in %s -- run process_brainstorm_data/process_meg_compute_mn_source_fem.m first.', OUTPUTS_DIR);
end
kernelData = load(kernelManifestFile);
KernelManifest = kernelData.KernelManifest;

for j = 1:numel(SnrValues) % outer loop: SnrFixed=3 for every subject before 1/5
    snr = SnrValues(j);

    for s = 1:numel(SubjectNames)
        subject = SubjectNames{s};

        if ~isKey(RestConditionManifest, subject)
            fprintf('Skipping %s: not preprocessed yet (run preprocess_meg_notch_highpass.m first).\n', subject);
            continue;
        end
        condition = RestConditionManifest(subject);

        outFile = fullfile(RMS_DIR, sprintf('%s_snr%d_idip_rms.csv', subject, snr));
        if exist(outFile, 'file')
            fprintf('Skipping %s SnrFixed=%d: already computed (%s).\n', subject, snr, outFile);
            continue;
        end

        kernelKey = sprintf('%s_snr%d', subject, snr);
        if ~isKey(KernelManifest, kernelKey)
            fprintf('Skipping %s SnrFixed=%d: no MN kernel yet (run process_meg_compute_mn_source_fem.m first).\n', subject, snr);
            continue;
        end
        kernel_file = KernelManifest(kernelKey);

        % process_import_data_time (process_brainstorm_data/epoch_meg_rest_1s_blocks.m)
        % materializes 1 s blocks into a SIBLING condition folder with the
        % same name minus the "@raw" prefix, not inside the continuous
        % @raw link's own folder -- and produces plain "data_blockNNN.mat"
        % filenames (no "_02": that suffix in the original hand-computed
        % analysis was just an artifact of running the import twice by
        % hand, not a required naming convention).
        blockCondition = regexprep(condition, '^@raw', '');

        % Load the kernel ONCE for this subject/SnrFixed (previously
        % re-resolved via in_bst_results on every one of the 120 blocks --
        % ImagingKernel alone is ~98 MB, so that was 120x redundant disk
        % reads of the same matrix for no reason).
        kernelData = load(fullfile(STUDIES, kernel_file), 'ImagingKernel', 'GoodChannel');

        fprintf('=== %s SnrFixed=%d: RMS I_dip across %d s ===\n', subject, snr, num_seconds);
        sum_sq_idip = [];
        for t = 1:num_seconds
            trial_num = sprintf('%03d', t);
            data_file = sprintf('%s/%s/data_block%s.mat', subject, blockCondition, trial_num);
            blockData = load(fullfile(STUDIES, data_file), 'F');
            amp = kernelData.ImagingKernel * blockData.F(kernelData.GoodChannel, :);

            nrow = size(amp, 1); % 3 x this subject's actual vertex count (unconstrained orientation)
            if isempty(sum_sq_idip)
                sum_sq_idip = zeros(nrow/3, 3);
            end
            curr_xyz = convertTo3DArray(amp, nrow, ncol);
            ss_curr = squeeze(sum(curr_xyz.^2, 2));
            sum_sq_idip = sum_sq_idip + ss_curr;

            % Per-iteration progress so run_primary_p_pipeline.sh's
            % stall-detector (kills MATLAB if its log stops growing for 10
            % min) doesn't mistake slow-but-working computation for a hang.
            if mod(t, 10) == 0 || t == num_seconds
                fprintf('  ...%d/%d s\n', t, num_seconds);
            end
        end
        ms_idip = sum_sq_idip ./ (ncol * num_seconds);
        rms_idip = sqrt(ms_idip);
        csvwrite(outFile, rms_idip);
        fprintf('  Saved %s\n', outFile);
    end
end
