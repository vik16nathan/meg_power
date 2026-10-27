%Companion to process_meg_compute_mn_source_fem.m: computes the SAME MN
%source estimation but with CONSTRAINED (fixed, normal-to-cortex) dipole
%orientation instead of free/unconstrained -- for a targeted sensitivity
%analysis on whether the orientation constraint affects the estimated I_dip
%RMS (calculate_parcellate_primary_p/rms_idip_rest_constrained.m +
%snr_sensitivity_idip_rms.R). Only SnrFixed=3 is computed here (the main
%pipeline's value) -- this is NOT part of the SnrFixed=[1,3,5] sweep.
%
%Prerequisites: same as process_meg_compute_mn_source_fem.m (that script's
%docstring). Run this AFTER process_meg_compute_mn_source_fem.m has already
%completed for a subject (this script reuses its noise-covariance
%computation, so if that script's manifest entries exist the noisecov
%already exists on disk too).
%
%NOT YET RUN/VERIFIED against a live Brainstorm session -- spot-check on
%one subject before batch-running across all of SubjectNames (see
%TEST_SINGLE_SUBJECT below).

bst_headless_init();
db_reload_database('current'); % see process_meg_compute_mn_source_fem.m for why this is needed

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

% Spot-check switch: restrict to a single subject for an initial background
% test before committing to the full 31-subject batch. Single-subject test
% (sub-0002) verified 2026-08-06 -- kernel computed cleanly with 'fixed'
% surface orientations confirmed in the log. Widened to the full cohort.
TEST_SINGLE_SUBJECT = false;
if TEST_SINGLE_SUBJECT
    SubjectNames = {'sub-0002'};
end

SnrFixed = 3; % KEEP fixed at 3 -- this is a constrained-vs-unconstrained analysis, not a SnrFixed sweep

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';

sFiles = [];
bst_report('Start', sFiles);

ProtocolInfo = bst_get('ProtocolInfo');

manifest = load(fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat'));
RestConditionManifest = manifest.RestConditionManifest;
NoiseConditionManifest = manifest.NoiseConditionManifest;

% Separate manifest from the unconstrained pipeline's mn_kernel_manifest.mat
% -- different orientation, different kernel shape (1x vs 3x vertex count),
% so keeping it in its own file avoids any risk of downstream scripts
% accidentally reading a constrained kernel through the unconstrained key
% scheme (or vice versa).
kernelManifestFile = fullfile(OUTPUTS_DIR, 'mn_kernel_manifest_constrained.mat');
if exist(kernelManifestFile, 'file')
    prevManifest = load(kernelManifestFile);
    KernelManifestConstrained = prevManifest.KernelManifestConstrained;
else
    KernelManifestConstrained = containers.Map();
end

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    if ~isKey(RestConditionManifest, subject)
        fprintf('Skipping %s: not preprocessed yet (run preprocess_meg_notch_highpass.m first).\n', subject);
        continue;
    end
    condition = RestConditionManifest(subject);

    existingHeadmodel = dir(fullfile(ProtocolInfo.STUDIES, subject, condition, '*duneuro*.mat'));
    if isempty(existingHeadmodel)
        fprintf('Skipping %s: no FEM head model yet (run build_duneuro_forward_models.m first).\n', subject);
        continue;
    end

    if isKey(KernelManifestConstrained, sprintf('%s_snr%d', subject, SnrFixed))
        fprintf('Skipping %s: constrained source estimation already done.\n', subject);
        continue;
    end

    % Noise covariance: identical to process_meg_compute_mn_source_fem.m
    % (same subject, same condition -- orientation constraint doesn't
    % change the noise covariance computation at all). Reuses whatever
    % noisecov*.mat already exists if that script already ran; computes it
    % fresh here otherwise so this script can also run standalone.
    existingNoiseCov = dir(fullfile(ProtocolInfo.STUDIES, subject, condition, 'noisecov*.mat'));
    if isempty(existingNoiseCov)
        if ~isKey(NoiseConditionManifest, subject)
            fprintf('Skipping %s: no same-session noise recording was preprocessed (check preprocess_meg_notch_highpass.m''s output for this subject).\n', subject);
            continue;
        end
        noiseCondition = NoiseConditionManifest(subject);
        sFilesNoise = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname', subject, ...
            'condition',   noiseCondition);
        if isempty(sFilesNoise)
            fprintf('Skipping %s: could not re-select noise recording %s.\n', subject, noiseCondition);
            continue;
        end
        bst_process('CallProcess', 'process_noisecov', sFilesNoise, [], ...
            'baseline',    [], ...
            'sensortypes', 'MEG', ...
            'target',      1, ...
            'dcoffset',    1, ...
            'identity',    0, ...
            'copycond',    1, ...
            'copysubj',    0, ...
            'copymatch',   0, ...
            'replacefile', 1);
    end

    sFilesRest = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname', subject, ...
        'condition',   condition);
    if isempty(sFilesRest)
        fprintf('Skipping %s: no resting-state recordings found in %s.\n', subject, condition);
        continue;
    end

    fprintf('=== %s: MN source estimation (FEM forward model, CONSTRAINED orientation, SnrFixed=%d) ===\n', subject, SnrFixed);
    sSrc = bst_process('CallProcess', 'process_inverse_2018', sFilesRest, [], ...
        'output',  2, ...  %Kernel only: one per file
        'inverse', struct(...
            'Comment',        sprintf('MN: MEG (SnrFixed=%d, Constrained)', SnrFixed), ...
            'InverseMethod',  'minnorm', ...
            'InverseMeasure', 'amplitude', ...
            'SourceOrient',   {{'fixed'}}, ... % CONSTRAINED: normal-to-cortex, nComponents=1 (vs 'free' in the main pipeline)
            'Loose',          0.2, ... % no-op for 'fixed' orientation, kept identical to the unconstrained script for a clean diff
            'UseDepth',       1, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'reg', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'fixed', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       SnrFixed, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'MEG'}}));

    if isempty(sSrc)
        fprintf('  FAILED %s (SnrFixed=%d, Constrained): process_inverse_2018 returned no files (check the exported report for the real cause).\n', ...
            subject, SnrFixed);
        continue;
    end
    KernelManifestConstrained(sprintf('%s_snr%d', subject, SnrFixed)) = sSrc(1).FileName;

    save(kernelManifestFile, 'KernelManifestConstrained');
end

%Save and export report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, '/export02/data/vikramn/brainstorm3/reports/');
