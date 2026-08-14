%Automates MEG minimum-norm (MN) source estimation with a FEM forward model,
%sweeping the fixed-SNR regularization parameter used by process_inverse_2018,
%instead of clicking through each subject/parameter combination manually in
%Brainstorm. Mirrors the "SOURCE ESTIMATION" section of
%brainstorm3/toolbox/script/tutorial_omega.m (the same script Brainstorm's
%own OMEGA tutorial is generated from), but replaces the overlapping-spheres
%MEG forward model with DUNEuro FEM (computed separately, see prerequisites).
%
%Prerequisites:
% 1) preprocess_meg_notch_highpass.m must have already run for a subject
%    (this script looks up that subject's filtered resting condition AND
%    same-session filtered noise condition from preproc_condition_manifest.mat,
%    and skips subjects not yet preprocessed).
% 2) build_duneuro_forward_models.m must have already computed that
%    subject's DUNEuro FEM head model (auto-detected below; subjects
%    without one yet are skipped).
%
%NOT YET RUN/VERIFIED against a live Brainstorm session -- spot-check on one
%subject before batch-running across all of SubjectNames.

bst_headless_init();

% build_duneuro_forward_models.m runs in a separate, earlier MATLAB session,
% and this fresh session doesn't reliably see the head models it wrote
% without an explicit reload -- process_inverse_2018 fails with "No head
% model available" otherwise, even though the file is really on disk (same
% class of staleness issue as the FEM-mesh one fixed in
% build_duneuro_forward_models.m, just at the study/head-model level here).
db_reload_database('current');

SubjectNames = { ...
    'sub-0002', 'sub-0008', 'sub-0009', 'sub-0010', ...
    'sub-0011', 'sub-0012', 'sub-0014', 'sub-0015', ...
    'sub-0016', 'sub-0018', 'sub-0019', 'sub-0020', ...
    'sub-0021', 'sub-0022', 'sub-0023', 'sub-0024', 'sub-0025', ...
    'sub-0026', 'sub-0028', 'sub-0029', 'sub-0030', 'sub-0031', ...
    'sub-0032', 'sub-0033', 'sub-0034', 'sub-0035', 'sub-0036', ...
    'sub-0037', 'sub-0039', 'sub-0040', 'sub-0041', ...
};

SnrValues = [1, 3, 5]; %sensitivity analysis on the MNE regularization parameter

OUTPUTS_DIR = '/export02/data/vikramn/hbm_manuscript_code/outputs/';

sFiles = [];
bst_report('Start', sFiles);

ProtocolInfo = bst_get('ProtocolInfo');

%RestConditionManifest/NoiseConditionManifest map subject -> the
%preprocessed (notch+high-pass) resting-state / same-session task-noise
%condition names, produced by preprocess_meg_notch_highpass.m. Both come
%from the SAME session for every subject, so noise covariance can be
%computed directly from each subject's own data with no cross-session or
%cross-subject date matching.
manifest = load(fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat'));
RestConditionManifest = manifest.RestConditionManifest;
NoiseConditionManifest = manifest.NoiseConditionManifest;

%KernelManifest maps '<subject>_snr<value>' -> the MN_MEG_KERNEL result file
%Brainstorm generates for that subject/SnrFixed combination, so downstream
%scripts (calculate_dipole_p_inputs.m) don't need to guess result-file timestamps.
%Reload any existing manifest so re-running this script doesn't lose entries
%for subjects processed in an earlier run.
kernelManifestFile = fullfile(OUTPUTS_DIR, 'mn_kernel_manifest.mat');
if exist(kernelManifestFile, 'file')
    prevManifest = load(kernelManifestFile);
    KernelManifest = prevManifest.KernelManifest;
else
    KernelManifest = containers.Map();
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

    if isKey(KernelManifest, sprintf('%s_snr%d', subject, SnrValues(end)))
        fprintf('Skipping %s: source estimation already done.\n', subject);
        continue;
    end

    % Compute noise covariance from this subject's own SAME-SESSION
    % filtered noise recording (preprocess_meg_notch_highpass.m only ever
    % picks a session that has both rest and noise data together, so this
    % is always available -- no cross-session/cross-subject date matching
    % needed). 'copycond',1 copies to all of this subject's folders
    % (reaches the rest condition regardless of exact folder name);
    % 'copysubj',0 keeps it within this subject; 'copymatch',0 since there's
    % no ambiguity to resolve (avoids a real bug in Brainstorm's
    % process_noisecov.m where combining copycond=1/copysubj=0/copymatch=1
    % crashes with "Unrecognized function or variable 'ProtocolStudies'",
    % since that variable is only initialized in the copycond&&copysubj
    % branch but referenced unconditionally in the isCopyMatch block).
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

    fprintf('=== %s: MN source estimation (FEM forward model, SnrFixed sweep) ===\n', subject);
    for j = 1:numel(SnrValues)
        snr = SnrValues(j);
        sSrc = bst_process('CallProcess', 'process_inverse_2018', sFilesRest, [], ...
            'output',  2, ...  %Kernel only: one per file
            'inverse', struct(...
                'Comment',        sprintf('MN: MEG (SnrFixed=%d)', snr), ...
                'InverseMethod',  'minnorm', ...
                'InverseMeasure', 'amplitude', ...
                'SourceOrient',   {{'free'}}, ... % matches the hand-computed original (Comment 'FEM(Unconstr)', nComponents=3)
                'Loose',          0.2, ...
                'UseDepth',       1, ...
                'WeightExp',      0.5, ...
                'WeightLimit',    10, ...
                'NoiseMethod',    'reg', ...
                'NoiseReg',       0.1, ...
                'SnrMethod',      'fixed', ...
                'SnrRms',         1e-06, ...
                'SnrFixed',       snr, ...
                'ComputeKernel',  1, ...
                'DataTypes',      {{'MEG'}}));

        % process_inverse_2018 reports failures through Brainstorm's report
        % system rather than a MATLAB exception, so a failure here comes
        % back as an empty result rather than throwing -- check before
        % indexing instead of crashing the whole script on one bad subject/SNR.
        if isempty(sSrc)
            fprintf('  FAILED %s (SnrFixed=%d): process_inverse_2018 returned no files (check the exported report for the real cause).\n', ...
                subject, snr);
            continue;
        end
        KernelManifest(sprintf('%s_snr%d', subject, snr)) = sSrc(1).FileName;
    end

    save(kernelManifestFile, 'KernelManifest');
end

%Save and export report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, '/export02/data/vikramn/brainstorm3/reports/');
