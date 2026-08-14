%Epochs each subject's preprocessed (notch+high-pass filtered) resting-state
%recording into fixed 1-second blocks, registered in the SAME condition
%folder as the continuous recording (e.g. data_block001.mat, data_block002.mat,
%...), for calculate_parcellate_primary_p/rms_idip_rest.m, which reads these
%per-second blocks to compute the RMS of each vertex's primary current dipole
%moment over the first 120 s of resting state.
%
%Uses Brainstorm's process_import_data_time ("Import MEG/EEG: Time") with
%split=1s -- confirmed by inspecting the original hand-computed analysis's
%own data_block*.mat files (TutorialOmega/data/sub-0002/..., History field:
%import_time=[0.000000, 0.999583], i.e. exactly a 1000ms window per block,
%with SSP projectors applied) that this is the process/parameters used.
%
%NOTE: the original hand-computed analysis's block files ended up named
%"..._02.mat" (e.g. data_block001_02.mat) purely because that import was
%run twice by hand in the GUI (an earlier partial attempt left a handful of
%unsuffixed files behind, so Brainstorm auto-suffixed the real, complete
%run "_02" to avoid overwriting them). Since this script runs cleanly once
%per subject with no prior partial attempt, it produces plain
%"data_blockNNN.mat" filenames (no "_02"). Any downstream script that
%hardcodes the "_02" suffix (e.g. calculate_parcellate_primary_p/
%rms_idip_rest.m, calculate_dipole_p_inputs.m's disabled/commented block)
%needs to drop it to read this script's output.
%
%Prerequisites:
% 1) preprocess_meg_notch_highpass.m must have already run for a subject
%    (this script looks up that subject's filtered resting condition from
%    preproc_condition_manifest.mat and skips subjects not yet preprocessed).
%
%NOT YET RUN/VERIFIED against a live Brainstorm session -- spot-check on one
%subject before batch-running across all of SubjectNames.

bst_headless_init();

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

sFiles = [];
bst_report('Start', sFiles);

ProtocolInfo = bst_get('ProtocolInfo');

manifest = load(fullfile(OUTPUTS_DIR, 'preproc_condition_manifest.mat'));
RestConditionManifest = manifest.RestConditionManifest;

for s = 1:numel(SubjectNames)
    subject = SubjectNames{s};

    if ~isKey(RestConditionManifest, subject)
        fprintf('Skipping %s: not preprocessed yet (run preprocess_meg_notch_highpass.m first).\n', subject);
        continue;
    end
    condition = RestConditionManifest(subject);

    % process_import_data_time materializes epochs from a continuous @raw
    % link into a NEW SIBLING condition folder with the same name minus the
    % "@raw" prefix (confirmed both here and in the original hand-computed
    % analysis's own files: TutorialOmega has
    % sub-0002_ses-01_task-rest_run-01_meg_notch_high/data_block*.mat
    % alongside, not inside, @rawsub-0002_..._meg_notch_high/) -- it does
    % NOT write into the same folder as the source continuous link.
    blockCondition = regexprep(condition, '^@raw', '');

    existingBlocks = dir(fullfile(ProtocolInfo.STUDIES, subject, blockCondition, 'data_block*.mat'));
    if ~isempty(existingBlocks)
        fprintf('Skipping %s: already epoched (%d block files in %s).\n', subject, numel(existingBlocks), blockCondition);
        continue;
    end

    sFilesRest = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname', subject, ...
        'condition',   condition);
    if isempty(sFilesRest)
        fprintf('Skipping %s: no resting-state recordings found in %s.\n', subject, condition);
        continue;
    end

    fprintf('=== %s: epoching %s into 1 s blocks ===\n', subject, condition);
    try
        bst_process('CallProcess', 'process_import_data_time', sFilesRest, [], ...
            'subjectname',   subject, ...
            'condition',     '', ...     % empty = let Brainstorm derive the sibling condition name
            'timewindow',    [], ...     % whole recording
            'split',         {1, 's', []}, ...  % 1-second blocks
            'ignoreshort',   1, ...      % drop a trailing partial block shorter than 1 s
            'usectfcomp',    1, ...
            'usessp',        1, ...      % apply the ECG SSP projector computed during preprocessing
            'freq',          [], ...
            'baseline',      [], ...
            'blsensortypes', 'MEG, EEG');
    catch err
        fprintf('  FAILED to epoch %s: %s\n', subject, err.message);
        continue;
    end

    % process_import_data_time is an "import"-type process that reports
    % failures through Brainstorm's report system rather than a MATLAB
    % exception, so verify the real output exists rather than trusting the
    % call above completed silently.
    newBlocks = dir(fullfile(ProtocolInfo.STUDIES, subject, blockCondition, 'data_block*.mat'));
    if isempty(newBlocks)
        fprintf('  FAILED %s: no data_block*.mat files on disk afterward (check the exported report for the real cause).\n', subject);
        continue;
    end
    fprintf('  %s: %d one-second blocks created in %s.\n', subject, numel(newBlocks), blockCondition);
end

db_save();

ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, '/export02/data/vikramn/brainstorm3/reports/');
